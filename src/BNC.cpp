/*
 * @file bnc.cpp
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 10/17/2014
 */

#include "../include/BNC.h"

ILOSTLBEGIN

BNC::BNC(Instance& instance) {
    this->i = new Instance("");
    this->i->operator=(instance);
}

BNC::~BNC() {
    delete i;
}

inline double max(double &lhs, double &rhs) {
    if (lhs > rhs)
        return lhs;
    return rhs;
}

inline double sum(double ** sol, std::set<int> Smin) {
    double s = 0;
    std::vector<int> vmin(Smin.begin(), Smin.end());
    for (unsigned i = 0; i < vmin.size(); i++) {
        for (unsigned j = i + 1; j < vmin.size(); j++) {
            s += sol[vmin[i]][vmin[j]];
        }
    }
    return s;
}

void BNC::maxBack(double ** sol, std::set<int>& Smin) {
    int dim = this->i->getDim();
    double *b = new double[dim];
    double cutmin = 0, cutval = 0;
    std::set<int> V, S;
    std::set<int>::iterator it, jt;

    int vertex, i, j, m;

    V.clear();
    S.clear();
    Smin.clear();

    // inicializacao algoritmo max-back
    for (i = 1; i < dim; i++) {
        V.insert(i);
        b[i] = 0.0;
    }
    b[0] = 0.0;

    S.insert(0);
    // soma dos pesos (X*) das arestas cruzando S0
    for (i = 1; i < dim; i++) {
        cutmin += sol[0][i];
        b[i] = sol[0][i];
    }

    cutval = cutmin;
    Smin = S;

    while (S.size() < dim) {
        // seleciona v fora de S de maior max-back
        m = -1.0;
        for (it = V.begin(); it != V.end(); ++it) {
            i = *it;
            if (b[i] > m) {
                m = b[i];
                vertex = i;
            }
        }

        S.insert(vertex); // S = S + v
        V.erase(V.find(vertex));
        cutval = cutval + 2 - (2 * b[vertex]);

        for (it = V.begin(); it != V.end(); ++it) {
            i = *it;
            if (i < vertex)
                b[i] += sol[i][vertex];
            else
                b[i] += sol[vertex][i];
        }

        if (cutval < (cutmin - bnc::eps)) {
            cutmin = cutval;
            Smin = S;
        }
    }

    //    vector<int> t(Smin.begin(), Smin.end());
    //    for(i=0;i<t.size();i++)
    //        cout << t[i] << " ";
    //    cout << endl;

    delete b;
}

int BNC::CPXPUBLIC mycutcallback(CPXCENVptr env, void *cbdata, int wherefrom, 
        void *cbhandle, int *useraction_p) {
    using namespace std;

    int i, j, m, dim;
    set<int> Smin;

    // lock
    pthread_mutex_t cs_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&cs_mutex);

    *useraction_p = CPX_CALLBACK_DEFAULT;

    BNC * bnc = (BNC * ) cbhandle;
    dim = bnc->i->getDim();

    int status = 0;

    double *X = new double[bnc->numCols];
    double **sol = new double * [dim];
    for (i = 0; i < dim; i++) {
        sol[i] = new double[dim];
    }

    CPXgetcallbacknodex(env, cbdata, wherefrom, X, 0, bnc->numCols - 1);
    bnc::storeLPSolution(env, bnc->model, bnc->numCols, X, sol);

    // max-back
    bnc->maxBack(sol, Smin);

    if (Smin.size() < dim) {
        vector<int> vmin(Smin.begin(), Smin.end());

        char varName[100];
        int colIndex, counter = 0;
        int cutnz = (Smin.size() * (Smin.size() - 1)) / 2;
        int * cutInd = new int[cutnz];
        double * cutVal = new double[cutnz];
        double left = 0, rhs = Smin.size() - 1;

        for (i = 0; i < vmin.size(); i++) {
            for (j = i + 1; j < vmin.size(); j++) {
                sprintf(varName, "X_%d_%d", vmin[i], vmin[j]);
                CPXgetcolindex(env, bnc->model, varName, &colIndex);
                cutInd[counter] = colIndex;
                cutVal[counter] = 1;
                left += sol[vmin[i]][vmin[j]];
                counter++;
            }
        }

        if (left + bnc::eps <= rhs) {
            cout << "Corte nao violado: " << left << " " << rhs << endl;
        }

        status = CPXcutcallbackadd(env, cbdata, wherefrom, cutnz, rhs, 'L', cutInd, cutVal, 1);

        if (status) {
            cout << "Falha ao adicionar corte." << endl;
        }

        delete[] cutInd;
        delete[] cutVal;

    }

    delete[] X;

    for (i = 0; i < dim; i++)
        delete[] sol[i];
    delete[] sol;

    // unlock
    pthread_mutex_unlock(&cs_mutex);


    return status;
}

void BNC::createLP(const int ** matrix, unsigned dim) {
    int i, j;
    try {
        IloEnv env;
        IloModel model(env);
        IloArray<IloIntVarArray> x(env, dim);

        char var[100];
        for (i = 0; i < dim; i++) {
            IloIntVarArray array(env, dim, 0, 1);
            x[i] = array;
            for (j = i + 1; j < dim; j++) {
                x[i][j].setBounds(0, 1);
                sprintf(var, "X_%d_%d", i, j);
                x[i][j].setName(var);
                model.add(x[i][j]);
            }
        }

        // Restricoes
        char c[100];
        for (i = 0; i < dim; i++) {
            IloExpr sum(env);
            for (j = i + 1; j < dim; j++) {
                sum += x[i][j];

            }

            for (j = 0; j < i; j++) {
                sum += x[j][i];
            }

            IloRange r = (sum == 2);
            sprintf(c, "c_i%d", i);
            r.setName(c);
            model.add(r);
        }

        // Funcao objetivo
        IloExpr cost(env);
        for (i = 0; i < dim; i++) {
            for (j = i + 1; j < dim; j++) {
                cost += matrix[i][j] * x[i][j];
            }
        }
        model.add(IloMinimize(env, cost));

        // Escreve modelo LP para arquivo
        IloCplex TSP(model);
        TSP.exportModel("TSP.lp");
        env.end();

    } catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;   
    }
    catch(...) {
        cerr  << " ERROR" << endl;   
    }
}

int BNC::initBranchAndCut(int ub, std::string instanceName) {
    createLP(this->i->getMatrix(), this->i->getDim());

    // Inicializa ambiente CPLEX
    int status = 0;
    CPXENVptr env = CPXopenCPLEX(&status);
    model = CPXcreateprob(env, &status, "TSP");

    CPXreadcopyprob(env, model, "TSP.lp", NULL);

    numCols = CPXgetnumcols(env, model); 

    // Parametros CPLEX
    CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
    CPXsetintparam(env, CPX_PARAM_MIPINTERVAL, 20);
    CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
    CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
    CPXsetintparam(env, CPX_PARAM_REDUCE, CPX_PREREDUCE_PRIMALONLY);
    CPXsetintparam(env, CPX_PARAM_THREADS, 1);
    CPXsetintparam(env, CPX_PARAM_PARALLELMODE, 1);
    CPXsetintparam(env, CPX_PARAM_VARSEL, 3);
    CPXsetdblparam(env, CPX_PARAM_TILIM, 86400);
    CPXsetdblparam(env, CPX_PARAM_CUTUP, ub);

    // cut callback
    CPXsetcutcallbackfunc(env, mycutcallback, this);

    CPXFILEptr fp = CPXfopen ("TSP.log", "w");

    // Le MIPStart
    //char mst[100];
    //sprintf(mst, "%s.mst", instanceName.c_str());
    //CPXreadcopymipstarts(env, model, mst);

    // Optimize
    double cplexTimeBefore, cplexTimeAfter;	
    CPXgettime(env, &cplexTimeBefore);
    CPXmipopt (env, model);
    CPXgettime(env, &cplexTimeAfter);
    double time = cplexTimeAfter - cplexTimeBefore;

    double objval;
    CPXgetobjval(env, model, &objval);
    ofstream f;
    f.open(i->getName().c_str());
    f << objval;
    f.close();
    

    bnc::printSolution(env, model, numCols, i->getDim());
    CPXfclose(fp);
    // Free
    CPXfreeprob(env, &model);
    CPXcloseCPLEX(&env);

    return 0;
}

