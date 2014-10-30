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

int BNC::CPXPUBLIC mycutcallback(CPXCENVptr env, void *cbdata, int wherefrom, 
        void *cbhandle, int *useraction_p) {
    using namespace std;

    int i, j, m, dim, inf = 9999;

    // lock
    pthread_mutex_t cs_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&cs_mutex);

    *useraction_p = CPX_CALLBACK_DEFAULT;

    BNC * bnc = (BNC * ) cbhandle;
    dim = bnc->i->getDim();

    set<int> V, Smin, S;
    set<int>::iterator it, jt;
    double cutmin = 0, cutval = 0;

    int status = 0;

    double *X = new double[bnc->numCols];
    double *b = new double[dim];
    double **sol = new double * [dim];
    for (i = 0; i < dim; i++) {
        sol[i] = new double[dim];
    }

    CPXgetcallbacknodex(env, cbdata, wherefrom, X, 0, bnc->numCols - 1);
    bnc::storeLPSolution(env, bnc->model, bnc->numCols, X, sol);

    // inicializacao algoritmo max-back
    for (i = 0; i < dim; i++) {
        V.insert(i);
        b[i] = 0;
    }

    S.insert(0);
    // soma dos pesos (X*) das arestas cruzando S0
    for (i = 1; i < dim; i++) {
        cutmin += sol[0][i];
        b[i]   += (sol[0][i] > 0) ? sol[0][i] : 0;
    }

    cutval = cutmin;
    Smin = S;

    while(S.size() < dim) {
        // seleciona v fora de S de maior max-back
        m = -1, i = -1;
        for (it = V.begin(); it != V.end(); ++it) {
            if (S.find(*it) == S.end()) { // v nao pertencente a S
                if (b[*it] > m) {
                    m = b[*it];
                    i = *it;
                }
            }
        }

        jt = V.find(i);
        S.insert(*jt); // S = S + v

        cutval += 2 - (2*b[*jt]);

        for (it = V.begin(); it != V.end(); ++it) {
            if (S.find(*it) == S.end()) { // v nao pertencente a S
                if (*jt < *it)
                    b[*it] += sol[*jt][*it];
                else
                    b[*it] += sol[*it][*jt];
            }
        }

        if (cutval < cutmin) {
            cutmin = cutval;
            Smin = S;
        }
    }

    cout << "Smin: (" << Smin.size() << ")" << endl;
    for (it = Smin.begin(); it != Smin.end(); ++it) {
        cout << *it << " ";
    }
    cout << endl;

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
            if (sol[vmin[i]][vmin[j]] > 0)
                left += sol[vmin[i]][vmin[j]];
            counter++;
        }
    }

    cout << "LEFT = " << left << endl;
    if (left <= rhs || Smin.size() == dim) {
        cout << "Corte nao violado: " << left << " " << rhs << endl;
    } else {
        if (CPXcutcallbackadd(env, cbdata, wherefrom, cutnz, rhs, 'L', cutInd, cutVal, 1)) {
            cout << "Falha ao adicionar corte." << endl;
        }
    }

    delete[] cutInd;
    delete[] cutVal;

    delete[] X;
    delete[] b;
    
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
    CPXsetintparam(env, CPX_PARAM_MIPINTERVAL, 100);
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

    CPXFILEptr fp;
    fp = CPXfopen ("TSP.log", "w");

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

    bnc::printSolution(env, model, numCols, i->getDim());
    CPXfclose(fp);

    // Free
    CPXfreeprob(env, &model);
    CPXcloseCPLEX(&env);

    return 0;
}

