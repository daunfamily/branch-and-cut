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

    double *X = new double[bnc->numCols];
    double *b = new double[dim];

    CPXgetcallbacknodex(env, cbdata, wherefrom, X, 0, bnc->numCols - 1);

    // inicializacao algoritmo max-back
    for (i = 0; i < dim; i++) {
        V.insert(i);
    }

    S.insert(0);
    // soma dos pesos (X*) das arestas cruzando S0
    for (i = 1; i < dim; i++) {
        cutmin += X[i];
        b[i]   += X[i*dim];
    }
    cutval = cutmin;
    S = Smin;

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
        S.insert(*jt); // S = S + [v]

        cutval += 2 - (2*b[*jt]);

        for (it = V.begin(); it != V.end(); ++it) {
            if (S.find(*it) == S.end()) { // v nao pertencente a S
                b[*it] += X[(*jt)*dim + *it];
            }
        }

        if (cutval < cutmin) {
            cutmin = cutval;
            Smin = S;
        }
    }

    // unlock
    pthread_mutex_unlock(&cs_mutex);

    delete[] X;
    delete[] b;

    return 0;
}

void BNC::createLP(const int ** matrix, unsigned dim) {
    int i, j;
    try {
        IloEnv env;
        IloModel model(env);
        IloArray<IloIntVarArray> x(env, dim);

        char var[100];
        for (i = 0; i < dim; i++) {
            IloIntVarArray array(env, dim + 1, 0, 2);
            x[i] = array;
            for (j = 0; j < dim; j++) {
                x[i][j].setBounds(0, 1);
                sprintf(var, "X_%d_%d", i, j);
                x[i][j].setName(var);
                model.add(x[i][j]);

            }
        }

        // Funcao objetivo
        IloExpr cost(env);
        for (i = 0; i < dim; i++) {
            for (j = i + 1; j < dim; j++) {
                cost += matrix[i][j] * x[i][j];
            }
        }
        model.add(IloMinimize(env, cost));

        // Restricoes
        char c[100];
        for (i = 0; i < dim; i++) {
            IloExpr sum(env);
            for (j = 0; j < dim; j++) {
                if (i != j) {
                    sum += x[i][j];
                }
            }

            IloRange r = (sum == 1);
            sprintf(c, "c_i%d", i);
            r.setName(c);
            model.add(r);
        }

        for (j = 0; j < dim; j++) {
            IloExpr sum(env);
            for (i = 0; i < dim; i++) {
                if (i != j) {
                    sum += x[i][j];
                }

            }
            IloRange r = (sum == 1);
            sprintf(c, "c_j%d", i);
            r.setName(c);
            model.add(r);
        }

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
    CPXLPptr  model = CPXcreateprob(env, &status, "TSP");

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

    // Imprime resultados
    //printResults(env, model, instanceName, time);
    //printNumberOfCuts();
    //printResultsToFile(env, model, instanceName, time);
    //printSolution(argv[1], env, model, cur_numcols);

    CPXfclose(fp);

    // Free
    CPXfreeprob(env, &model);
    CPXcloseCPLEX(&env);

    return 0;
}

