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

void adjustSol(double ** sol, int dim) {
    for (int i = 0; i < dim; i++) {
        for (int j = i + 1; j < dim; j++) {
            if (sol[i][j] < bnc::eps)
                sol[i][j] = 0;
        }
    }
}

void BNC::maxBack(double ** sol, std::set<int>& Smin) {
    int vertex, i, j, m, dim = this->i->getDim();
    double *b = new double[dim];
    double cutmin = 0, cutval = 0;
    std::set<int> V, S;
    std::set<int>::iterator it, jt;

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

    adjustSol(sol, dim);

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

    delete b;
}

void BNC::minCut(double ** sol, std::set<int>& Smin, double& mc) {
    using namespace std;
    using namespace bnc;

    int a = 1, c, vertex, i, j, m, dim = this->i->getDim();
    list<list<int>> graph, V, S;
    double ** bkp = sol;
    double cutphasemin = inf, cutphase = 0;

    list<int> l1, l2, l3;
    list<int>::iterator it, jt;
    list<list<int>>::iterator git, gjt, g1, g2;
    list<list<int>>::iterator rgt, lgt;

    adjustSol(sol, dim);
    Smin.clear();
    for (i = 0; i < dim; i++) {
        list<int> l;
        graph.push_back(l);
    }

    i = 0;
    for (git = graph.begin(); git != graph.end(); ++git) {
        git->push_back(i);
        i++;
    }

    // ------------------------------
    // grafo exemplo
    /*     dim = 8;
     *     sol = new double * [dim];
     *     graph.clear();
     *     for (i = 0; i < dim; i++) {
     *         sol[i] = new double [dim];
     *         for (j = 0; j < dim; j++)
     *             sol[i][j] = 0.0;
     * 
     *         list<int> l;
     *         l.push_back(i);
     *         graph.push_back(l);
     *     }
     * 
     *     sol[0][1] = 2;
     *     sol[0][4] = 3;
     *     sol[1][2] = 3;
     *     sol[1][4] = 2;
     *     sol[1][5] = 2;
     *     sol[2][3] = 4;
     *     sol[2][6] = 2;
     *     sol[3][6] = 2;
     *     sol[3][7] = 2;
     *     sol[4][5] = 3;
     *     sol[5][6] = 1;
     *     sol[6][7] = 3;
     * 
     * 
     */
    // -----------------------------

    while (graph.size() > 1) {
        V.clear();
        S.clear();

        V.insert(V.begin(), graph.begin(), graph.end());
        S.push_back(*(findNode(V, a)));
        deleteNode(V, a);

        while (S.size() < graph.size()) {
            // print(S);
            // seleciona v fora de S de maior max-back
            m = -1.0;
            for (gjt = V.begin(); gjt != V.end(); ++gjt) {
                c = 0;
                j = gjt->front();

                for (git = S.begin(); git != S.end(); ++git) {
                    i = git->front();
                    if (i < j) {
                        c += sol[i][j];
                    } else {
                        c += sol[j][i];
                    }
                }

                if (c > m) {
                    m = c, vertex = j;
                }
            }

            S.push_back(*(findNode(V, vertex)));
            deleteNode(V, vertex);
        }

        //        print(S);

        // ultimos dois adicionados 
        l1 = S.back();
        S.pop_back();
        l2 = S.back();
        S.pop_back();

        //cout << "dois ultimos adicionados em S: " << endl;
        //        print(l1);
        //        print(l2);

        cutphase = evalCutPhase(sol, graph, l1);    // valor do cut-of-the-phase

        //        cout << "cut-of-the-phase: " << cutphase << endl;

        updateEdges(sol, graph, l1, l2);            // atualizacao das arestas

        // "encolhimento" do grafo
        l3.clear();                                 // l3 contem novo no
        l3.insert(l3.begin(), l1.begin(), l1.end());
        l3.insert(l3.begin(), l2.begin(), l2.end());

        deleteNode(graph, l1.front());              // remove l1 do grafo
        deleteNode(graph, l2.front());
        graph.push_front(l3);

        //cout << "Graph: " << endl;
        //        print(graph);
        //        cin >> i;

        if (cutphase < cutphasemin) {
            cutphasemin = cutphase;
            Smin.clear();
            Smin.insert(l1.begin(), l1.end());
            if (cutphase + eps == eps)
                break;
        }
    }
    mc = cutphasemin;

    //    cout << "minimum cut : " << cutphasemin << endl;
    //    cout << "Smin: " << endl;
    //    for (auto x : Smin) {
    //        cout << x << " ";
    //    }
    //    cout << endl;

    // ----------------------
    // grafo exemplo
    /*     for (i = 0; i < dim; i++) 
     *         delete[] sol[i];
     *     delete[] sol;
     *     sol = bkp;
     */


    // ----------------------
}

int BNC::CPXPUBLIC mycutcallback(CPXCENVptr env, void *cbdata, int wherefrom, 
        void *cbhandle, int *useraction_p) {
    using namespace std;

    bool frac;                  // solucao fracionaria
    double *X = NULL, ** sol = NULL, mc;
    int i, j, m, dim, ndepth, status = 0;

    set<int> Smin;

    // lock
    pthread_mutex_t cs_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&cs_mutex);
    *useraction_p = CPX_CALLBACK_DEFAULT;

    CPXgetcallbacknodeinfo (env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &ndepth);

    BNC * bnc = (BNC * ) cbhandle;
    dim = bnc->i->getDim();

    status = 0;

    X = new double[bnc->numCols];
    sol = new double * [dim];
    for (i = 0; i < dim; i++) {
        sol[i] = new double[dim];
    }

    CPXgetcallbacknodex(env, cbdata, wherefrom, X, 0, bnc->numCols - 1);
    bnc::storeLPSolution(env, bnc->model, bnc->numCols, X, sol);

    frac = bnc::isFrac(sol, dim);

//    if (ndepth <= 7 || !frac) {
//        bnc->maxBack(sol, Smin);
//        if (Smin.size() < 2) {
            bnc->minCut(sol, Smin, mc);
//        }
//    } else { 
//        bnc->minCut(sol, Smin, mc);
//    }
//
    if (Smin.size() < dim && mc + bnc::eps < 2) {
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
            cerr << "Corte nao violado: " << left << " " << rhs << endl;
        }

        if (CPXcutcallbackadd(env, cbdata, wherefrom, cutnz, rhs, 'L', cutInd, cutVal, 1)) {
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

/**
 * Cria modelo LP
 */
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
    CPXsetintparam(env, CPX_PARAM_MIPINTERVAL, 1);
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

    // Optimize
    double cplexTimeBefore, cplexTimeAfter;	
    CPXgettime(env, &cplexTimeBefore);
    CPXmipopt (env, model);
    CPXgettime(env, &cplexTimeAfter);
    double time = cplexTimeAfter - cplexTimeBefore;

    // saida do valor da funcao objetiva ao termino 
    // (vazio caso solucao nao encontrada)
    double objval;
    CPXgetobjval(env, model, &objval);
    ofstream f;
    std::string iOpt = "results/";
    iOpt += i->getName();
    f.open(iOpt.c_str());
    f << objval;
    f.close();

    //    bnc::printSolution(env, model, numCols, i->getDim());
    CPXfclose(fp);
    // Free
    CPXfreeprob(env, &model);
    CPXcloseCPLEX(&env);

    return 0;
}

