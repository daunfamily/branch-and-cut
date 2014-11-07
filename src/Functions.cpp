/*
 * @file Functions.h
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 10/22/2014
 */

#include "../include/Functions.h"

namespace bnc {

    bool operator==(std::list<int> l, std::list<int> r) {
        l.sort();
        r.sort();

        for (auto i : l) {
            for (auto j : r) {
                if (i != j)
                    return false;
            }   
        }
        return true;
    }

    bool isFrac(double ** sol, int dim) {

        bool flagFrac = false;

        for (int i = 0; i < dim + 1; i++) {
            for (int j = i + 1; j < dim + 1; j++) {
                if ((sol[i][j] > eps && sol[i][j] < 1 - eps) || 
                        (sol[i][j] > 1 + eps && sol[i][j] < 2 - eps)) {
                    flagFrac = true;
                }
            }
        }
        return flagFrac;
    }

    double updateEdges(double ** sol, std::list<std::list<int>>& g, std::list<int>& l1, std::list<int>& l2) {
        double e = 0.0, d, i;
        std::list<std::list<int>>::iterator git, itP1, itP2;
        std::list<int>::iterator it, jt;

        itP1 = findNode(g, l1.front());
        itP2 = findNode(g, l2.front());

        for (git = g.begin(); git != g.end(); ++git) {
            if (git == itP1 || git == itP2) {
                continue;
            }

            if (git->front() < itP1->front())
                e = sol[git->front()][itP1->front()];
            else
                e = sol[itP1->front()][git->front()];

            for (it = git->begin(); it != git->end(); ++it) {
                for (jt = itP2->begin(); jt != itP2->end(); ++jt) {
                    if (*it < *jt)
                        d = sol[*it][*jt] += e;
                    else
                        d = sol[*jt][*it] += e;
                }

                for (jt = itP1->begin(); jt != itP1->end(); ++jt) {
                    if (*it < *jt)
                        sol[*it][*jt] = e;
                    else
                        sol[*jt][*it] = e;
                }
            }
        }
    }

    std::list<std::list<int>>::iterator findNode(std::list<std::list<int>>& V, int a) {
        std::list<int>::iterator it;
        std::list<std::list<int>>::iterator git;

        for (git = V.begin(); git != V.end(); ++git) {
            for (it = git->begin(); it != git->end(); ++it) {
                if (*it == a)
                    return git;
            }
        }
        return V.end();
    }

    std::list<std::list<int>>::iterator deleteNode(std::list<std::list<int>>& V, int a) {
        std::list<std::list<int>>::iterator git = findNode(V, a);

        if (git != V.end()) {
            return V.erase(git);
        }
        return V.end();
    }

    double evalCutPhase(double ** sol, std::list<std::list<int>>& g, std::list<int> p) {
        double cut = 0.0;
        std::list<std::list<int>>::iterator git;

        for (git = g.begin(); git != g.end(); ++git) {
            if (git->front() < p.front())
                cut += sol[git->front()][p.front()];
            else
                cut += sol[p.front()][git->front()];
        }
        return cut;
    }

    std::string getColName(CPXCENVptr cpxEnv, CPXLPptr cpxModel, int colInd) {
        int surplus;
        CPXgetcolname (cpxEnv, cpxModel, NULL, NULL, 0, &surplus, colInd, colInd);
        int storespace = - surplus;

        char* aux = new char[storespace];
        char** colname = new char*[1];
        strcat(colname[0],"");

        int largestNameLength = 100;
        CPXgetcolname(cpxEnv,cpxModel,colname,aux,largestNameLength,&surplus,colInd,colInd);

        std::string result(colname[0]);

        delete[] colname;
        delete[] aux;
        return result;
    }

    void storeLPSolution(CPXCENVptr env, CPXLPptr model, int numcols, double *x, double **sol) {
        char name[100];
        for (int k = 0; k < numcols; k++) { 
            std::string s = getColName(env, model, k);
            int i,j;
            if (s.at(0) == 'X') {
                strncpy(name,s.c_str(),100); 
                sscanf(name,"X_%d_%d",&i,&j);
                sol[i][j] = x[k];
            }
        }
    }

    void printSolution(CPXCENVptr env, CPXLPptr model, int cur_numcols, int N) {
        using namespace std;

        double *x = new double[cur_numcols];
        CPXgetx (env, model, x, 0, cur_numcols-1);

        double **sol = new double*[N+1];

        for (int i = 0; i < N + 1; i++) {
            sol[i] = new double[N + 1];
            for (int j = 0; j < N + 1; j++) {
                sol[i][j] = 0;
            }
        }

        char varName[100];
        cout << "\n\nSolution:\n " << endl;
        for (int k = 0; k < cur_numcols; k++) { 
            string s = getColName(env, model, k);
            int i,j;
            if (s.at(0) == 'X') {
                strncpy(varName,s.c_str(),100); 
                sscanf(varName,"X_%d_%d", &i, &j);
                sol[i][j] = x[k];
                if (sol[i][j] > 0.01) {
                    cout << "X_" << i << "_" << j << " " << sol[i][j] << endl;
                }
            }
        }

        delete[] x;
        for (int i = 0; i < N + 1; i++) {
            delete[] sol[i];
        }
        delete[] sol;
    }

    /**
     * Imprime grafo
     */
    void print(std::list<std::list<int>> g) {
        for (auto gi : g) {
            for (auto i : gi) {
                std::cout << i << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    void print(std::list<int> p) {
        for (auto i : p) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }

};
