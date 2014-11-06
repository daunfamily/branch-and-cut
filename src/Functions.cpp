/*
 * @file Functions.h
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 10/22/2014
 */

#include "../include/Functions.h"

namespace bnc {

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


};
