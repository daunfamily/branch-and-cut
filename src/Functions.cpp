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
                sscanf(name,"X(%d,%d)",&i,&j);
                sol[i][j] = x[k];
            }
        }
    }
};
