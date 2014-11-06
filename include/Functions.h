/*
 * @file Functions.h
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 10/22/2014
 */

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <string>
#include <list>
#include <ilcplex/ilocplex.h>

namespace bnc {
    const double eps = 0.001;
    const int inf = 99999;
    std::string getColName(CPXCENVptr cpxEnv, CPXLPptr cpxModel, int colInd);
    void storeLPSolution(CPXCENVptr env, CPXLPptr model, int numcols, double *x, double **sol);
    void printSolution(CPXCENVptr env, CPXLPptr model, int cur_numcols, int N);
}; 

#endif
