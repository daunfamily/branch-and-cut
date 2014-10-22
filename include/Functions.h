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
#include <ilcplex/ilocplex.h>

namespace bnc {
    std::string getColName(CPXCENVptr cpxEnv, CPXLPptr cpxModel, int colInd);
    void storeLPSolution(CPXCENVptr env, CPXLPptr model, int numcols, double *x, double **sol);
}; 

#endif
