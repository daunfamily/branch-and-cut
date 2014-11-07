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
    const double eps = 0.000001;
    const int inf = 999999;

    bool operator==(std::list<int>, std::list<int>);
    bool isFrac(double **sol, int dim);

    double updateEdges(double **, std::list<std::list<int>>&, std::list<int>&, std::list<int>&);
    double evalCutPhase(double **, std::list<std::list<int>>&, std::list<int>);

    std::list<std::list<int>>::iterator deleteNode(std::list<std::list<int>>& V, int a);
    std::list<std::list<int>>::iterator findNode(std::list<std::list<int>>&, int);

    std::string getColName(CPXCENVptr cpxEnv, CPXLPptr cpxModel, int colInd);

    void storeLPSolution(CPXCENVptr env, CPXLPptr model, int numcols, double *x, double **sol);
    void printSolution(CPXCENVptr env, CPXLPptr model, int cur_numcols, int N);
    void print(std::list<std::list<int>> g);
    void print(std::list<int> p);
}; 

#endif
