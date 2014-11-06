/*
 * @file BNC.h
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 10/17/2014
 */

#ifndef BNC_H
#define BNC_H

#include <algorithm>
#include <list>
#include <pthread.h>
#include <set>
#include <string>

#include <ilcplex/ilocplex.h>

#include "Functions.h"
#include "Instance.h"

class BNC {
    public:
        BNC() = delete;
        BNC(Instance&);
        virtual ~BNC();
        int initBranchAndCut(int ub, std::string instanceName);
        void maxBack(double ** sol, std::set<int>& Smin);
        void minCut(double ** sol, std::set<int>& Smin);
        void createLP(const int ** matrix, unsigned dim);

    private:
        CPXENVptr env;
        CPXLPptr  model;
        Instance* i;

        int numCols;

        static int CPXPUBLIC mycutcallback(CPXCENVptr env, void *cbdata, int wherefrom, 
                void *cbhandle, int *useraction_p);

};

#endif
