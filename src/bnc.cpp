/*
 * @file bnc.cpp
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 10/17/2014
 */

#include "../include/bnc.h"

ILOSTLBEGIN

typedef IloArray<IloNumArray>    FloatMatrix;
typedef IloArray<IloNumVarArray> NumVarMatrix;

int initBranchAndCut(const int ** matrix, unsigned dim) {
    int i, j;
    try {
        IloEnv env;
        IloModel model(env);
        IloArray<IloIntVarArray> x(env, dim + 1);
        IloExpr cost(env);

        for (i = 0; i < dim; i++) {
            IloIntVarArray array(env, dim + 1, 0, 2);
            x[i] = array;
            for (j = 0; j < dim; j++) {
                x[i][j].setBounds(0, 1);
            }
        }

        char var[100];
        for (i = 0; i < dim; i++) {
            for (j = 0; j < dim; j++) {
                sprintf(var, "X_%d_%d", i, j);
                x[i][j].setName(var);
                model.add(x[i][j]);
                cost += matrix[i][j] * x[i][j];
            }
        }

        model.add(IloMinimize(env, cost));

    } catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;   
    }
    catch(...) {
        cerr  << " ERROR" << endl;   
    }


}

