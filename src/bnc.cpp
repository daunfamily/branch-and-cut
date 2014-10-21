/*
 * @file bnc.cpp
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 10/17/2014
 */

#include "../include/bnc.h"

ILOSTLBEGIN

void createLP(const int ** matrix, unsigned dim) {
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

int initBranchAndCut(const int ** matrix, unsigned dim) {
    createLP(matrix, dim);


    return 0;
}

