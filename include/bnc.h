/*
 * @file bnc.h
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 10/17/2014
 */

#ifndef BNC_H
#define BNC_H

#include <ilcplex/ilocplex.h>

int initBranchAndCut(const int ** matrix, unsigned dim);

#endif
