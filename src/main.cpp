/*
 * @file main2.cpp
 *      
 * @author: Fabio Cruz B. de Albuquerque
 * 
 * @date: 10/20/2014
 */

#include <exception>
#include <iostream>
#include <sstream>
#include <string>

#include "../include/bnc.h"

#include "Instance.h"
#include "Util.h"

int main(int argc, char** argv) {
    std::string file = "";
    std::stringstream ss;

    if (argc >= 2) {
        file.assign(argv[1]);
    } else {
        std::cout << "Uso " << argv[0] << " arquivo_entrada" << std::endl;
        return 1;
    }

    srand(time(NULL));
    Instance instance(file);

    try {
        instance.readInfo(tsp::TSP);
        instance.printInfo();
        instance.initMatrix();
        instance.readMatrixData();
        instance.closeFile();
    } catch (std::exception &e) {
        std::cout << e.what() << std::endl;
        std::cout << "Falha na inicializacao... encerrando." << std::endl;
        return 1;
    }

    double bef = tsp::cpuTime();
    initBranchAndCut(instance.getMatrix(), instance.getDim());
    double aft = tsp::cpuTime();
    std::cout << (aft-bef) << "ms" << std::endl;




    return 0;
}
