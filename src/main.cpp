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

#include "../include/BNC.h"

#include "Instance.h"
#include "Util.h"

std::string inline getInstanceName(std::string name) {
    size_t s = name.find("/") + 1;
    size_t e = name.find(".tsp") - 1;
    return name.substr(s, (e-s));
}

int main(int argc, char** argv) {
    std::string file = "";
    std::stringstream ss;
    int ub = 999999;

    if (argc >= 2) {
        file.assign(argv[1]);
        if (argc >= 3) {
            ss.str(argv[2]);
            ss >> ub;
        }
    } else {
        std::cout << "Uso " << argv[0] << " arquivo_entrada [ub]" << std::endl;
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
    BNC bnc(instance);
    bnc.initBranchAndCut(ub, getInstanceName(file));
    double aft = tsp::cpuTime();
    std::cout << (aft-bef) << "ms" << std::endl;

    return 0;
}
