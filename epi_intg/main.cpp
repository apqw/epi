/*
 * main.cpp
 *
 *  Created on: 2016/09/20
 *      Author: yasu7890v
 */



#include <iostream>
#include "Params.h"
#include "global.h"
#include "parser.h"
#include "cmdline/cmdline.h"
#include "calc/cpu2/atomics.h"
#include "calc/cpu2/CellManager.h"
#include "calc/cpu2/cell_connection.h"

int main(int argc,char** argv){
    cmdline::parser cp;
    cp.add("cpu", 'c', "Calculate with CPU.");
    cp.add("gpu", 'g', "Calculate with GPU.");
    cp.add("visualize", 'v', "Visualize the result.");

    cp.add<std::string>("param", '\0', "Specify an parameter file path for calculation.");

    

    cp.parse_check(argc, argv);

    bool is_gpu = cp.exist("gpu");
    bool is_vis = cp.exist("visualize");
    bool is_cpu = cp.exist("cpu") || (!is_gpu && !is_vis);

    Lockfree_push_stack_dyn<size_t> oooo(100);
#pragma omp parallel for
    for (int i = 0; i < 50; i++) {
        oooo.push_back(i);
        std::cout << i << std::endl;
    }

    oooo.test_realloc();
    std::cout << oooo.max_size() << ":"<<oooo.size()<<std::endl;
    for (int i = 0; i < 41; i++) {
        oooo.push_back(1);
    }
    oooo.test_realloc();
    std::cout << oooo.max_size() << std::endl;
    if (is_cpu || is_gpu) {
        if (!cp.exist("param")) {
            std::cerr << "A parameter file must be specified." << std::endl;
            std::exit(1);
        }

        try {
            pm = std::unique_ptr<const Params>(new const Params(cp.get<std::string>("param")));
        }
        catch (std::exception& e) {
            std::cerr << e.what() << std::endl;
            std::exit(1);
        }

        CellManager cman(100000);
        try {
            cman.load("../epi5_MT/input_2layer_nfix16.dat");
        }
        catch (std::exception& e) {
            std::cerr << e.what() << std::endl;
            std::exit(1);
        }
        connect_cell(cman);
        
    }




	std::cout<<"ok"<<std::endl;
}
