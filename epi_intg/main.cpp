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
#include "vis/visualize.h"
#include "vis/VisParams.h"

int main(int argc,char** argv){
    cmdline::parser cp;

    cp.add<std::string>("mode", 'm', "Specify a mode shown below. \n\t\tc:Calculate with CPU.\n\t\tg:Calculate with GPU.\n\t\tv:Visualize the result.");
    /*
    cp.add("cpu", 'c', "Calculate with CPU.");
    cp.add("gpu", 'g', "Calculate with GPU.");
    cp.add("visualize", 'v', "Visualize the result.");
    */
    cp.add<std::string>("cparam", '\0', "Specify an parameter file path for calculation.");
    cp.add<std::string>("vparam", '\0', "Specify an parameter file path for visualization.");
    cp.parse_check(argc, argv);
    if (!cp.exist("mode")) {
        std::cerr << "Mode must be specified." << std::endl;
        std::exit(1);
    }
    auto mode = cp.get<std::string>("mode");
    bool is_gpu = mode == "g";
    bool is_vis = mode == "v";
    bool is_cpu = mode == "c";
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
        if (!cp.exist("cparam")) {
            std::cerr << "A parameter file must be specified." << std::endl;
            std::exit(1);
        }

        try {
            pm = std::unique_ptr<const CalcParams>(new const CalcParams(cp.get<std::string>("cparam")));
        }
        catch (std::exception& e) {
            std::cerr << e.what() << std::endl;
            std::exit(1);
        }
        pm->generate_paramfile("opop");
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

    if (is_vis) {
        if (!cp.exist("cparam")) {
            std::cerr << "A parameter file must be specified." << std::endl;
            std::exit(1);
        }
        try {
            pm = std::unique_ptr<const CalcParams>(new const CalcParams(cp.get<std::string>("cparam")));
        }
        catch (std::exception& e) {
            std::cerr << e.what() << std::endl;
            std::exit(1);
        }
       // pm->generate_paramfile("opop2");

        if (!cp.exist("vparam")) {
            std::cerr << "A vis parameter file must be specified." << std::endl;
            std::exit(1);
        }
        
        try {
            VisParams vp(cp.get<std::string>("vparam"));
            init_visualizer(&argc, argv,vp);
            visualize();
        }
        catch (std::exception& e) {
            std::cerr << e.what() << std::endl;
            std::exit(1);
        }
    }




	std::cout<<"ok"<<std::endl;
}
