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
    const std::string udelim = "\n\t\t";
    cp.add<std::string>("mode", 'm', 
        "Specify a mode shown below."_s+udelim+
        "c:Calculate with CPU." + udelim +
        "g:Calculate with GPU." + udelim +
        "v:Visualize the result." + udelim + 
        "d:Diff parameter file." + udelim + 
        "i:Generate default parameter files.");
    /*
    cp.add("cpu", 'c', "Calculate with CPU.");
    cp.add("gpu", 'g', "Calculate with GPU.");
    cp.add("visualize", 'v', "Visualize the result.");
    */
    cp.add<std::string>("cparam", '\0', "Specify an parameter file path for calculation.",false);
    cp.add<std::string>("cparam2", '\0', "Additional parameter file (for diff)",false);
    cp.add<std::string>("vparam", '\0', "Specify an parameter file path for visualization.",false);
    cp.parse_check(argc, argv);
    /*
    if (!cp.exist("mode")) {
        std::cerr << "Mode must be specified." << std::endl;
        std::cout<<cp.usage() << std::endl;
        std::exit(1);
    }
    */
    auto mode = cp.get<std::string>("mode");
    bool is_gpu = mode == "g";
    bool is_vis = mode == "v";
    bool is_cpu = mode == "c";
    bool is_diff = mode == "d";
    bool is_gen = mode == "i";
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
    if (is_gen) {
        if (cp.exist("cparam")) {
            std::cout << "Generate cparam as " << cp.get<std::string>("cparam")  << std::endl;
            CalcParams().generate_paramfile(cp.get<std::string>("cparam"));
        }

        if (cp.exist("vparam")) {
            std::cout << "Generate vparam as " << cp.get<std::string>("vparam") << std::endl;
            VisParams().generate_paramfile(cp.get<std::string>("vparam"));
        }

    }
    if (is_diff) {
        if (!cp.exist("cparam")|| !cp.exist("cparam2")) {
            std::cerr << "--cparam and --cparam2 must be specified." << std::endl;
            std::exit(1);
        }
        try {
            const CalcParams cp1 = CalcParams(cp.get<std::string>("cparam"));
            const CalcParams cp2 = CalcParams(cp.get<std::string>("cparam2"));
            std::cout << "Compare " << cp.get<std::string>("cparam") << " and " << cp.get<std::string>("cparam2") <<std::endl;
            CalcParams::diff(cp1, cp2);
        }
        catch (std::exception& e) {
            std::cerr << e.what() << std::endl;
            std::exit(1);
        }
    }
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
        CellManager cman;
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
