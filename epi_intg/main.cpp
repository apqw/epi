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
#include "calc/calc.h"
#include "calc/cpu2/atomics.h"
#include "calc/cpu2/CellManager.h"
#include "calc/cpu2/cell_connection.h"
#include "vis/visualize.h"
#include "vis/VisParams.h"
#include "calc/init_gen.h"
#include "util/vec/Vec.h"
struct ErrorHandling{


	cmdline::parser& cp;

	template<typename First>
	std::string checkstr(First f){
		return f;
	}
	template<typename First,typename Second,typename...Rest>
	std::string checkstr(First f,Second s,Rest...r){
		return std::string(f)+" and "+checkstr(s,r...);
	}

	template<typename First>
	bool checkexist(First f){
		return cp.exist(f);
	}
	template<typename First,typename Second,typename...Rest>
	bool checkexist(First f,Second s,Rest...r){
		return cp.exist(f)&&checkexist(s,r...);
	}
	template<typename...Str>
	void check(Str... r){

		if(!checkexist(r...)){
			std::cerr<<checkstr(r...)<<" must be specified."<<std::endl;
			std::exit(1);
		}
	}
	ErrorHandling(cmdline::parser&_cp):cp(_cp){}
};

void set_param(const cmdline::parser& cp){
    try {
         pm = std::unique_ptr<const CalcParams>(new const CalcParams(cp.get<std::string>("cparam")));
     }
     catch (std::exception& e) {
         std::cerr << e.what() << std::endl;
         std::exit(1);
     }
}

int main(int argc,char** argv){

    cmdline::parser cp;
    const std::string udelim = "\n\t\t";
    cp.add<std::string>("mode", 'm',
        "Specify a mode shown below."_s+udelim+
        "c:Calculate with CPU. Need --cparam,--input" + udelim +
        "g:Calculate with GPU. Need --cparam,--input" + udelim +
        "v:Visualize the result.Need --cparam,--vparam" + udelim +
        "d:Diff calc parameter file. Need --cparam,--cparam2" + udelim +
        "gp:Generate default parameter files. Use --cparam,--vparam for output."+udelim+
        "i:Generate initial cell data file.Need --cparam,--output");
    cp.add<std::string>("input", 'i', "Cell data file.", false);
    cp.add<std::string>("output",'o', "Cell data output file(for init)",false);
    cp.add<std::string>("cparam", '\0', "Specify an parameter file path for calculation.",false);
    cp.add<std::string>("cparam2", '\0', "Additional parameter file (for diff)",false);
    cp.add<std::string>("vparam", '\0', "Specify an parameter file path for visualization.",false);
    cp.parse_check(argc, argv);
    auto mode = cp.get<std::string>("mode");
    bool is_gpu = mode == "g";
    bool is_vis = mode == "v";
    bool is_cpu = mode == "c";
    bool is_diff = mode == "d";
    bool is_gen = mode == "gp";
    bool is_init = mode == "i";


    auto EH = ErrorHandling(cp);

    if(is_init){
    	EH.check("cparam","output");
    	set_param(cp);
    	try{
    	init_gen(cp.get<std::string>("output"),16,8);
    	}catch(std::exception& e){
    		  std::cerr << e.what() << std::endl;
    		            std::exit(1);
    	}

    }
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
    	EH.check("cparam","cparam2");
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

    	EH.check("cparam","input");
    	set_param(cp);
        CellManager cman;
        try {
            cman.load(cp.get<std::string>("input"));
        }
        catch (std::exception& e) {
            std::cerr << e.what() << std::endl;
            std::exit(1);
        }
        //connect_cell(cman);
        calc_with_cpu(cman);
    }

    if (is_vis) {
       EH.check("cparam","vparam");
       set_param(cp);
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




	std::cout<<"Done."<<std::endl;
}
