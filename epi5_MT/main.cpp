#include <iostream>
#include "codetest.h"
#include "proc.h"
#include "define.h"

#include <string>
using namespace std;

int main(int argc, char *argv[])
{
    if(argc<=1){
        std::cout<<"No input file specified."<<std::endl;
        exit(1);
    }

    if(argc<=2){
        std::cout<<"No last_data_mode specified."<<std::endl;
        exit(1);
    }

    if (argc <= 3) {
        std::cout << "No cornif_mode specified." << std::endl;
        exit(1);
    }

    atomic_double_test();
    lfpstack_test();
	cell_test();
    proc(argv[1],std::stoi(argv[3])==1,std::stoi(argv[2])==1,last_data_uvp_name,last_data_w_name,last_data_a_name,last_data_B_name );
#ifdef _WIN32
	system("pause");
#endif
    return 0;
}
