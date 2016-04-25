#include <iostream>
#include "define.h"
#include "update.h"
using namespace std;

int main(int argc, char *argv[])
{
    ENV* ei = new ENV(),*ei2=new ENV();
    ei->current_cell_num=100000;
    ei2->current_cell_num=100000;


    update_all(&ei,&ei2);
    std::cout<<"end"<<std::endl;
    return 0;
}
