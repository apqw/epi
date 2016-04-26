#include <iostream>
#include "define.h"
#include "update.h"
using namespace std;

int main(int argc, char *argv[])
{
alignas(32)    ENV* ei = new ENV(),*ei2=new ENV();
    ei->current_cell_num=3000;
    ei2->current_cell_num=3000;

	for (int i = 0; i < CONST::NUM_ITR; i++) {
		printf("%d/%d\n", i, CONST::NUM_ITR);
		update_all(&ei, &ei2);
	}
    std::cout<<"end"<<std::endl;
    return 0;
}
