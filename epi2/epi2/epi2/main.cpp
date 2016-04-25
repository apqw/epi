#include <iostream>
#include "define.h"
#include "update.h"

int main() {
	ENV* newe = new ENV();
	ENV* ee=new ENV(), *ee2 = new ENV();
	newe->Ca2P_value[0][0][0] = 19191;
	ee->current_cell_num = 1000;
	ee2->current_cell_num = 1000;
	update_cell_internal(*ee, *ee2);
	std::cout << newe->Ca2P_value[0][0][0] <<std::endl;
	system("pause");
}