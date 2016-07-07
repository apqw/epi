/*
 * main.cpp
 *
 *  Created on: 2016/07/07
 *      Author: yasu7890v
 */

#include "code_test.h"
#include "stdio.h"
#include "CellManager.h"
#include "cell_connection.h"
int main(int argc,char** argv){

	codetest();

	CellManager cm;
	cm.init_with_file(argv[1]);
	printf("loop start\n");
	for(int i=0;i<100000;i++){
		connect_cell(&cm);
		printf("tesuya %d\n",i);
	}
}
