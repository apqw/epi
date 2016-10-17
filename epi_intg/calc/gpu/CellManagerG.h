/*
 * CellManagerG.h
 *
 *  Created on: 2016/10/14
 *      Author: yasu7890v
 */

#ifndef CELLMANAGERG_H_
#define CELLMANAGERG_H_
#include<thrust/device_vector.h>
#include<thrust/device_ptr.h>
#include "../../cuda_define.h"
#include "../../misc/swapdata.h"

struct GCellConnectionData{
	int connect_num;
	int connect_index[G_MAX_CONNECT_CELL_NUM];
};

struct CellData_phys{
	//CELL_STATE state;
	GCellConnectionData connect_data;
   real radius;
};

struct CellData_chem{
    //real ca2p;
    real ca2p_avg;
    //real IP3;
    real ex_inert;
    real agek;
    real ageb;
    real ex_fat;
    real in_fat;
    real spr_nat_len;

};
class GCellManager {
	thrust::device_ptr<int> cell_count;
	thrust::device_ptr<int> nder;
	thrust::device_ptr<int> nmemb;

	SwapData<thrust::device_vector<real4>> pos;

};

#endif /* CELLMANAGERG_H_ */
