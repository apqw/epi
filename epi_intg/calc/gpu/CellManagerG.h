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
typedef
class CellManagerG {
	thrust::device_ptr<int> cell_count;
	thrust::device_ptr<int> nder;
	thrust::device_ptr<int> nmemb;
	SwapData<thrust::device_vector<real4>> pos;

};

#endif /* CELLMANAGERG_H_ */
