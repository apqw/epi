/*
 * init_gen.cpp
 *
 *  Created on: 2016/10/04
 *      Author: yasu7890v
 */




#include "init_gen.h"
#include "../global.h"
#include "../define.h"
#include "cpu2/CellManager.h"
#include <iostream>
#include <functional>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


static void relaxation(CellManager& cman,int nfix,double Z_LEVEL){
	const size_t NB=cman.size()-nfix;
	double distlj;
	for(int j=NB;j<NB+nfix;j++){
		Cell* fixp=cman[j];
		if(fixp->state!=FIX){
			throw std::logic_error("Non-FIX cell found at last section.");
		}
double distmin=10000.0;
		while(1){
			for(size_t l=0;l<NB;l++){
				if(cman[l]->state!=MEMB)break;
				distlj=sqrt(p_cell_dist_sq(fixp,cman[l]));
				if(distlj<distmin)distmin=distlj;
			}
			if(distmin<1.1*(fixp->radius+pm->R_der)){
				break;
			}else{
				if(fixp->z()<Z_LEVEL-pm->R_max){
					throw std::logic_error("zz too low: zz="_s
							+std::to_string(fixp->z())+" zz-zlevel="+std::to_string(fixp->z()-Z_LEVEL));
				}
				fixp->z._set(fixp->z()-0.001);
			}
		}
	}
}
void init_gen(const std::string& filename,int nfix,int der){

	CellManager cman;
	using namespace std::placeholders;
	auto cgen=std::bind(&CellManager::create_resizable,&cman,_1,_2,_3,_4,_5,_6,0,0,0,0,0,0,0,0,0,0,0);
	double pratio=2*pm->P_MEMB;
	const double Z_LEVEL = der*2.0*pm->R_der;

	if(!pm->USE_TRI_MEMB){
		for(unsigned int k=0;k<pm->MEMB_NUM_Y;k++){
			for(unsigned int j=0;j<pm->MEMB_NUM_X;j++){
				const double xx=0.5*pm->R_memb+j*pm->R_memb*pratio;
				const double yy=0.5*pm->R_memb+k*pm->R_memb*pratio;
				const double zz=Z_LEVEL+pm->R_memb;
				cgen(MEMB,0,xx,yy,zz,pm->R_memb);
			}
		}
	}else{
		for(unsigned int k=0;k<pm->MEMB_NUM_Y;k++){
			for(unsigned int j=0;j<pm->MEMB_NUM_X;j++){
				const double xx=0.5*pm->R_memb+j*pm->R_memb*pratio+(k%2==0?0.0:pm->R_memb*pratio*0.5);
				const double yy=0.5*pm->R_memb+k*pm->R_memb*pratio/pm->y_tri_comp_ratio;
				const double zz=Z_LEVEL+pm->R_memb;
				cgen(MEMB,0,xx,yy,zz,pm->R_memb);
			}
		}
	}


	const int NX_DER=pm->LX/(2.0*pm->R_der);
	const int NY_DER=pm->LY/(2.0*pm->R_der);
	for(int l=0;l<der;l++){
		for(unsigned int k=0;k<NY_DER;k++){
			for(unsigned int j=0;j<NX_DER;j++){
				const double xx=pm->R_der+2.0*j*pm->R_der;
				const double yy=pm->R_der+2.0*k*pm->R_der;
				const double zz=pm->R_der+2.0*l*pm->R_der;
				cgen(DER,0,xx,yy,zz,pm->R_der);
			}
		}
	}

	for(int l=0;l<nfix;l++){
		const double xx=pm->LX/2.0+(pm->LX/4.0)*cos(2.0*M_PI*l/nfix);
		const double yy=pm->LX/2.0+(pm->LX/4.0)*sin(2.0*M_PI*l/nfix);
		const double zz=pm->R_max+pm->R_memb+Z_LEVEL+pm->R_memb;
		cgen(FIX,l,xx,yy,zz,pm->R_max);
	}

	relaxation(cman,nfix,Z_LEVEL);

	cman.output(filename);

}
