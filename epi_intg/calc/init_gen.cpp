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
	for(size_t j=NB;j<NB+nfix;j++){
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
				fixp->z._set(fixp->z()-real(0.001));
			}
		}
	}
}
CellManager init_gen(int nfix,int der){

	CellManager cman;
	using namespace std::placeholders;
	auto cgen=std::bind(&CellManager::create_resizable,&cman,_1,_2,_3,_4,_5,_6,
        real(0.0), real(0.0), real(0.0), real(0.0), real(0.0), real(0.0), real(0.0), real(0.0), real(0.0),0,false);
	double pratio=2*pm->P_MEMB;
	const double Z_LEVEL = der*2.0*pm->R_der;

	unsigned int nmemb=0,nder=0;
	if(!pm->USE_TRI_MEMB){
		for(unsigned int k=0;k<pm->MEMB_NUM_Y;k++){
			for(unsigned int j=0;j<pm->MEMB_NUM_X;j++){
				const real xx=real(0.5*pm->R_memb+j*pm->R_memb*pratio);
				const real yy= real(0.5*pm->R_memb+k*pm->R_memb*pratio);
				const real zz= real(Z_LEVEL+pm->R_memb);
				cgen(MEMB,0,xx,yy,zz,pm->R_memb);
			}
		}
	}else{
		for(unsigned int k=0;k<pm->MEMB_NUM_Y;k++){
			for(unsigned int j=0;j<pm->MEMB_NUM_X;j++){
				const real xx= real(0.5*pm->R_memb+j*pm->R_memb*pratio+(k%2==0?0.0:pm->R_memb*pratio*0.5));
				const real yy= real(0.5*pm->R_memb+k*pm->R_memb*pratio/pm->y_tri_comp_ratio);
				const real zz= real(Z_LEVEL+pm->R_memb);
				cgen(MEMB,0,xx,yy,zz,pm->R_memb);
			}
		}
	}
	cman.nmemb=pm->MEMB_NUM_Y*pm->MEMB_NUM_X;

	const unsigned int NX_DER=static_cast<unsigned int>(pm->LX/(2.0*pm->R_der));
    const unsigned int NY_DER = static_cast<unsigned int>(pm->LY / (2.0*pm->R_der));
	for(int l=0;l<der;l++){
		for(unsigned int k=0;k<NY_DER;k++){
			for(unsigned int j=0;j<NX_DER;j++){
				const real xx= real(pm->R_der+2.0*j*pm->R_der);
				const real yy= real(pm->R_der+2.0*k*pm->R_der);
				const real zz= real(pm->R_der+2.0*l*pm->R_der);
				cgen(DER,0,xx,yy,zz,pm->R_der);
			}
		}
	}
	cman.nder=der*NY_DER*NX_DER;
	for(int l=0;l<nfix;l++){
		const real xx= real(pm->LX/2.0+(pm->LX/4.0)*cos(2.0*M_PI*l/nfix));
		const real yy= real(pm->LY/2.0+(pm->LY/4.0)*sin(2.0*M_PI*l/nfix));
		const real zz= real(pm->R_max+pm->R_memb+Z_LEVEL+pm->R_memb);
		cgen(FIX,l,xx,yy,zz,pm->R_max);
	}

	relaxation(cman,nfix,Z_LEVEL);

	return cman;
	//cman.output(filename);

}

void init_gen_output(const std::string& filename,int nfix,int der){
	init_gen(nfix,der).output(filename);
}
