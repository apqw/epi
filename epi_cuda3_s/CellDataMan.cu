#include "CellDataMan.h"


CellDataMan::CellDataMan() :sw(0), next_uid(0), ncell(0), nder(0)
{
#ifdef DBG
	printf("CellDataMan ctor called.\n");
#endif
}


CellDataMan::~CellDataMan()
{
#ifdef DBG
	printf("CellDataMan dtor called.\n");
#endif
}


devPtr<CellPos> CellDataMan::current_pos_arr(){
	return pos[swt.current()];
}

devPtr<CellPos> CellDataMan::next_pos_arr(){
	return pos[swt.next()];
}