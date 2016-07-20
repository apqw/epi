#include "switcher.h"


switcher::switcher() :_swt(0){}
void switcher::switch_value(){
	_swt = 1 - _swt;
}

int switcher::current()const{
	return _swt;
}

int switcher::next()const{
	return 1 - _swt;
}