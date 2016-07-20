#pragma once

#include <iostream>

void pause();
void make_dir(const char* dir);
void exit_after_pause();

template<class F>
void report_error(F&& lmbd){
	std::cerr 
		<< std::endl 
		<< "/////////////////////////////// ERROR REPORT ///////////////////////////////" 
		<< std::endl
		<< std::endl;
	lmbd();
	std::cerr 
		<< std::endl 
		<< "/////////////////////////// END OF ERROR REPORT ////////////////////////////" 
		<< std::endl 
		<< std::endl;
}

template<class F>
void report_warn(F&& lmbd){
	std::cerr
		<< std::endl
		<< "WARN:"
		<< std::endl
		<< std::endl;
	lmbd();
	std::cerr
		<< std::endl
		<< std::endl;
}