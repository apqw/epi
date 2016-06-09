#pragma once
#include <tbb/parallel_for.h>
#include <tbb/blocked_range3d.h>
#include <array>
#include <string>
#include <iostream>
#include <ostream>
#include <fstream>
template<typename T,size_t X,size_t Y,size_t Z>
class Field
{
	using Arr = T[X][Y][Z];
	Arr _fdata;
	
public:

	Field():_fdata{}{}
	void parallel_body_proc(T&& lmbd) {
		tbb::parallel_for(tbb::blocked_range3d<double>(0,X-1,0,Y-1,0,Z-1), lmbd);
	}

    inline Arr& operator()() {
		return _fdata;
	}

	inline const Arr& operator()()const {
		return _fdata;
	}

    void output_binary(const std::string& filename)const{
        std::ofstream ofs;
        ofs.open(filename,std::ios_base::out|std::ios_base::binary|std::ios_base::trunc);
        if(!ofs){
            std::cout<<"Field data output error. Filename:"<<filename<<std::endl;
            exit(1);
        }
        ofs.write(reinterpret_cast<const char*>(_fdata),sizeof(T)*X*Y*Z);
        ofs.close();
    }

    void read_binary(const std::string& filename){
        std::ifstream ifs;
        ifs.open(filename,std::ios_base::in|std::ios_base::binary);
        if(!ifs){
            std::cout<<"Field data input error. Filename:"<<filename<<std::endl;
            exit(1);
        }
        ifs.read(reinterpret_cast<char*>(_fdata),sizeof(T)*X*Y*Z);
        ifs.close();
    }

	//~Field();
};

