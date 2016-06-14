#include "cell_init.h"
#include "cell.h"
#include "cell_connection.h"
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <cinttypes>
#include <cstdio>


void _cman_value_init(CellManager & cman)
{
	using namespace cont;
	cman.all_foreach_parallel_native([](Cell*&c) {
	
		c->ca2p._set(ca2p_init);
		c->ca2p_avg = c->ca2p();
		c->ex_inert=ex_inert_init;
		c->IP3._set(IP3_init);
		c->connected_cell.foreach([&](Cell* conn) {
            c->gj.emplace(conn,gj_init);
		});
	});
	cman.other_foreach_parallel_native([](Cell*&c) {
	
		if (c->state == DEAD) {
			c->ca2p_avg = 0;
			c->ca2p._set(0.0);
		}
	});


}

void load_uvp(CellManager& cman,const std::string& ld_uvp){
    std::ifstream fuvp(ld_uvp);
    if(!fuvp){
        std::cout<<"Last uvp data open error. Filename:"<<ld_uvp<<std::endl;
        exit(1);
    }
    std::string ln;
    size_t lcount=0;
    int cidx=0;
    double ca2p,ex_inert,IP3;
    while(std::getline(fuvp,ln)){
        std::sscanf(ln.c_str(),"%d %lf %lf %lf",&cidx,&ca2p,&ex_inert,&IP3);
        Cell* c = cman[cidx];
        c->ca2p._set(ca2p);
        c->ex_inert=ex_inert;
        c->IP3._set(IP3);
        lcount++;
    }
    if(lcount!=cman.size()){
        std::cout<<"Num of cell mismatch. num in file ="<<lcount<<",current="<<cman.size()<<std::endl;
        exit(1);
    }


}


void load_w(CellManager& cman,const std::string& ld_w){
    std::ifstream fw_alt(ld_w);
    if(!fw_alt){
        std::cout<<"Last w data open error. Filename:"<<ld_w<<std::endl;
        exit(1);
    }

std::string ln;
    size_t lcount=0;int cidx=0;
    int gj_len=0;
    while(std::getline(fw_alt,ln)){
        if(std::sscanf(ln.c_str(),"%d %d",&cidx,&gj_len)!=2){
            std::cout<<"bad GJ data format (index)"<<std::endl;
            exit(1);
        }
        Cell* c=cman[cidx];
        int gj_idx=0;double gj_val=0;
        for(int i=0;i<gj_len;i++){
            if(!std::getline(fw_alt,ln)){
                std::cout<<"bad GJ data format"<<std::endl;
                exit(1);
            }
            std::sscanf(ln.c_str(),"%d %lf",&gj_idx,&gj_val);
            c->gj[cman[gj_idx]]()=gj_val;
        }
        lcount++;
    }
    if(lcount!=cman.size()){
        std::cout<<"Num of cell mismatch. num in file ="<<lcount<<",current="<<cman.size()<<std::endl;
        exit(1);
    }
}

void load_from_last(CellManager& cman,const std::string& ld_uvp,const std::string& ld_w){
    load_uvp(cman,ld_uvp);
    load_w(cman,ld_w);
}

//////////////////////////////////////////////////////////////////

void cman_init(CellManager& cells,const std::string& init_data_path,
               bool use_last,const std::string& ld_uvp,const std::string& ld_w) {
	cells.init_internal(init_data_path);
	connect_cell(cells);
	_cman_value_init(cells);
    if(use_last){
        load_from_last(cells,ld_uvp,ld_w);
    }
}
