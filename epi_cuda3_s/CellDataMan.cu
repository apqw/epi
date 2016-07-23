#include "CellDataMan.h"
#include <cassert>
#include <typeinfo>
#include <iostream>
#include <fstream>
#include "fsys.h"
#include "cell_connection.h"
#include <device_functions.h>

void reset_index_map(IndexMap3D* imap){
	imap->_memset(MEMSET_NEGATIVE);
}

void reset_field_mask(FieldMask3D* fmask){
	fmask->_memset(0);
}

__global__ void cell_connection_device_init(devPtr<CellConnectionData> data){
	const int gid = blockDim.x*blockIdx.x + threadIdx.x;
	if (gid < MAX_CELL_NUM){
		CellConnectionData& cptr = thrust::raw_reference_cast(data[gid]);
		memset(cptr.connect_index,MEMSET_NEGATIVE, sizeof(CellIndex)*MAX_CONNECT_CELL_NUM); //negative value
		cptr.gj.init();
		//printf("initted %d\n", gid);
	}
}


__device__ void CellPackedInfo::set_info(CellPos cp, CellIndex i){
	x = __float2half(cp.x);
	y = __float2half(cp.y);
	z = __float2half(cp.z);
	idx = i;
}

__device__ float3 CellPackedInfo::get_pos()const{
	return make_float3(__half2float(x), __half2float(y), __half2float(z));
}

CellDataMan::CellDataMan() :ATP(ca2p_swt, _ATP), ca2p(ca2p_swt, _ca2p), IP3(ca2p_swt,_IP3)
, pos(swt, _pos), ext_stim(swt,_ext_stim)
{
	dbgprint("CellDataMan ctor called.\n");
}


CellDataMan::~CellDataMan()
{
	dbgprint("CellDataMan dtor called.\n");
}



void CellDataMan::cell_phase_switch(){
	return swt.switch_value();
}

void CellDataMan::ca2p_phase_switch(){
	return ca2p_swt.switch_value();
}

template<typename First>
bool __init_check_fn(int counter,const First& f){
	if (!f.has_initialized()){
		report_error([&]{
			std::cout
				<< "ERROR:Init Item No." << counter << " has not been initialized." << std::endl
				<< "Object Type:\t\t" << typeid(First).name() << std::endl
				<< "Object Address:\t\t" << &f << std::endl;
		});
		assert(f.has_initialized());
		exit_after_pause();
	}
	return f.has_initialized();
}

template<typename First,typename Second,typename... Rest>
bool __init_check_fn(int counter,const First& f, const  Second& s, const  Rest&... rest){
	return __init_check_fn(counter,f) && __init_check_fn(counter+1,s, rest...);
}
bool CellDataMan::check_initialized()const{
	return __init_check_fn(0,
		state,				//0
		fix_origin,			//1
		connection_data,	//2
		pair_index,			//3
		ca2p_avg,			//4
		ex_inert,			//5
		agek,				//6
		ageb,				//7
		ex_fat,				//8
		in_fat,				//9
		spr_nat_len,		//10
		//uid,				
		rest_div_times,		//11
		dermis_index,		//12
		zmax,				
		ncell,
		nder,
		nmemb,
		sw,
		//next_uid,
		ATP.current(),
		ATP.next(),
		ca2p.current(),
		ca2p.next(),
		IP3.current(),
		IP3.next(),
		pos.current(),
		pos.next(),
		ext_stim.current(),
		ext_stim.next(),
		field_mask,
		index_map);
}

__device__ __host__ void CellDataMan::_set_cell(
	CellIndex index,
	CELL_STATE cstate,
	int cfix_origin,
	real cx, real cy, real cz,
	//rad,
	real cca2p,
	real cca2p_avg,
	real cIP3,
	real cex_inert,
	real cagek,
	real cageb,
	real cex_fat,
	real cfat,
	real cspr_len,
	int cdiv_times,
	bool cis_malig
	){
	state[index] = cstate;
	fix_origin[index] = cstate==DER||cstate==MEMB?-1:cfix_origin;
	pos.set(index, make_cell_pos(cx, cy, cz));
	ca2p.set(index, cca2p);
	ca2p_avg[index] = cca2p_avg;
	IP3.set(index, cIP3);
	ex_inert[index] = cex_inert;
	agek[index] = cagek;
	ageb[index] = cageb;
	ex_fat[index] = cex_fat;
	in_fat[index] = cfat;
	spr_nat_len[index] = cspr_len;
	rest_div_times[index] = cstate==FIX?div_max:cdiv_times;
	is_malignant[index] = cis_malig;

}

void CellDataMan::_load_from_file(const std::string& path){
	std::ifstream dstrm;
	dstrm.open(path, std::ios::in);

	if (!dstrm) {
		report_error([&]{
			std::cerr << "Error: Failed to open cell data (in _load_from_file)." << std::endl
				<< "File path:\t\t" << path << std::endl;
		});
		assert(dstrm);
		exit_after_pause();
	}

	CELL_STATE _state = UNUSED;
	std::string line;
	int div_times = 0, touch = 0, pair_cell_id = 0, stem_orig_id = 0;
	real rad, _ageb, _agek, x, y, z, _fat, _spr_len, _ex_fat, _ca2p_avg, _ca2p;
	int id_count = 0;
	unsigned int phase = 0;
	int _nmemb = 0;
	int _nder = 0;

	while (std::getline(dstrm, line)) {

		sscanf(line.c_str(), "%*d " CELL_STATE_FMT " " R_FMT " " R_FMT " " R_FMT " " R_FMT " " R_FMT " " R_FMT " " R_FMT " " R_FMT " %d " R_FMT " " R_FMT " %d " R_FMT " %d %d",
			(CELL_STATE_t*)&_state, &rad, &_ageb, &_agek, &_ca2p, &x, &y, &z, &_ca2p_avg, &div_times, &_ex_fat, &_fat, &touch, &_spr_len, &pair_cell_id, &stem_orig_id);

		/*
		BLANK‚É‚½‚Ç‚è’…‚¢‚½‚çI—¹
		*/
		if (_state == BLANK)break;

		/*
		validation
		*/
		if (SYSTEM == BASAL && (_state == ALIVE || _state == DEAD || _state == AIR)) {
			printf(" input date must not contain ALIVE or DEAD in case of BASAL\n");
			exit(1);
		}
		if (_state == DER && rad != R_der) {
			printf("radii of DER not consistent with param.h\n");
			exit(1);
		}
		if (_state == MEMB && rad != R_memb) {
			printf("radii of DER not consistent with param.h\n");
			exit(1);
		}
		if (phase == 0 && _state != MEMB) {
			assert(_state == DER);
			phase++;
		}

		if (phase == 1 && _state != DER) {
			//assert(state == DER);
			phase++;
		}
		if (phase > 0 && _state == MEMB) {
			printf("non phase0 memb\n");
			exit(1);
		}
		if (phase > 1 && _state == DER) {
			printf("non phase1 der\n");
			exit(1);
		}

		//if (state == FIX)printf("FIX\n");

		if (_state == MEMB)_nmemb++;
		if (_state == DER)_nder++;

		_set_cell(id_count,
			_state,
			stem_orig_id,
			x, y, z,
			_ca2p, _ca2p_avg,
			IP3_init,
			ex_inert_init,
			_agek, _ageb,
			_ex_fat, _fat,
			_spr_len,
			div_times,
			stem_orig_id < MALIG_NUM
			);

		pair_index[id_count] = pair_cell_id >= 0 ? pair_cell_id : -1;
		printf("Phase %d  Cell loaded:%d\n", phase, id_count++);

	}
	nmemb = _nmemb; nmemb.mark_as_initialized();
	nder = _nder; nder.mark_as_initialized();
	ncell = id_count; ncell.mark_as_initialized();
}



void CellDataMan::init(const std::string& init_data_path, bool use_last, const std::string& init_uvp_data, 
	const std::string& init_w_data, const std::string& init_ATP_data, const std::string& init_ext_stim_data){
	
	state.fill(UNUSED);
	fix_origin.fill(-1);
	pos.foreach([] __device__ __host__ (decltype(pos.current())& _p){
		_p.fill(make_cell_pos(ZERO, ZERO, ZERO));
	});

	ca2p.foreach([] __device__ __host__ (decltype(ca2p.current())& _p){
		_p.fill(ca2p_init);
	});
	ca2p_avg.fill(ca2p_init);
	IP3.foreach([] __device__ __host__(decltype(IP3.current())& _p){
		_p.fill(IP3_init);
	});
	ex_inert.fill(ex_inert_init);
	agek.fill(ZERO);
	ageb.fill(ZERO);
	ex_fat.fill(ZERO);
	in_fat.fill(ZERO);
	spr_nat_len.fill(ZERO);
	rest_div_times.fill(0);
	is_malignant.fill(false);
	pair_index.fill(-1);
	//uid.fill(-1);

	_load_from_file(init_data_path);


	state.mark_as_initialized();
	fix_origin.mark_as_initialized();
	pos.foreach([]__device__ __host__(decltype(pos.current())& _p){
		_p.mark_as_initialized();
	});
	ca2p.foreach([]__device__ __host__(decltype(ca2p.current())& _p){
		_p.mark_as_initialized();
	});
	ca2p_avg.mark_as_initialized();
	IP3.foreach([]__device__ __host__(decltype(IP3.current())& _p){
		_p.mark_as_initialized();
	});
	ex_inert.mark_as_initialized();
	agek.mark_as_initialized();
	ageb.mark_as_initialized();
	ex_fat.mark_as_initialized();
	in_fat.mark_as_initialized();
	spr_nat_len.mark_as_initialized();
	rest_div_times.mark_as_initialized();
	is_malignant.mark_as_initialized();
	pair_index.mark_as_initialized();

	zmax = -REAL_MAX; zmax.mark_as_initialized();
	sw = 0; sw.mark_as_initialized();
	if (use_last){
		ATP.foreach([&] __device__ __host__(decltype(ATP.current())& _p){
			_p.read_binary(init_ATP_data);
		});

		ext_stim.foreach([&] __device__ __host__(decltype(ext_stim.current())& _p){
			_p.read_binary(init_ext_stim_data);
		});
	}
	else{
		ATP.foreach([] __device__ __host__(decltype(ATP.current())& _p){
			_p.fill(ATP_init);
		});

		ext_stim.foreach([] __device__ __host__(decltype(ext_stim.current())& _p){
			_p.fill(ext_stim_init);
		});
	}
	ATP.foreach([] __device__ __host__(decltype(ATP.current())& _p){
		_p.mark_as_initialized();
	});
	ext_stim.foreach([] __device__ __host__(decltype(ext_stim.current())& _p){
		_p.mark_as_initialized();
	});
	reset_field_mask(&field_mask);
	reset_index_map(&index_map);
	field_mask.mark_as_initialized();
	index_map.mark_as_initialized();

	
		cell_connection_device_init << <DEFAULT_THB_ALL_CELL >> >(connection_data);
		gpuErrchk(cudaDeviceSynchronize());
		//test
		for (int k = 0; k < 10000; k++){
			connect_cell(this);
			if (k % 100 == 0)printf("done %d\n", k);
		}
	connection_data.mark_as_initialized();
	//dbgprint("conn v:%d\n", ((CellConnectionData)connection_data[0]).connect_index[0]);
	pause();
}



