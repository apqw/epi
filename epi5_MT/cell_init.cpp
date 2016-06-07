#include "cell_init.h"
#include "cell.h"
#include "cell_connection.h"
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <cinttypes>



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

//////////////////////////////////////////////////////////////////

void cman_init(CellManager& cells, std::string init_data_path) {
	cells.init_internal(init_data_path);
	connect_cell(cells);
	_cman_value_init(cells);
}
