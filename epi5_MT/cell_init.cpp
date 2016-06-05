#include "cell_init.h"
#include "cell.h"
#include "cell_connection.h"
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <cinttypes>



void _cman_memb_init(CellManager& cells)
{
	using namespace cont;
	cells.memb_foreach_with_index([&](CellPtr& cptr, size_t j) {
		int jj = j%NMX;
		int kk = j / NMX;
		if (jj == 0) {
			cptr->md.memb_l = cells[j + NMX - 1];
		}
		else {
			cptr->md.memb_l = cells[j - 1];
		}

		if (jj <= 1) {
			cptr->md.memb_ll = cells[j + NMX - 2];
		}
		else {
			cptr->md.memb_ll = cells[j - 2];
		}

		if (jj == NMX - 1) {
			cptr->md.memb_r = cells[j - (NMX - 1)];
		}
		else {
			cptr->md.memb_r = cells[j + 1];
		}

		if (kk == 0) {
			cptr->md.memb_b = cells[j + NMX*NMY - NMX];
		}
		else {
			cptr->md.memb_b = cells[j - NMX];
		}

		if (kk <= 1) {
			cptr->md.memb_bb = cells[j + NMX*NMY - 2 * NMX];
		}
		else {
			cptr->md.memb_bb = cells[j - 2 * NMX];
		}

		if (kk == NMY - 1) {
			cptr->md.memb_u = cells[j - (NMX*NMY - NMX)];
		}
		else {
			cptr->md.memb_u = cells[j + NMX];
		}

		cptr->connected_cell.push_back(cptr->md.memb_l);
		cptr->connected_cell.push_back(cptr->md.memb_r);
		cptr->connected_cell.push_back(cptr->md.memb_b);
		cptr->connected_cell.push_back(cptr->md.memb_u);
		//test
		/*
		cptr->gj._emplace(cptr->md.memb_l, 0);
		cptr->gj._emplace(cptr->md.memb_r, 0);
		cptr->gj._emplace(cptr->md.memb_b, 0);
		cptr->gj._emplace(cptr->md.memb_u, 0);
		*/
	});
}

void _cman_value_init(CellManager & cman)
{
	using namespace cont;
	cman.all_foreach_parallel([&](size_t i) {
		auto&c = cman[i];
		c->ca2p._set(ca2p_init);
		c->ca2p_avg = c->ca2p();
		c->ex_inert._set(ex_inert_init);
		c->IP3._set(IP3_init);
		c->connected_cell.foreach([&](Cell* conn) {
			c->gj._set(conn, gj_init);
		});
	});
	cman.other_foreach_parallel([&](size_t i) {
		auto&c = cman[i];
		if (c->state == DEAD) {
			c->ca2p_avg = 0;
			c->ca2p._set(0.0);
		}
	});


}

//////////////////////////////////////////////////////////////////

void cman_init(CellManager& cells, std::string init_data_path) {
	cman_load_from_file(cells, init_data_path);
	_cman_memb_init(cells);
	connect_cell(cells);
	_cman_value_init(cells);
}