#include "data_IO.h"
#include "cell.h"
#include "cellmanager.h"
#include <iostream>
#include <fstream>
void cell_bin_output(std::string path)
{
	std::ofstream fout;
	fout.open(path, std::ios::out | std::ios::binary);
	assert(fout);

}
