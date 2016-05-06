
#include <memory>
#include <iostream>
#ifdef _WIN32
#include <Windows.h>
#else
#include "util_func.h"
void OutputDebugString(const char* str) {
	std::cout << "No input specified." << std::endl;
}
#endif

#include "codetest.h"
#include "fsys.h"
int main(int argc, char *argv[]) {
   // system("pause");
	Vec3Test();
	arr_map_test();
	cell_test();
	conv_state_num_test();
	cell_man_test();
	//CAS_test();

	//static_error_check();
	using std::cout;
	using std::endl;
	//using std::cerr;
	if (argc <= 1) {
		OutputDebugString("No input specified.");
		//printf("No input specified.\n");
		exit(1);
	}

	make_dir(OUTPUTDIR);

	std::ifstream data;

	data.open(argv[1], std::ios::in);//read-only
	if (!data) {
		cout << "Couldn't open the file " << argv[1] << " ." << endl;
		exit(1);
	}

	auto fld = std::make_unique<Field>();
	fld->init_with_file(data);
genrand_init();
	fld->main_loop();
	system("pause");
}
