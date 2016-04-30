#define NDEBUG
#include <iostream>
#include <fstream>
#include <memory>
#include "fsys.h"
#include "const.h"
#include "cell.h"
#include "field.h"
#ifdef _WIN32
#include <Windows.h>
#else
void OutputDebugString(const char* str) {
	std::cout << "No input specified." << std::endl;
}
#endif

int main(int argc, char *argv[]) {
	static_error_check();
	using std::cout;
	using std::endl;
	//using std::cerr;
	if (argc <= 1) {
		OutputDebugString("No input specified.");
		exit(1);
	}

	make_dir(OUTPUTDIR);

	std::ifstream data;

	data.open(argv[1], std::ios::in);//read-only
	if (!data) {
		cout << "Couldn't open the file " << argv[1] << " ." << endl;
		exit(1);
	}

	auto field = std::unique_ptr<Field>(new Field());
	cout << "Loading file..." << endl;
	field->init_with_file(data);
	cout << "done." << endl;

	field->main_loop(null_coef());
	system("pause");
}