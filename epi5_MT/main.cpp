#include <iostream>
#include "codetest.h"
#include "proc.h"
#include <string>
using namespace std;

int main(int argc, char *argv[])
{
    cout << argv[1] << endl;
    atomic_double_test();
    lfpstack_test();
	cell_test();
	proc(argv[1], argv[2]);
#ifdef _WIN32
	system("pause");
#endif
    return 0;
}
