#include <iostream>
#include "codetest.h"
using namespace std;

int main(int argc, char *argv[])
{
    cout << "Hello World!" << endl;
    atomic_double_test();
    lfpstack_test();
#ifdef _WIN32
	system("pause");
#endif
    return 0;
}
