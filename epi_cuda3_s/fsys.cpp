#include "fsys.h"
#include <cstdlib>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

void pause(){
#ifdef _WIN32
	system("pause");
#else
	system("read -p \"Press Enter to proceed...\"");
#endif
}

void make_dir(const char* dir){
#ifdef _WIN32
	_mkdir(dir);
#else
	mkdir(dir, 0777);
#endif
}

void exit_after_pause(){
	pause();
	exit(1);
}