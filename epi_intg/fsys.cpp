
#include "fsys.h"
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif
void make_dir(const char* path) {
#ifdef _WIN32
    _mkdir(path);
#else
    mkdir(path,0777);
#endif
}
