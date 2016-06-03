TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    cell.cpp \
    codetest.cpp

HEADERS += \
    cell.h \
    atomics.h \
    codetest.h \
    define.h \
    swapdata.h
QMAKE_INCDIR+=/opt/intel/tbb44_20160413oss/include
QMAKE_CXXFLAGS_RELEASE-=-O2
QMAKE_CXXFLAGS+=-std=c++14 -O3 -ipo -xCORE-AVX2
QMAKE_LFLAGS+=-O3  -ipo -xCORE-AVX2 -L/opt/intel/tbb44_20160413oss/intel64/gcc4.4 -L/home/yasu7890v/opt/intel/lib/intel64 -ltbb
