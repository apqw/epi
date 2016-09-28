TEMPLATE = app
CONFIG += console
#CONFIG -=c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    cell.cpp \
    codetest.cpp \
    cell_connection.cpp \
    cell_init.cpp \
    cell_interaction.cpp \
    cell_state_renew.cpp \
    cellmanager.cpp \
    DualValue.cpp \
    ext_stim.cpp \
    Field.cpp \
    map.cpp \
    proc.cpp \
    Random_gen.cpp \
    utils.cpp \
    ca2p.cpp \
    fsys.cpp \
    vec3.cpp

HEADERS += \
    cell.h \
    atomics.h \
    codetest.h \
    define.h \
    swapdata.h \
    cell_connection.h \
    cell_init.h \
    cell_interaction.h \
    cell_state_renew.h \
    cellmanager.h \
    DualValue.h \
    ext_stim.h \
    Field.h \
    map.h \
    proc.h \
    Random_gen.h \
    utils.h \
    cell_conn_value.h \
    ca2p.h \
    fsys.h \
    vec3.h
QMAKE_INCDIR+=/opt/intel/tbb44_20160526oss/include
QMAKE_CXXFLAGS_RELEASE-=-O2
#QMAKE_CXXFLAGS+=-std=c++11 -O3 -ipo -xHost -no-prec-div
#QMAKE_LFLAGS+=-std=c++11 -O3 -ipo -xHost -no-prec-div -L/opt/intel/tbb44_20160526oss/intel64/gcc4.4  -L/home/yasu7890v/opt/intel/lib/intel64 -ltbb -ltbbmalloc
QMAKE_CXXFLAGS+=-O3
QMAKE_LFLAGS+=-Ofast -L/opt/intel/tbb44_20160526oss/intel64/gcc4.4 -L/home/yasu7890v/opt/intel/lib/intel64 -ltbb -ltbbmalloc
