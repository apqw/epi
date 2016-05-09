TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    Cell.cpp \
    codetest.cpp \
    Field.cpp \
    general_func.cpp \
    primitive_func.cpp \
    util_func.cpp

HEADERS += \
    Cell.h \
    codetest.h \
    component.h \
    define.h \
    Field.h \
    fsys.h \
    general_func.h \
    primitive_func.h \
    util_func.h
QMAKE_INCDIR+=/opt/intel/tbb44_20160413oss/include
QMAKE_CXXFLAGS_RELEASE-=-O2
QMAKE_CXXFLAGS+=-std=c++14 -Ofast -flto -march=haswell
QMAKE_LFLAGS+=-Ofast -flto -march=haswell -L/opt/intel/tbb44_20160413oss/intel64/gcc4.4 -L/home/yasu7890v/opt/intel/lib/intel64 -ltbb
