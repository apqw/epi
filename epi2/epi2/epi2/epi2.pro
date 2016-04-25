TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    update.cpp \
    SFMT.c

HEADERS += \
    define.h \
    update.h \
    SFMT-params.h \
    SFMT-params607.h \
    SFMT.h \
    SFMT-params19937.h

QMAKE_CXXFLAGS += -fopenmp -mavx
QMAKE_LFLAGS += -fopenmp -DMEXP=607
LIBS+= ../epi2/SFMT.c
