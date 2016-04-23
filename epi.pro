TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    epi/define.cpp \
    epi/main.cpp \
    VCL/dispatch_example.cpp \
    VCL/instrset_detect.cpp \

HEADERS += \
    epi/define.h \
    VCL/instrset.h \
    VCL/vectorclass.h \
    VCL/vectorf128.h \
    VCL/vectorf256.h \
    VCL/vectorf256e.h \
    VCL/vectorf512.h \
    VCL/vectorf512e.h \
    VCL/vectori128.h \
    VCL/vectori256.h \
    VCL/vectori256e.h \
    VCL/vectori512.h \
    VCL/vectori512e.h \
    VCL/vectormath_common.h \
    VCL/vectormath_exp.h \
    VCL/vectormath_hyp.h \
    VCL/vectormath_lib.h \
    VCL/vectormath_trig.h

QMAKE_CXXFLAGS += -mavx -fabi-version=0

LIBS+= -lsvml -limf
