TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    field.cpp \
    cell.cpp

HEADERS += \
    field.h \
    define.h \
    cell.h

QMAKE_CXXFLAGS_RELEASE-=-O2
QMAKE_CXXFLAGS+= -std=c++14 -O3 -ipo
QMAKE_LFLAGS+= -O3 -ipo
