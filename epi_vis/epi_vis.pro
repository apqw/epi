TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

QMAKE_INCDIR += include
LIBS += -L../epi_vis/lib/ -lvoro++
