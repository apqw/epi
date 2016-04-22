TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    epi/define.cpp

HEADERS += \
    epi/define.h

QMAKE_CXXFLAGS += -mavx
