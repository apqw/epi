TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    21-der.c \
    auto21-1.c \
    ca_dynamics.c \
    cell_dynamics.c \
    func21-2.c \
    SFMT.c

HEADERS += \
    funcs.h \
    main.h \
    param.h \
    SFMT-params.h \
    SFMT-params607.h \
    SFMT.h
QMAKE_CFLAGS += -DMEXP=607
QMAKE_LFLAGS += -DMEXP=607
