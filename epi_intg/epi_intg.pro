TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    calc/cpu2/ca2p.cpp \
    calc/cpu2/cell.cpp \
    calc/cpu2/cell_connection.cpp \
    calc/cpu2/cell_interaction.cpp \
    calc/cpu2/cell_state_renew.cpp \
    calc/cpu2/CellManager.cpp \
    calc/cpu2/ext_stim.cpp \
    calc/cpu2/map.cpp \
    calc/calc.cpp \
    calc/CalcParams.cpp \
    calc/init_gen.cpp \
    lodepng/lodepng.cpp \
    util/rand/Random_gen.cpp \
    util/vec/Vec.cpp \
    vis/color.cpp \
    vis/VisParams.cpp \
    vis/visualize.cpp \
    fsys.cpp \
    global.cpp \
    Params.cpp \
    parser.cpp \
    utils.cpp

DISTFILES += \
    cmdline/LICENSE

HEADERS += \
    calc/cpu2/atomics.h \
    calc/cpu2/ca2p.h \
    calc/cpu2/cell.h \
    calc/cpu2/cell_connection.h \
    calc/cpu2/cell_interaction.h \
    calc/cpu2/cell_state_renew.h \
    calc/cpu2/CellManager.h \
    calc/cpu2/DualValue.h \
    calc/cpu2/ext_stim.h \
    calc/cpu2/map.h \
    calc/calc.h \
    calc/CalcParams.h \
    calc/init_gen.h \
    cmdline/cmdline.h \
    lodepng/lodepng.h \
    misc/DynArr.h \
    misc/swapdata.h \
    util/rand/Random_gen.h \
    util/vec/Vec.h \
    vis/color.h \
    vis/vis_global.h \
    vis/VisParams.h \
    vis/visualize.h \
    define.h \
    fsys.h \
    global.h \
    Params.h \
    parser.h \
    utils.h
