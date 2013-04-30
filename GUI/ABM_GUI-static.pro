CONFIG += release

INCLUDEPATH += C:/VTK-VS8/include/vtk-5.10
INCLUDEPATH += c:/qwt-5.2.1/src

FORMS         = ABM_GUI.ui SimpleView2DUI.ui SimpleView3DUI.ui
HEADERS       = mainwindow.h qmylabel.h params.h plot.h log.h myvtk.h misc.h \
                libpara32.h result_set.h graphs.h \
                SimpleView2DUI.h SimpleView3DUI.h ImageSave.h
RESOURCES     += icons.qrc
SOURCES       = main.cpp mainwindow.cpp params.cpp plot.cpp \
                myvtk.cpp misc.cpp lognormal.cpp graphs.cpp \
                SimpleView2DUI.cxx SimpleView3DUI.cxx ImageSave.cpp

# See cmake_link_command.txt for the full list of libraries that CMake links
LIBS += -LC:/VTK-VS8/lib/vtk-5.10 -lQVTK -lvtkRendering -lvtkGraphics -lvtkImaging -lvtkIO -lvtkFiltering -lvtkCommon \
-lvtkpng -lvtktiff -lvtkjpeg -lvtkexpat -lvfw32 -lopengl32 -lwsock32 -lvtksys -lws2_32 -lvtkexoIIc -lvtkNetCDF \
-lvtklibxml2 -lvtkzlib -lvtkalglib -lgdi32 -lkernel32 -luser32 -lgdi32 -lwinspool -lshell32 -lole32 -loleaut32 \
-luuid -lcomdlg32 -ladvapi32
# -lVPIC -lCosmo -lpthread

#LIBS += -LC:\users\gib\abm\build32msvs\release -lpara32
LIBS += -LC:/bin -lqwt5
LIBS += C:/bin/para-v2.lib

QT           += network

#QMAKE_CXXFLAGS += -Wno-deprecated -Wno-write-strings

# install
target.path = .
sources.files = $$SOURCES $$HEADERS $$RESOURCES $$FORMS ABM_GUI.pro icons
sources.path = .
INSTALLS += target sources

DEFINES += _CRT_SECURE_NO_WARNINGS
