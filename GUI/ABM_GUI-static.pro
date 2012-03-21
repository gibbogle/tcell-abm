#CONFIG      += uitools
CONFIG += release

INCLUDEPATH +=  C:\VTK-build-4.4.0-static C:\VTK-build-4.4.0-static\rendering \
                C:\VTK-src\GUISupport\Qt C:\VTK-src\common C:\VTK-src\rendering \
                C:\VTK-src\graphics C:\VTK-src\filtering C:\VTK-src\IO \
                C:\VTK-src\imaging \
INCLUDEPATH  += c:/qt/qwt-5.2.1/src

FORMS         = ABM_GUI.ui
HEADERS       = mainwindow.h qmylabel.h params.h plot.h log.h myvtk.h misc.h \
				libpara32.h result_set.h graphs.h
RESOURCES     += icons.qrc
SOURCES       = main.cpp mainwindow.cpp params.cpp plot.cpp \
                myvtk.cpp misc.cpp lognormal.cpp \
    graphs.cpp

# See cmake_link_command.txt for the full list of libraries that CMake links
LIBS += -L"c:\VTK-build-4.4.0-static\bin" -lQVTK -lvtkRendering -lvtkGraphics -lvtkImaging -lvtkIO -lvtkFiltering -lvtkCommon \
-lvtkpng -lvtktiff -lvtkjpeg -lvtkexpat -lvfw32 -lopengl32 -lVPIC -lCosmo \
-lwsock32 -lvtksys -lws2_32 -lvtkexoIIc -lvtkNetCDF \
-lvtklibxml2 -lvtkzlib -lpthread -lvtkalglib \
-lgdi32 -lkernel32 -luser32 -lgdi32 -lwinspool -lshell32 -lole32 -loleaut32 -luuid -lcomdlg32 -ladvapi32

#LIBS += -LC:\users\gib\abm\build32msvs\release -lpara32
LIBS += -LC:\qt\qwt-5.2.1\lib -lqwt5

DEPENDPATH   += c:/qt/qwt-5.2.1/lib

QT           += network

QMAKE_CXXFLAGS += -Wno-deprecated -Wno-write-strings

# install
target.path = .
sources.files = $$SOURCES $$HEADERS $$RESOURCES $$FORMS ABM_GUI.pro icons
sources.path = .
INSTALLS += target sources

#DEFINES += __GFORTRAN_DLL__
DEFINES += __MSVS_DLL__
LIBS += C:\windows\system\libpara32_ms.dll

DEFINES += __COMPILETIME_LOADING__
# Note: On Windows compile-time DLL loading apparently does not permit OpenMP to use multiple threads
