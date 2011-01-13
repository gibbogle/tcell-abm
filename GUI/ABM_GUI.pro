#CONFIG      += uitools
CONFIG += release

INCLUDEPATH += C:\VTK C:\VTK\rendering \
                C:\VTK-src\GUISupport\Qt C:\VTK-src\common C:\VTK-src\rendering \
                C:\VTK-src\graphics C:\VTK-src\filtering C:\VTK-src\IO \
                C:\VTK-src\imaging \
INCLUDEPATH  += c:/qt/qwt-5.2.1/src

FORMS         = ABM_GUI.ui
HEADERS       = mainwindow.h qmylabel.h params.h plot.h log.h myvtk.h misc.h \
    			libpara32.h \
				transfer.h
RESOURCES     += icons.qrc
SOURCES       = main.cpp mainwindow.cpp params.cpp plot.cpp \
                myvtk.cpp misc.cpp lognormal.cpp

LIBS += -L"c:\VTK\bin" -lvtkCommon -lvtkGraphics \
		-lvtkFiltering -lvtkIO -lvtkImaging -lvtkRendering -lQVTK
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
#LIBS += C:\windows\system\libpara32.dll
DEFINES += __MSVS_DLL__
LIBS += C:\windows\system\libpara32_ms.dll

DEFINES += __COMPILETIME_LOADING__

#NOTE: When the DEFINES are changed it is necessary to do Clean then Build
#NOTE: Library symbol info in libpara32.h
