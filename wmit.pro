TARGET = wmit
TEMPLATE = app

QT += opengl xml

INCLUDEPATH += src src/basic src/formats src/ui src/widgets 3rdparty/GLee

HEADERS += \
    src/formats/WZM.hpp \
    src/formats/Pie.hpp \
    src/formats/OBJ.hpp \
    src/formats/Mesh.hpp \
    src/ui/UVEditor.hpp \
    src/ui/TransformDock.hpp \
    src/ui/TeamColoursDock.hpp \
    src/ui/MainWindow.hpp \
    src/ui/ImportDialog.hpp \
    src/ui/ExportDialog.hpp \
    src/Util.hpp \
    src/Generic.hpp \
    src/basic/VectorTypes.hpp \
    src/basic/Vector.hpp \
    src/basic/Polygon.hpp \
    src/basic/IGLTextureManager.hpp \
    src/basic/IGLRenderable.hpp \
    src/basic/IAnimatable.hpp \
    src/basic/GLTexture.hpp \
    3rdparty/GLee/GLee.h \
    src/widgets/QWZM.hpp \
    src/widgets/QtGLView.hpp \
    src/wmit.h \
    src/basic/IGLTexturedRenderable.hpp \
    src/basic/IGLShaderManager.h \
    src/basic/IGLShaderRenderable.h \
    src/ui/TextureDialog.h \
    src/ui/TexConfigDialog.hpp
    
SOURCES += \
    src/formats/WZM.cpp \
    src/formats/Pie_t.cpp \
    src/formats/Pie.cpp \
    src/formats/Mesh.cpp \
    src/ui/UVEditor.cpp \
    src/ui/TransformDock.cpp \
    src/ui/TeamColoursDock.cpp \
    src/ui/MainWindow.cpp \
    src/ui/ImportDialog.cpp \
    src/ui/ExportDialog.cpp \
    src/Util.cpp \
    src/main.cpp \
    src/Generic.cpp \
    src/basic/Polygon_t.cpp \
    src/basic/GLTexture.cpp \
    3rdparty/GLee/GLee.c \
    src/widgets/QWZM.cpp \
    src/widgets/QtGLView.cpp \
    src/ui/TextureDialog.cpp \
    src/ui/TexConfigDialog.cpp
    
FORMS += \
    src/ui/UVEditor.ui \
    src/ui/TransformDock.ui \
    src/ui/TeamColoursDock.ui \
    src/ui/MainWindow.ui \
    src/ui/ImportDialog.ui \
    src/ui/ExportDialog.ui \
    src/ui/TextureDialog.ui \
    src/ui/TexConfigDialog.ui
    
OTHER_FILES += \
    TODO.txt \
    COPYING.txt \
    HACKING.txt \
    COPYING.nongpl \
    data/shaders/pie3.vert \
    data/shaders/pie3.frag \
    data/images/notex.png \
    WMIT.xcodeproj/project.pbxproj \
    Resources/Warzone\\ Model\\ Importer\\ Tool-Info.plist \
    Resources/GenericFramework-Info.plist \
    configs/FetchPrebuilt.sh \
    configs/FetchSource.sh \
    configs/lib3ds-All.xcconfig \
    configs/lib3ds-Debug.xcconfig \
    configs/lib3ds-Release.xcconfig \
    configs/Project-All.xcconfig \
    configs/QGLViewer-All.xcconfig \
    configs/QGLViewer-Debug.xcconfig \
    configs/QGLViewer-Release.xcconfig \
    configs/WMIT-All.xcconfig \
    configs/WMIT-Debug.xcconfig \
    configs/WMIT-Release.xcconfig

CONFIG(debug, debug|release) {
    DEFINES *= DEBUG _DEBUG
} else {
    DEFINES *= NDEBUG
}

# turn off c++0x support on win32 till decent compiler is available out-of-the-box
!win32 {
    QMAKE_CXXFLAGS += -std=c++0x
    DEFINES += CPP0X_AVAILABLE
}

LIBS += -lm

RESOURCES += \
    resources.qrc

# If your system uses different paths for QGLViewer, create a file named
# "config.pri" and override the necessary variables below (with "=").
QGLVIEWER_INCL = /usr/include/QGLViewer
QGLVIEWER_LIBS = -lQGLViewer



UI_DIR = ui
MOC_DIR = moc
OBJECTS_DIR = bin

include("config.pri")

INCLUDEPATH += $$QGLVIEWER_INCL
LIBS += $$QGLVIEWER_LIBS
