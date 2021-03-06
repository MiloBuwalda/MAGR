# File generated by kdevelop's qmake manager. 
# ------------------------------------------- 
# Subdir relative project main directory: ./../viewers/widgets
# Target is a library:  

FORMS += GLGeometryViewerDesigner.ui \
         LogoViewerDesigner.ui \
         GLGeometryViewer3DDesigner.ui 
HEADERS += GeoXGLWidget.h \
           GLGeometryViewer.h \
           LogoViewer.h \
           VisTypes.h \
           GeoXGLWidget3D.h \
           GLGeometryViewer3D.h 
SOURCES += GeoXGLWidget.cpp \
           GLGeometryViewer.cpp \
           LogoViewer.cpp \
           GeoXGLWidget3D.cpp \
           GLGeometryViewer3D.cpp 
INCLUDEPATH += ../../windows \
../../viewers/widgets \
../../viewers/cameraControl \
../../system/basics \
../../system/streaming \
../../system/properties \
../../system/misc \
../../system/gui \
../../math \
$QTDIR/include/
CONFIG += debug \
warn_on \
qt \
staticlib
TEMPLATE = lib
