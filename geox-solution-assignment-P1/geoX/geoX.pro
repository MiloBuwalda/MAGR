# File generated by kdevelop's qmake manager. 
# ------------------------------------------- 
# Subdir relative project main directory: ./../geoX
# Target is a library:  

FORMS += GeoXMainWindowDesigner.ui

INCLUDEPATH += ../experiments \
  ../viewers/widgets \
  ../viewers/cameraControl \
  ../geoX \
  ../math \
  ../windows \
  ../system/basics \
  ../system/streaming \
  ../system/properties \
  ../system/misc \
  ../system/gui \
               $QTDIR/include/ 


CONFIG += debug \
          warn_on \
          qt \
          staticlib 
TEMPLATE = lib 

HEADERS += GeoXMainWindow.h
SOURCES += GeoXMainWindow.cpp

LIBS += ../system/basics/libbasics.a
