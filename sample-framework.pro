# -*- mode: Makefile -*-

TARGET = sample-framework
CONFIG *= qt opengl debug
QT += gui widgets opengl xml

HEADERS = MyWindow.h MyViewer.h MyViewer.hpp
SOURCES = MyWindow.cpp MyViewer.cpp main.cpp vector.cc bspline.cc rmf.cc

LIBS *= -lQGLViewer -L/usr/lib/OpenMesh -lOpenMeshCored -lGL -lGLU

RESOURCES = sample-framework.qrc

QMAKE_CXXFLAGS += -std=c++11
