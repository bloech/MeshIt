# MeshIt - a 3D mesh generator for fractured reservoirs
#
# Copyright (C) 2020
#
# Mauro Cacace (GFZ, cacace@gfz-potsdam.de),
# Guido Blöcher (GFZ, bloech@gfz-potsdam.de),
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version, complemented with
# the following provision:
# For the scientific transparency and verification of results obtained
# and communicated to the public after using a modified version of the
# work, You (as the recipient of the source code and author of this
# modified version, used to produce the published results in scientific
# communications) commit to make this modified source code available in
# a repository that is easily and freely accessible for a duration of
# five years after the communication of the obtained results.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

TEMPLATE = app
TARGET = MeshIt
ICON = ./resources/images/app_logo.icns
CONFIG += console
CONFIG += warn_off
QT += widgets opengl
INCLUDEPATH += ./include

# Include system-dependent variables.
!include(meshit.pri) {
    EXODUS_LIBMESH = false
    EXODUS_LIBRARY = false
}

# Mac
macx {
    QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.15
}

# Linux
unix:!macx {
    LIBS += -lGLU
}

# Windows - MinGW
win32-g++ {
    LIBS += libopengl32 libglu32
}

# Windows - Microsoft Visual C++
win32-msvc* {
    LIBS += opengl32.lib glu32.lib
}

if($$EXODUS_LIBMESH) {
    INCLUDEPATH += . $$LIBMESH/include
    INCLUDEPATH += . $$LIBMESH/include/libmesh
    INCLUDEPATH += . $$LIBMESH/../include
    LIBS += -L$$LIBMESH/lib -lmesh_opt
    QMAKE_LFLAGS += -Wl,-rpath,$$LIBMESH
}else:if($$EXODUS_LIBRARY) {
    INCLUDEPATH += . $$EXODUS_PATH/cbind/include
    LIBS += $$EXODUS_PATH/build/cbind/Release/exoIIv2c.lib
    INCLUDEPATH += . $$NETCDF_PATH/include
    LIBS += -L$$NETCDF_PATH/lib -lnetcdf
}else{
    DEFINES += NOEXODUS
}

# Configuration of Triangle library.
DEFINES += TRILIBRARY EXTERNAL_TEST

HEADERS += include/c_vector.h \
           include/geometry.h \
           include/glwidget.h \
           include/commandline.h \
           include/intersections.h \
           include/mainwindow.h \
           include/tetgen.h \
           include/triangle.h \
           include/feflow.h \
           include/exodus.h \
           include/core.h
SOURCES += src/geometry.cpp \
           src/glwidget.cpp \
           src/commandline.cpp \
           src/main.cpp \
           src/mainwindow.cpp \
           src/predicates.cxx \
           src/tetgen.cxx \
           src/triangle.c \
           src/feflow.cpp \
           src/exodus.cpp \
           src/core.cpp
RESOURCES += resources/MeshIT.qrc
