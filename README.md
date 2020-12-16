MeshIt 2010-2020
================

The tool MeshIT uses TRIANGLE <http://www.cs.cmu.edu/~quake/triangle.html> and
TETGEN <http://wias-berlin.de/software/tetgen> to generate a quality
tetrahedral mesh based on structural geological information.

This procedure is fully automatized and needs at least scattered data points
as input!

Main developers: Mauro Cacace (<mailto:cacace@gfz-potsdam.de>) and
                 Guido Blöcher (<mailto:bloech@gfz-potsdam.de>).

Some extensions were added by PERFACCT (www.perfacct.eu) by the following developers:
				Johannes Spazier
				Nihed Boussaidi
				Danny Puhan

The source can be compiled with some limitations on Windows, Linux and Apple.
The limitation is the export of exodus file format used by Moose (https://github.com/idaholab/moose) which is only
available for WIN64 (netcdf (https://www.unidata.ucar.edu/software/netcdf/) and exodus libraries (https://github.com/gsjaardema/seacas) are required). On all other
platforms the following requirements are suggested and tested:

Windows:

    -  Qt 5.9.9 for Windows 64-bit (VS 2015, VS 2017)

Linux:

    -  Qt 5.9.9 for Linux 64-bit

Mac:

    -  Qt 5.14.1 for OS X

Preparation:
*	Create a copy of file `meshit.pri.tmpl` and save it as `meshit.pri`.
	Set either the EXODUS_LIBMESH variable or the EXODUS_LIBRARY variable
	to `true` if you want to enable the exodus export option.

*	In case the exodus export (EXODUS_LIBMESH) was enabled, define the path
	to the root directory of a `libmesh` installation using variable LIBMESH
	(a compatible `libmesh` library is usually already installed as part of the Moose framework).

*	In case the exodus export (EXODUS_LIBRARY) was enabled, define the path
	to the root directory of a `netcdf` installation and 'exodusII' installation.

Building:
*  `qmake meshit.pro`
*  `make/nmake/mingw32-make`
