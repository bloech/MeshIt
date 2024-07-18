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

The source can be compiled as it comes on Windows, Linux and MacOS by running the building options below.
For exporting the exodus file format used by Moose (https://github.com/idaholab/moose) the user has to specify some internal flags in the project file meshit.pro:
 
*	For Linux and Mac we suggest to link the static exodus library which comes along with libMesh (https://libmesh.github.io/) provided by Moose framework installation:
	*	Set the EXODUS_LIBMESH variable to "true"
	*	Define the path to the root directory of the `libmesh` installation using variable LIBMESH

*	For Windows you have to link the dynamic exodus library which will be provided by the package mingw-w64-ucrt-x86_64-libexodus provided by the MSYS2 (https://www.msys2.org/) installation:
	*	Set the EXODUS_LIBRARY variable to "true"
	*	Define the path the rootdirectory of 'exodusII' installation
	
On all platforms the following requirements are suggested and tested:

Windows:

    -  Qt 5.15.2 for Windows 64-bit

Linux:

    -  Qt 5.9.9 for Linux 64-bit

Mac:

    -  Qt 5.14.1 for OS X

Preparations:
*	For Windows users, please follow the following steps for a proper "exodusII" dynamic library installation:
	*	Download the "msys2-x86_64" installer (https://www.msys2.org/)
	*	Run the installer. Installing MSYS2 requires 64 bit Windows 10 or newer.
	*	Enter your desired Installation Folder (short ASCII-only path on a NTFS volume, no accents, no spaces, no symlinks, no subst or network drives, no FAT).
	*	When done, click Finish.
	*	Now MSYS2 is ready for you and a terminal for the UCRT64 environment will launch.
	*	Install "mingw-w64-ucrt-x86_64-libexodus". Run the following command:
	
	pacman -S mingw-w64-ucrt-x86_64-libexodus

	*	To enable the dynamic dependencies of the "exodusII" library add the root folder to your environment variable path:
		*	Open the Start Search, type in “env”, and choose “Edit the system environment variables”.
		*	Click the “Environment Variables…” button.
		*	Under the “User Variables” section, find the row with “Path” in the first column, and click edit.
		*	The “Edit environment variable” UI will appear. Here, you can click “New” and type in the new path (default installation path C:\msys64\ucrt64) you want to add.
		*	Dismiss all of the dialogs by choosing “OK”. Your changes are saved!
		*	You will probably need to restart apps for them to pick up the change. Restarting the machine would ensure all apps are run with the PATH change.

*	For macOS users, please check the version of macOS specified in the
	`meshit.pro` file (line `QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.15`). By default, the version of macOS specified is 10.15 (macOS Catalina)

Building:
*	`qmake meshit.pro`
*	`make/nmake/mingw32-make`

Cite:
*	https://doi.org/10.5281/zenodo.4327281
*	Cacace, M., Blöcher, G. MeshIt - a software for three dimensional volumetric meshing of complex faulted reservoirs. Environ Earth Sci 74, 5191–5209 (2015). https://doi.org/10.1007/s12665-015-4537-x
