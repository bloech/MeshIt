# This file is a template for the definition of system-dependend variables.
# Please create a copy of this file under the name "meshit.pri" and fill in
# the missing variables.

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

# Set either EXODUS_LIBMESH or EXODUS_LIBRARY to true if you want the exodus export option.
EXODUS_LIBMESH = true
EXODUS_LIBRARY = false

# If EXODUS_LIBMESH is true, insert below the path to the root directory of a libmesh installation.
LIBMESH= /home/mauro/projects/moose_libmesh/libmesh/installed

# If EXODUS_LIBRARY is true, insert below the path to the root directory of the exodusII library and netcdf.
#NETCDF_PATH = <path-to-netcdf>
#EXODUS_PATH = <path-to-exodusII>
