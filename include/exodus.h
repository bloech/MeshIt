/* MeshIt - a 3D mesh generator for fractured reservoirs
 *
 * Copyright (C) 2020
 *
 * Mauro Cacace (GFZ, cacace@gfz-potsdam.de),
 * Guido Blöcher (GFZ, bloech@gfz-potsdam.de),
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version, complemented with
 * the following provision:
 * For the scientific transparency and verification of results obtained
 * and communicated to the public after using a modified version of the
 * work, You (as the recipient of the source code and author of this
 * modified version, used to produce the published results in scientific
 * communications) commit to make this modified source code available in
 * a repository that is easily and freely accessible for a duration of
 * five years after the communication of the obtained results.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _EXODUS_H_
#define _EXODUS_H_

#include <QtCore/QtCore>
#include <iostream>

class C_ExodusTri
{
public:
	C_ExodusTri();
	C_ExodusTri(int, int, int, int, int);
	~C_ExodusTri();

	QList <int> nodes;
	int index_side;
	int index_blk;
};

class C_Exodus
{
public:
	C_Exodus(int, int, int, int, int, int);
	~C_Exodus(){};
	void deallocate();
	void generateAllTriangles(int, QList<long>);
	int side(int, int, int, int, int, int, int);
	int generateMarkedTriangles(int, int, QList<int>, QList <long>);
	int findTriangle(C_ExodusTri, int);

	QList <C_ExodusTri> allTriangles;
	QList <C_ExodusTri> markedTriangles;
	int num_dim;
	int num_nodes;
	int num_elem;
	int num_surfaces;
	int num_elem_blk;
	int num_side_sets;
	int num_node_sets;
	int *num_elem_in_blk;
	int *ebids;
	int *num_nodes_per_elem;
	int *node_num_map;
	int *elem_num_map;
	int *elem_ss_map;
	int *elem_list;
	int *side_list;
	int *num_side_ss;
	QList<int> *connections;
    double *x;
    double *y;
    double *z;
	const char *coord_names[3];
};

#endif	// _EXODUS_H_
