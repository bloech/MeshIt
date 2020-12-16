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

#ifndef _FEFLOW_H_
#define _FEFLOW_H_

#include <QtCore/QtCore>

class C_FeFlowEdg
{
public:
	C_FeFlowEdg();
	C_FeFlowEdg(int, int, int);
	~C_FeFlowEdg();

	QList <int> nodes;
	long index;
};

class C_FeFlowTri
{
public:
	C_FeFlowTri();
	C_FeFlowTri(int, int, int, int);
	~C_FeFlowTri();

	QList <int> nodes;
	long index;
};

class C_FeFlow
{
public:
	~C_FeFlow();
	void generateAllTriangles(int numberoftetrahedra, QList <long> tetrahedronlist);
	void generateUndefinedTriangles(int marker, int numberoftriangles, QList <int> trianglemarkerlist, QList <long> trianglelist);
	void generateDefinedTriangles();
	void generateAllEdges(int numberoftetrahedra, QList <long> tetrahedronlist);
	void generateUndefinedEdges(int marker, int numberofedges, QList <int> edgemarkerlist, QList <long> edgelist);
	void generateDefinedEdges();

	QList <C_FeFlowTri> allTriangles;
	QList <C_FeFlowTri> allTrianglesWithoutDuplicates;
	QList <C_FeFlowTri> undefinedTriangles;
	QList <C_FeFlowTri> definedTriangles;
	QList <C_FeFlowEdg> allEdges;
	QList <C_FeFlowEdg> allEdgesWithoutDuplicates;
	QList <C_FeFlowEdg> undefinedEdges;
	QList <C_FeFlowEdg> definedEdges;
};

#endif	// _FEFLOW_H_
