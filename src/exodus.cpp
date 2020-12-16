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

#include "exodus.h"

/********** Commons **********/
bool
sortTrianglesByNodes(const C_ExodusTri &d1, const C_ExodusTri &d2)
{
	if (d1.nodes[0] == d2.nodes[0] && d1.nodes[1] == d2.nodes[1] && d1.nodes[2] == d2.nodes[2])
		return true;
	if (d1.nodes[0] == d2.nodes[0] && d1.nodes[1] == d2.nodes[1])
		return d1.nodes[2] < d2.nodes[2];
	if (d1.nodes[0] == d2.nodes[0])
		return d1.nodes[1] < d2.nodes[1];
	return d1.nodes[0] < d2.nodes[0];
}


/********** Class C_ExodusTri **********/
C_ExodusTri::C_ExodusTri()
{
	nodes.append(0);
	nodes.append(0);
	nodes.append(0);
	index_blk = 0;
	index_side = 0;
}

C_ExodusTri::C_ExodusTri(int n1, int n2, int n3, int side, int blk)
{
	nodes.append(n1);
	nodes.append(n2);
	nodes.append(n3);
	std::sort(nodes.begin(), nodes.end());
	index_side = side;
	index_blk = blk;
}

C_ExodusTri::~C_ExodusTri()
{
	nodes.clear();
}

/********** Class C_Exodus **********/
C_Exodus::C_Exodus(int dim, int nodes, int elem, int surfaces, int num_materials, int numberoftetrahedra)
{
	num_dim = dim;
	num_nodes = nodes;
	num_elem = elem;
	num_surfaces = surfaces;
	num_elem_blk = num_materials;
	coord_names[0] = "xcoor";
	coord_names[1] = "ycoor";
	coord_names[2] = "zcoor";

	num_side_ss = new int[num_surfaces];
	x = new float[num_nodes];
	y = new float[num_nodes];
	z = new float[num_nodes];
	node_num_map = new int[num_nodes];
	elem_num_map = new int[num_elem];
	num_nodes_per_elem = new int[num_elem_blk];
	ebids = new int[num_elem_blk];
	connections = new QList<int>[num_elem_blk];
	num_elem_in_blk = new int[num_elem_blk];
	for (int n = 0; n != num_elem_blk; n++)
		num_elem_in_blk[n] = 0;
	elem_ss_map = new int[numberoftetrahedra];
}

void
C_Exodus::deallocate()
{
	delete[] num_nodes_per_elem;
	delete[] ebids;
	delete[] num_side_ss;
	delete[] x;
	delete[] y;
	delete[] z;
	delete[] node_num_map;
	delete[] elem_num_map;
	delete[] connections;
}

void
C_Exodus::generateAllTriangles(int numberoftetrahedra, QList<long> tetrahedronlist)
{
	allTriangles.clear();
	int c1, c2, c3, c4;
	int t1, t2, t3;
	int count_blk = 0;
	for (int b = 0; b != num_elem_blk; b++)
		if (num_nodes_per_elem[b] == 4)
			for (int c = 0; c < connections[b].length(); c += 4)
			{
				c1 = connections[b][c + 0];
				c2 = connections[b][c + 1];
				c3 = connections[b][c + 2];
				c4 = connections[b][c + 3];
				count_blk++;
				allTriangles.append(C_ExodusTri(c1, c2, c3, side(c1, c2, c3, c1, c2, c3, c4), count_blk));
				allTriangles.append(C_ExodusTri(c1, c2, c4, side(c1, c2, c4, c1, c2, c3, c4), count_blk));
				allTriangles.append(C_ExodusTri(c2, c3, c4, side(c2, c3, c4, c1, c2, c3, c4), count_blk));
				allTriangles.append(C_ExodusTri(c3, c1, c4, side(c3, c1, c4, c1, c2, c3, c4), count_blk));
			}
	qSort(allTriangles.begin(), allTriangles.end(), sortTrianglesByNodes);
}

int
C_Exodus::side(int t1, int t2, int t3, int c1, int c2, int c3, int c4)
{
	if (t1 != c1 && t2 != c1 && t3 != c1)
		return 2;
	if (t1 != c2 && t2 != c2 && t3 != c2)
		return 3;
	if (t1 != c3 && t2 != c3 && t3 != c3)
		return 1;
	if (t1 != c4 && t2 != c4 && t3 != c4)
		return 4;
	return -1;
}

int
C_Exodus::generateMarkedTriangles(int marker, int numberoftriangles, QList <int> trianglemarkerlist, QList <long> trianglelist)
{
	int num_side = 0;
	if(!markedTriangles.isEmpty())
		markedTriangles.clear();
	for (int t = 0; t != numberoftriangles; t++)
		if (trianglemarkerlist[t] == marker)
		{
			markedTriangles.append(C_ExodusTri(trianglelist[t * 3 + 0] + 1, trianglelist[t * 3 + 1] + 1, trianglelist[t * 3 + 2] + 1, 0, 0));
			num_side++;
		}
	qSort(markedTriangles.begin(), markedTriangles.end(), sortTrianglesByNodes);
	return num_side;
}

int
C_Exodus::findTriangle(C_ExodusTri tri, int start)
{
	for (int t = start; t != allTriangles.length(); t++)
		if ((allTriangles[t].nodes[0] == tri.nodes[0] || allTriangles[t].nodes[0] == tri.nodes[1] || allTriangles[t].nodes[0] == tri.nodes[2]) &&
			(allTriangles[t].nodes[1] == tri.nodes[0] || allTriangles[t].nodes[1] == tri.nodes[1] || allTriangles[t].nodes[1] == tri.nodes[2]) &&
			(allTriangles[t].nodes[2] == tri.nodes[0] || allTriangles[t].nodes[2] == tri.nodes[1] || allTriangles[t].nodes[2] == tri.nodes[2]))
			return t;
	return -1;
}
