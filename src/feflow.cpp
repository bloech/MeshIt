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

#include "feflow.h"

/********** Commons **********/
bool
sortEdgesByNodes(const C_FeFlowEdg &d1, const C_FeFlowEdg &d2)
{
	if (d1.nodes[0] == d2.nodes[0])
		return d1.nodes[1] < d2.nodes[1];
	return d1.nodes[0] < d2.nodes[0];
}

bool
sortEdgesByIndex(const C_FeFlowEdg &d1, const C_FeFlowEdg &d2)
{
	return d1.index < d2.index;
}

bool
sortTrianglesByNodes(const C_FeFlowTri &d1, const C_FeFlowTri &d2)
{
	if (d1.nodes[0] == d2.nodes[0] && d1.nodes[1] == d2.nodes[1] && d1.nodes[2] == d2.nodes[2])
		return d1.index < d2.index;
	if (d1.nodes[0] == d2.nodes[0] && d1.nodes[1] == d2.nodes[1])
		return d1.nodes[2] < d2.nodes[2];
	if (d1.nodes[0] == d2.nodes[0])
		return d1.nodes[1] < d2.nodes[1];
	return d1.nodes[0] < d2.nodes[0];
}

/********** Class C_FeflowEdge **********/
C_FeFlowEdg::C_FeFlowEdg()
{
	this->nodes.append(0);
	this->nodes.append(0);
	this->index = 0;
}

C_FeFlowEdg::C_FeFlowEdg(int n1, int n2, int index)
{
	this->nodes.append(n1);
	this->nodes.append(n2);
	std::sort(nodes.begin(), nodes.end());
	this->index = index;
}

C_FeFlowEdg::~C_FeFlowEdg()
{
	this->nodes.clear();
}

/********** Class C_FeflowTri **********/
C_FeFlowTri::C_FeFlowTri()
{
	this->nodes.append(0);
	this->nodes.append(0);
	this->nodes.append(0);
	this->index = 0;
}

C_FeFlowTri::C_FeFlowTri(int n1, int n2, int n3, int index)
{
	this->nodes.append(n1);
	this->nodes.append(n2);
	this->nodes.append(n3);
	std::sort(nodes.begin(), nodes.end());
	this->index = index;
}

C_FeFlowTri::~C_FeFlowTri()
{
	this->nodes.clear();
}

bool sortTrianglesByIndex(const C_FeFlowTri &d1, const C_FeFlowTri &d2)
{
	return d1.index < d2.index;
}

/********** Class C_Feflow **********/
C_FeFlow::~C_FeFlow()
{
	this->allEdges.clear();
	this->allEdgesWithoutDuplicates.clear();
	this->allTriangles.clear();
	this->allTrianglesWithoutDuplicates.clear();
	this->definedEdges.clear();
	this->definedTriangles.clear();
	this->undefinedEdges.clear();
	this->undefinedTriangles.clear();
}

void
C_FeFlow::generateAllTriangles(int numberoftetrahedra, QList <long> tetrahedronlist)
{
	this->allTriangles.clear();
	// Insert all possible triangles based on all tetrahedrons to the List allTriangles. The constructor C_FeFlowTri(n1,n2,n3,index) will sort the nodeIDs n1<=n2<=n3
	for (int t = 0; t != numberoftetrahedra; t++)
	{
		allTriangles.append(C_FeFlowTri(tetrahedronlist[t * 4 + 0], tetrahedronlist[t * 4 + 1], tetrahedronlist[t * 4 + 2], this->allTriangles.length()));
		allTriangles.append(C_FeFlowTri(tetrahedronlist[t * 4 + 0], tetrahedronlist[t * 4 + 1], tetrahedronlist[t * 4 + 3], this->allTriangles.length()));
		allTriangles.append(C_FeFlowTri(tetrahedronlist[t * 4 + 1], tetrahedronlist[t * 4 + 2], tetrahedronlist[t * 4 + 3], this->allTriangles.length()));
		allTriangles.append(C_FeFlowTri(tetrahedronlist[t * 4 + 2], tetrahedronlist[t * 4 + 0], tetrahedronlist[t * 4 + 3], this->allTriangles.length()));
	}
	// Sort allTriangles first by n1 than by n2 then by n3 then by their index. The index is the previous position in the list, which gives the general order of triangles
	qSort(allTriangles.begin(), allTriangles.end(), sortTrianglesByNodes);
	this->allTrianglesWithoutDuplicates.clear();
	// Transfer all induvidual triangles to alltrianglesWithoutDublicates. Ff triangles are duplicated, the one with the lowest index will be transfered only
	for (int t = 0; t != this->allTriangles.length(); t++)
		if (allTrianglesWithoutDuplicates.length() == 0 ||
			allTriangles[t].nodes[0] != allTrianglesWithoutDuplicates.last().nodes[0] ||
			allTriangles[t].nodes[1] != allTrianglesWithoutDuplicates.last().nodes[1] ||
			allTriangles[t].nodes[2] != allTrianglesWithoutDuplicates.last().nodes[2])
			allTrianglesWithoutDuplicates.append(allTriangles[t]);
	// Resort all triangles by index, the previous order is recaptured. This means no ascending list with removed duplicates
	qSort(allTrianglesWithoutDuplicates.begin(), allTrianglesWithoutDuplicates.end(), sortTrianglesByIndex);
	// Reindex the right ordered list (based on the order). The new index is now the feflow index
	for (int t = 0; t != this->allTrianglesWithoutDuplicates.length(); t++)
		this->allTrianglesWithoutDuplicates[t].index = t;
	// Sort the list without duplicates (nodes and index are ascending --> this is required the find a triangle faster)
	qSort(allTrianglesWithoutDuplicates.begin(), allTrianglesWithoutDuplicates.end(), sortTrianglesByNodes);
}

void
C_FeFlow::generateUndefinedTriangles(int marker, int numberoftriangles, QList <int> trianglemarkerlist, QList <long> trianglelist)
{
	this->undefinedTriangles.clear();
	// Store all triangles of tetgen having m as marker m (undefinedTriangles list). The index of all triangles are set to ZERO.
	for (int t = 0; t != numberoftriangles; t++)
		if (trianglemarkerlist[t] == marker)
			undefinedTriangles.append(C_FeFlowTri(trianglelist[t * 3 + 0], trianglelist[t * 3 + 1], trianglelist[t * 3 + 2], 0));
	// Sort undefinedTriangles by nodes and index. This supports the search of each induvidual triangle in the allTrianglesWihoutDuplicates list.
	qSort(undefinedTriangles.begin(), undefinedTriangles.end(), sortTrianglesByNodes);
}

void
C_FeFlow::generateDefinedTriangles()
{
	//Search for each undefined triangle in allTrianglesWithoutDuplicates (this list has the required feflow index). Found -> the triangle with the feflow index is copied to the definedTriangles list
	this->definedTriangles.clear();
	int lastPos = 0;
	for (int u = 0; u != this->undefinedTriangles.length(); u++)
		for (int a = lastPos; a != this->allTrianglesWithoutDuplicates.length(); a++)
			if (this->undefinedTriangles[u].nodes[0] == this->allTrianglesWithoutDuplicates[a].nodes[0] &&
				this->undefinedTriangles[u].nodes[1] == this->allTrianglesWithoutDuplicates[a].nodes[1] &&
				this->undefinedTriangles[u].nodes[2] == this->allTrianglesWithoutDuplicates[a].nodes[2]){
				this->definedTriangles.append(this->allTrianglesWithoutDuplicates[a]);
				lastPos = a;
				break;
			}
	// The list with the defined triangles is sorted by the index (to have a not required but sorted output)
	qSort(definedTriangles.begin(), definedTriangles.end(), sortTrianglesByIndex);
}

void
C_FeFlow::generateAllEdges(int numberoftetrahedra, QList <long> tetrahedronlist)
{
	this->allEdges.clear();
	// Insert all possible edges based on all tetrahedrons to the List allTriangles (C_FeFlowTri(n1,n2,index) will sort the nodeIDs n1<=n2)
	for (int t = 0; t != numberoftetrahedra; t++)
	{
		allEdges.append(C_FeFlowEdg(tetrahedronlist[t * 4 + 0], tetrahedronlist[t * 4 + 1], this->allEdges.length()));
		allEdges.append(C_FeFlowEdg(tetrahedronlist[t * 4 + 0], tetrahedronlist[t * 4 + 2], this->allEdges.length()));
		allEdges.append(C_FeFlowEdg(tetrahedronlist[t * 4 + 0], tetrahedronlist[t * 4 + 3], this->allEdges.length()));
		allEdges.append(C_FeFlowEdg(tetrahedronlist[t * 4 + 1], tetrahedronlist[t * 4 + 2], this->allEdges.length()));
		allEdges.append(C_FeFlowEdg(tetrahedronlist[t * 4 + 1], tetrahedronlist[t * 4 + 3], this->allEdges.length()));
		allEdges.append(C_FeFlowEdg(tetrahedronlist[t * 4 + 2], tetrahedronlist[t * 4 + 3], this->allEdges.length()));
	}
	// Sort the allEdges first by n1 than by n2 than by index. The index is the previous position in the list, which gives the general order of edges.
	qSort(allEdges.begin(), allEdges.end(), sortEdgesByNodes);
	this->allEdgesWithoutDuplicates.clear();
	// Transfer all individual edges to allEdgesWithoutDuplicates. If edges are duplicated, only the one with the lowest index will be transfered
	for (int e = 0; e != this->allEdges.length(); e++)
		if (allEdgesWithoutDuplicates.length() == 0 ||
			allEdges[e].nodes[0] != allEdgesWithoutDuplicates.last().nodes[0] ||
			allEdges[e].nodes[1] != allEdgesWithoutDuplicates.last().nodes[1])
			allEdgesWithoutDuplicates.append(allEdges[e]);
	// Reindex the right ordered list based on the order. The new index is now the feflow index
	for (int e = 0; e != this->allEdgesWithoutDuplicates.length(); e++)
		this->allEdgesWithoutDuplicates[e].index = e;
}

void
C_FeFlow::generateUndefinedEdges(int marker, int numberofedges, QList <int> edgemarkerlist, QList <long> edgelist)
{
	this->undefinedEdges.clear();
	// Store all edges of tetgen having m as marker m in undefinedEdges list. The index of all edges is set to ZERO.
	for (int e = 0; e != numberofedges; e++)
		if (edgemarkerlist[e] == marker)
			undefinedEdges.append(C_FeFlowEdg(edgelist[e * 2 + 0], edgelist[e * 2 + 1],0));
	// Sort undefinedEdges by nodes and index. This support the search of each individual edge in the allEdgesWihoutDuplicates list.
	qSort(undefinedEdges.begin(), undefinedEdges.end(), sortEdgesByNodes);
}

void
C_FeFlow::generateDefinedEdges()
{
	// Search for each undefined edge in allEdgesWithoutDuplicates (this list has the required feflow index). Found -> the edge with the feflow index is copied to the definedEdges list
	this->definedEdges.clear();
	int lastPos = 0;
	for (int u = 0; u != this->undefinedEdges.length(); u++)
		for (int a = lastPos; a != this->allEdgesWithoutDuplicates.length(); a++)
			if (this->undefinedEdges[u].nodes[0] == this->allEdgesWithoutDuplicates[a].nodes[0] &&
				this->undefinedEdges[u].nodes[1] == this->allEdgesWithoutDuplicates[a].nodes[1])
			{
				this->definedEdges.append(this->allEdgesWithoutDuplicates[a]);
				lastPos = a;
				break;
			}
	// The list with the defined edges is sorted by the index (to have a not required but sorted output)
	qSort(definedEdges.begin(), definedEdges.end(), sortEdgesByIndex);
}
