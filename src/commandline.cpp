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

#include "commandline.h"
#include "geometry.h"

C_Model CmdModel;

C_CommandLine::C_CommandLine(QCommandLineParser * parser)
{
	if (parser->isSet("input"))
	{
		CmdModel.FileNameModel = parser->value("input");
		if (QFileInfo(CmdModel.FileNameModel).suffix() != "pvd")
			CmdModel.FileNameModel = QFileInfo(CmdModel.FileNameModel).path() + "/" + QFileInfo(CmdModel.FileNameModel).baseName() + ".pvd";
		CmdModel.Open();
	}
	if (parser->isSet("p"))
	{
		this->preMeshJob();
		for (int s = 0; s != CmdModel.Surfaces.length(); s++)
			for (int c = 0; c != CmdModel.Surfaces[s].Constraints.length(); c++)
				CmdModel.Surfaces[s].Constraints[c].Type = "SEGMENTS";
		for (int p = 0; p != CmdModel.Polylines.length(); p++)
			for (int c = 0; c != CmdModel.Polylines[p].Constraints.length(); c++)
				CmdModel.Polylines[p].Constraints[c].Type = "SEGMENTS";
	}
	if (parser->isSet("m"))
		this->MeshJob();
	if (parser->isSet("output"))
	{
		CmdModel.FileNameModel = parser->value("output");
		if (QFileInfo(CmdModel.FileNameModel).suffix() != "pvd")
			CmdModel.FileNameModel = QFileInfo(CmdModel.FileNameModel).path() + "/" + QFileInfo(CmdModel.FileNameModel).baseName() + ".pvd";
		CmdModel.Save();
	}
	if (parser->isSet("export-vtu"))
	{
		CmdModel.FileNameTmp = parser->value("exportvtu");
		if (QFileInfo(CmdModel.FileNameTmp).suffix() != "vtu")
			CmdModel.FileNameTmp = QFileInfo(CmdModel.FileNameTmp).path() + "/" + QFileInfo(CmdModel.FileNameTmp).baseName() + ".vtu";
		CmdModel.ExportVTU3D();
	}
}

C_CommandLine::~C_CommandLine()
{}

void
C_CommandLine::preMeshJob()
{
	//QDateTime startdate, enddate;
	//startdate = QDateTime::currentDateTime();
	//std::cout << ">Start Time: " << startdate.toString().toUtf8().constData() << std::endl;;
	// convex hull
	int currentStep = 0, totalSteps = CmdModel.Surfaces.length();
	for (int s=0;s!=CmdModel.Surfaces.length();s++)
	{
		C_CmdTask *task = new C_CmdTask(this, "CONVEXHULL", s, 0, ++currentStep, totalSteps);
		QThreadPool::globalInstance()->start(task);
	}
	QThreadPool::globalInstance()->waitForDone();
	// segments
	currentStep = 0;
	totalSteps = CmdModel.Polylines.length();
	for (int p=0;p!=CmdModel.Polylines.length();p++)
	{
		C_CmdTask *task = new C_CmdTask(this, "SEGMENTS", p, 0, ++currentStep, totalSteps);
		QThreadPool::globalInstance()->start(task);
	}
	QThreadPool::globalInstance()->waitForDone();
	//triangulation (coarse)
	currentStep = 0;
	totalSteps = CmdModel.Surfaces.length();
	for (int s=0;s!=CmdModel.Surfaces.length();s++)
	{
		C_CmdTask *task = new C_CmdTask(this, "TRIANGLES", s, 0, ++currentStep, totalSteps);
		QThreadPool::globalInstance()->start(task);
	}
	QThreadPool::globalInstance()->waitForDone();
	// intersection: surface-surface
	CmdModel.Intersections.clear();
	currentStep = 0;
	totalSteps = CmdModel.Surfaces.length()*(CmdModel.Surfaces.length() - 1) / 2;
	if (totalSteps>0)
	{
		for (int s1=0;s1!=CmdModel.Surfaces.length()-1;s1++)
			for (int s2=s1+1;s2!=CmdModel.Surfaces.length();s2++)
			{
				C_CmdTask *task = new C_CmdTask(this, "INTERSECTION_MESH_MESH", s1, s2, ++currentStep, totalSteps);
				QThreadPool::globalInstance()->start(task);
			}
		QThreadPool::globalInstance()->waitForDone();
	}
	// intersection: surface-polyline
	currentStep = 0;
	totalSteps = CmdModel.Polylines.length()*CmdModel.Surfaces.length();
	if (totalSteps>0)
	{
		for (int p=0;p!=CmdModel.Polylines.length();p++)
			for (int s=0;s!=CmdModel.Surfaces.length();s++)
			{
				C_CmdTask *task = new C_CmdTask(this, "INTERSECTION_POLYLINE_MESH", p, s, ++currentStep, totalSteps);
				QThreadPool::globalInstance()->start(task);
			}
		QThreadPool::globalInstance()->waitForDone();
	}

	CmdModel.calculate_size_of_intersections();

	//intersection: triple points
	CmdModel.TPs.clear();
	currentStep = 0;
	totalSteps = CmdModel.Intersections.length()*(CmdModel.Intersections.length() - 1) / 2;
	if (totalSteps>0)
	{
		for (int i1=0;i1!=CmdModel.Intersections.length()-1;i1++)
			for (int i2 = i1 + 1; i2 != CmdModel.Intersections.length(); i2++)
			{
				C_CmdTask * task = new C_CmdTask(this, "INTERSECTION_TRIPLEPOINTS", i1, i2, ++currentStep, totalSteps);
				QThreadPool::globalInstance()->start(task);
			}
		QThreadPool::globalInstance()->waitForDone();
	}
	CmdModel.insert_int_triplepoints();
	// aligning convex hull to intersection
	for (int s=0;s!=CmdModel.Surfaces.length();s++)
		CmdModel.Surfaces[s].alignIntersectionsToConvexHull();
	// constraints
	for (int s=0;s!=CmdModel.Surfaces.length();s++)
		CmdModel.Surfaces[s].calculate_Constraints();
	for (int p=0;p!=CmdModel.Polylines.length();p++)
		CmdModel.Polylines[p].calculate_Constraints();
	CmdModel.calculate_size_of_constraints();
	//enddate = QDateTime::currentDateTime();
}

void
C_CommandLine::MeshJob()
{
//	QDateTime startdate, enddate;
//	startdate = QDateTime::currentDateTime();
	//emit progress_append(">Start Time: " + startdate.toString() + "\n");
	int currentStep, totalSteps;
	// segments - fine
	currentStep = 0;
	totalSteps = CmdModel.Polylines.length();
	for (int p = 0; p != CmdModel.Polylines.length(); p++)
	{
		C_CmdTask *task = new C_CmdTask(this, "SEGMENTS_FINE", p, 0, ++currentStep, totalSteps);
		QThreadPool::globalInstance()->start(task);
	}
	QThreadPool::globalInstance()->waitForDone();
	// triangulation - fine
	currentStep = 0;
	totalSteps = CmdModel.Surfaces.length();
	for (int s = 0; s != CmdModel.Surfaces.length(); s++)
	{
		C_CmdTask *task = new C_CmdTask(this, "TRIANGLES_FINE", s, 0, ++currentStep, totalSteps);
		QThreadPool::globalInstance()->start(task);
	}
	QThreadPool::globalInstance()->waitForDone();
	// tetrahedralization
	CmdModel.calculate_tets("pq1.2AY");
	//	enddate = QDateTime::currentDateTime();
}

void
C_CommandLine::runThreadPool(QString Attribute, int Object1, int Object2, int currentStep, int totalSteps)
{
	if (Attribute == "CONVEXHULL")
	{
		CmdModel.Surfaces[Object1].calculate_normal_vector();
		CmdModel.Surfaces[Object1].rotate(true);
		CmdModel.Surfaces[Object1].calculate_convex_hull();
		CmdModel.Surfaces[Object1].interpolation("ConvexHull", CmdModel.intAlgorythm);
		CmdModel.Surfaces[Object1].rotate(false);

	}
	if (Attribute == "SEGMENTS")
	{
		CmdModel.Polylines[Object1].calculate_segments(false);
		CmdModel.Polylines[Object1].Intersections.clear();
		CmdModel.Polylines[Object1].Path.calculate_min_max();
	}
	if (Attribute == "SEGMENTS_FINE")
	{
		CmdModel.Polylines[Object1].calculate_segments(true);
		CmdModel.Polylines[Object1].Path.calculate_min_max();
	}
	if (Attribute == "TRIANGLES")
	{
		CmdModel.Surfaces[Object1].rotate(true);
		CmdModel.Surfaces[Object1].calculate_triangles(false);
		CmdModel.Surfaces[Object1].interpolation("Mesh", CmdModel.intAlgorythm);
		CmdModel.Surfaces[Object1].rotate(false);
		CmdModel.Surfaces[Object1].Intersections.clear();
		CmdModel.Surfaces[Object1].calculate_min_max();
		for (int t = 0; t != CmdModel.Surfaces[Object1].Ts.length(); t++)
		{
			CmdModel.Surfaces[Object1].Ts[t].calculate_min_max();
			CmdModel.Surfaces[Object1].Ts[t].setNormalVector();
		}
	}
	if (Attribute == "TRIANGLES_FINE")
	{
		CmdModel.Surfaces[Object1].calculate_normal_vector();
		CmdModel.Surfaces[Object1].rotate(true);
		CmdModel.Surfaces[Object1].separate_Constraints();
		CmdModel.Surfaces[Object1].calculate_triangles(true);
		CmdModel.Surfaces[Object1].interpolation("Mesh", CmdModel.intAlgorythm);
		CmdModel.Surfaces[Object1].rotate(false);
		for (int t = 0; t != CmdModel.Surfaces[Object1].Ts.length(); t++)
		{
			CmdModel.Surfaces[Object1].Ts[t].calculate_min_max();
			CmdModel.Surfaces[Object1].Ts[t].setNormalVector();
		}
	}
	if (Attribute == "INTERSECTION_POLYLINE_MESH")
		CmdModel.calculate_int_point(Object1, Object2);
	if (Attribute == "INTERSECTION_MESH_MESH")
		CmdModel.calculate_int_polyline(Object1, Object2);
	if (Attribute == "INTERSECTION_TRIPLEPOINTS")
		CmdModel.calculate_int_triplepoints(Object1, Object2);
}
