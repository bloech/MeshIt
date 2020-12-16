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

#include <sstream>

#include "core.h"
#include "geometry.h"
#include "intersections.h"

#define SQUAREROOTTWO 1.4142135623730950488016887242096980785696718753769480732
#define MY_PI 3.141592653589793238462643383279502884197169399375105820974944592308

/* This mutex is required to protected concurrent accesses on shared list
 * objects. Use a global mutex for now to avoid more extensive changes to the
 * code base. */
QMutex mutex;

/********** Class C_Eigenvalue **********/

void
C_Eigenvalue::ComputeEigenvalue()
{
	Tridiagonal();

	if( QLAlgorithm() )
	{
		DecreasingSort();
		GuaranteeRotation();
	}
	else
		exit(EXIT_FAILURE);
}

void
C_Eigenvalue::Tridiagonal()
{
	double SM_00 = Element[0][0];
	double SM_01 = Element[0][1];
	double SM_02 = Element[0][2];
	double SM_11 = Element[1][1];
	double SM_12 = Element[1][2];
	double SM_22 = Element[2][2];

	Diag[0] = SM_00;
	Subd[2] = 0;
	
	if( SM_02!=0.0 )
	{
		double Length = std::sqrt(SM_01*SM_01 + SM_02*SM_02);
		double InvLength = 1.0 / Length;

		SM_01 *= InvLength;
		SM_02 *= InvLength;
		
		double Q = 2.0*SM_01*SM_12 + SM_02*(SM_22 - SM_11);
		
		Diag[1] = SM_11 + SM_02*Q;
		Diag[2] = SM_22 - SM_02*Q;
		
		Subd[0] = Length;
		Subd[1] = SM_12 - SM_01*Q;
		
		Element[0][0] = 1.0;
		Element[0][1] = 0.0;
		Element[0][2] = 0.0;
		Element[1][0] = 0.0;
		Element[1][1] = SM_01;
		Element[1][2] = SM_02;
		Element[2][0] = 0.0;
		Element[2][1] = SM_02;
		Element[2][2] = -SM_01;
		
		IsRotation = false;
	}
	else
	{
		Diag[1] = SM_11;
		Diag[2] = SM_22;
		
		Subd[0] = SM_01;
		Subd[1] = SM_12;
	  
		Element[0][0] = 1.0;
		Element[0][1] = 0.0;
		Element[0][2] = 0.0;
		Element[1][0] = 0.0;
		Element[1][1] = 1.0;
		Element[1][2] = 0.0;
		Element[2][0] = 0.0;
		Element[2][1] = 0.0;
		Element[2][2] = 1.0;

		IsRotation = true;
	}
}

bool
C_Eigenvalue::QLAlgorithm()
{
	const int iMaxIter = 32;

	for(int i0=0; i0!=3; i0++)
	{
		int i1;
		for(i1=0;i1!=iMaxIter;i1++)
		{
			int i2;
			for(i2=i0;i2<=(3-2);i2++)
			{
				double Tmp = FABS(Diag[i2]) + FABS(Diag[i2 + 1]);
				if( FABS(Subd[i2]) + Tmp == Tmp )
					break;
			}
			
			if( i2 == i0 )
				break;
			
			double G = (Diag[i0 + 1] - Diag[i0])/(2.0*Subd[i0]);
			double R = std::sqrt(G*G + 1.0);

			if( G < 0.0 )
				G = Diag[i2] - Diag[i0] + Subd[i0]/(G - R);
			else
				G = Diag[i2] - Diag[i0] + Subd[i0]/(G + R);
			
			double Sin = 1.0, Cos = 1.0, P = 0.0;
			
			for(int i3=i2-1; i3>=i0; i3--)
			{
				double F = Sin*Subd[i3];
				double B = Cos*Subd[i3];

				if( FABS(F)>=FABS(G) )
				{
					Cos = G/F;
					R = std::sqrt(Cos*Cos + 1.0);
					Subd[i3 + 1] = F*R;
					Sin = 1.0/R;
					Cos *= Sin;
				}
				else
				{
					Sin = F/G;
					R = std::sqrt(Sin*Sin + 1.0);
					Subd[i3 + 1] = G*R;
					Cos = 1.0/R;
					Sin *= Cos;
				}
				
				G = Diag[i3 + 1] - P;
				R = (Diag[i3] - G)*Sin + 2.0*B*Cos;
				P = Sin*R;
				Diag[i3 + 1] = G + P;
				G = Cos*R - B;
				
				for(int i4=0;i4!=3;i4++)
				{
					F = Element[i4][i3 + 1];
					Element[i4][i3 + 1] = Sin*Element[i4][i3] + Cos*F;
					Element[i4][i3] = Cos*Element[i4][i3] - Sin*F;
				}
			}
			Diag[i0] -= P;
			Subd[i0] = G;
			Subd[i2] = 0.0;
		}
		if( i1 == iMaxIter )
			return false;
	}
	return true;
}

void
C_Eigenvalue::DecreasingSort()
{
	/*Sort eigenvalues in decreasing order*/
	for(int i0=0,i1;i0<=3-2;i0++)
	{
		/*locate maximum eigenvalue*/
		i1 = i0;
		double Max = Diag[i1];
		int i2;
		for(i2=i0+1;i2!=3;i2++)
		{
			if( Diag[i2]>Max )
			{
				i1 = i2;
				Max = Diag[i1];
			}
		}
		if( i1!=i0 )
		{
			/*swap eigenvalues*/
			Diag[i1] = Diag[i0];
			Diag[i0] = Max;
			/*swap eigenvectors*/
			for(i2=0;i2!=3;i2++)
			{
				double Tmp = Element[i2][i0];
				Element[i2][i0] = Element[i2][i1];
				Element[i2][i1] = Tmp;
				IsRotation = !IsRotation;
			}
		}
	}
}

void
C_Eigenvalue::GuaranteeRotation()
{
	if( !IsRotation )
	{
		/*change sign on the first column*/
		for(int iRow=0;iRow!=3;iRow++)
			Element[iRow][0] = -Element[iRow][0];
	}
}

/********** Class C_Colors **********/

C_Colors::C_Colors()
{
	/*RED*/
	this->Red[0] = 1.0f; this->Red[1] = 0.0f; this->Red[2] = 0.0f; this->Red[3] = 1.0f;
	this->LightRed[0]=1.0f; this->LightRed[1]=0.75f; this->LightRed[2]=0.5f; this->LightRed[3]=1.0f; 
	this->RedTrans[0] = 1.0f; this->RedTrans[1] = 0.0f; this->RedTrans[2] = 0.0f; this->RedTrans[3] = 0.8f;
	/*GREEN*/
	this->Green[0] = 0.0f; this->Green[1] = 1.0f; this->Green[2] = 0.0f; this->Green[3] = 1.0f;
	this->LightGreen[0]=0.75f; this->LightGreen[1]=1.0f; this->LightGreen[2]=0.5f; this->LightGreen[3]=1.0f;
	this->GreenTrans[0] = 0.0f; this->GreenTrans[1] = 1.0f; this->GreenTrans[2] = 0.0f; this->GreenTrans[3] = 0.8f;
	/*BLUE*/
	this->Blue[0] = 0.0f; this->Blue[1] = 0.0f; this->Blue[2] = 1.0f; this->Blue[3] = 1.0f;
	this->LightBlue[0]=0.5f; this->LightBlue[1]=0.75f; this->LightBlue[2]=1.0f; this->LightBlue[3]=1.0f;
	this->BlueTrans[0] = 0.0f; this->BlueTrans[1] = 0.0f; this->BlueTrans[2] = 1.0f; this->BlueTrans[3] = 0.8f;
	/*MAGENTA*/
	this->Magenta[0] = 1.0f; this->Magenta[1] = 0.0f; this->Magenta[2] = 1.0f; this->Magenta[3] = 1.0f;
	this->LightMagenta[0]=1.0f; this->LightMagenta[1]=0.5f; this->LightMagenta[2]=1.0f; this->LightMagenta[3]=1.0f;
	this->MagentaTrans[0] = 1.0f; this->MagentaTrans[1] = 0.0f; this->MagentaTrans[2] = 1.0f; this->MagentaTrans[3] = 0.8f;
	/*YELLOW*/
	this->Yellow[0]=1.0f;this->Yellow[1]=1.0f;this->Yellow[2]=0.0f;this->Yellow[3]=1.0f;
	this->YellowTrans[0] = 1.0f; this->YellowTrans[1] = 1.0f; this->YellowTrans[2] = 0.0f; this->YellowTrans[3] = 0.8f;
	/*CYAN*/
	this->Cyan[0]=0.0f;this->Cyan[1]=1.0f;this->Cyan[2]=1.0f;this->Cyan[3]=1.0f;
	this->CyanTrans[0] = 0.0f; this->CyanTrans[1] = 1.0f; this->CyanTrans[2] = 1.0f; this->CyanTrans[3] = 0.8f;
	/*WHITE*/
	this->White[0]=1.0f;this->White[1]=1.0f;this->White[2]=1.0f;this->White[3]=1.0f;
	this->WhiteTrans[0] = 1.0f; this->WhiteTrans[1] = 1.0f; this->WhiteTrans[2] = 1.0f; this->WhiteTrans[3] = 0.1f;
	/*GREY*/
	this->Grey[0]=0.2f;this->Grey[1]=0.2f;this->Grey[2]=0.2f;this->Grey[3]=1.0f;
	/*ORANGE*/
	this->Orange[0] = 1.0f; this->Orange[1] = 0.5f; this->Orange[2] = 0.0f; this->Orange[3] = 1.0f;

	/* Mint Tulip */
	this->MintTulip[0] = 191.f / 255.f;
	this->MintTulip[1] = 247.f / 255.f;
	this->MintTulip[2] = 226.f / 255.f;
	this->MintTulip[3] = 1.f;
}

/********** Class C_Material **********/

C_Material::C_Material()
{
	this->drawMatFaces = false;
	this->drawMatEdges = false;
}


/********** Class C_VTU **********/

void
C_VTU::clear()
{
	this->NumberOfPoints = 0;
	this->NumberOfCells = 0;

	this->Points.clear();
	this->connectivity.clear();
	this->offsets.clear();
	this->types.clear();
	this->matType.clear();
	this->matR.clear();
	this->matG.clear();
	this->matB.clear();
	this->Object1.clear();
	this->Object2.clear();
	this->pointData.clear();
}

void
C_VTU::write(QString FileName)
{
	QFile file(FileName);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
		return;
	QTextStream out(&file);
	out.setRealNumberPrecision(24);
	/*Write the header*/
	out << "<?xml version=\"1.0\"?>\n";
	out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out << "  <UnstructuredGrid>\n";
	out << "    <Piece NumberOfPoints=\"" << this->NumberOfPoints << "\" NumberOfCells=\"" << this->NumberOfCells << "\">\n";
	/*Points block*/
	out << "      <Points>\n";
	out << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";
	for (int p = 0; p != this->Points.length(); p++)
		if (p % 10 == 0)
			out << endl << "          " << this->Points[p].x() << " " << this->Points[p].y() << " " << this->Points[p].z();
		else
			out << " " << this->Points[p].x() << " " << this->Points[p].y() << " " << this->Points[p].z();
	out << endl;
	out << "        </DataArray>\n";
	out << "      </Points>\n";
	/*Cells block*/
	out << "      <Cells>\n";
	out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
	out << "          ";
	for (int c = 0; c != this->connectivity.length(); c++)
		out << this->connectivity[c] << " ";
	out << "\n";
	out << "        </DataArray>\n";
	out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	out << "          ";
	for (int o = 0; o != this->offsets.length(); o++)
		out << this->offsets[o] << " ";
	out << "\n";
	out << "        </DataArray>\n";
	out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	out << "          ";
	for (int t = 0; t != this->types.length(); t++)
		out << this->types[t] << " ";
	out << "\n";
	out << "        </DataArray>\n";
	out << "      </Cells>\n";
	/*PointData block*/
	out << "      <PointData>\n";
	/*PointType is Set for ConvexHull or Intersections*/
	if (this->pointData.length() != 0)
	{
		out << "        <DataArray type=\"Int32\" Name=\"pointData\" format=\"ascii\">\n";
		out << "          ";
		for (int p = 0; p != this->pointData.length(); p++)
			out << this->pointData[p] << " ";
		out << "\n";
		out << "        </DataArray>\n";
	}
	out << "      </PointData>\n";
	/*CellData block*/
	out << "      <CellData>\n";
	/*Here the intersections belonging to a mesh are stored*/
	if (this->Object1.length() != 0 || this->Object2.length() != 0)
	{
		out << "        <DataArray type=\"Int32\" Name=\"Object1\" format=\"ascii\">" << "\n";
		out << "          ";
		for (int m = 0; m != this->Object1.length(); m++)
			out << this->Object1[m] << " ";
		out << "\n";
		out << "        </DataArray>" << "\n";
		out << "        <DataArray type=\"Int32\" Name=\"Object2\" format=\"ascii\">" << "\n";
		out << "          ";
		for (int m = 0; m != this->Object2.length(); m++)
			out << this->Object2[m] << " ";
		out << "\n";
		out << "        </DataArray>" << "\n";
	}
	/*matType is defined for constraints and for Tets*/
	if (this->matType.length() != 0)
	{
		out << "        <DataArray type=\"Int32\" Name=\"matType\" format=\"ascii\">" << "\n";
		out << "          ";
		for (int m = 0; m != this->matType.length(); m++)
			out << this->matType[m] << " ";
		out << "\n";
		out << "        </DataArray>" << "\n";
	}
	/*MatRGB and size are defined for constraints*/
	if (this->matR.length() != 0 || this->matG.length() != 0 || this->matB.length() != 0)
	{
		out << "        <DataArray type=\"Int32\" Name=\"red\" format=\"ascii\">" << "\n";
		out << "          ";
		for (int m = 0; m != this->matR.length(); m++)
			out << this->matR[m] << " ";
		out << "\n";
		out << "        </DataArray>" << "\n";
		out << "        <DataArray type=\"Int32\" Name=\"green\" format=\"ascii\">" << "\n";
		out << "          ";
		for (int m = 0; m != this->matG.length(); m++)
			out << this->matG[m] << " ";
		out << "\n";
		out << "        </DataArray>" << "\n";
		out << "        <DataArray type=\"Int32\" Name=\"blue\" format=\"ascii\">" << "\n";
		out << "          ";
		for (int m = 0; m != this->matB.length(); m++)
			out << this->matB[m] << " ";
		out << "\n";
		out << "        </DataArray>" << "\n";
	}
	out << "      </CellData>\n";
	out << "    </Piece>\n";
	out << "  </UnstructuredGrid>\n";
	out << "</VTKFile>\n";

	file.close();
	this->clear();
}

void
C_VTU::read(QString FileName)
{
	QString line;
	int connectivities;
	double d_value;
	int i_value;

	this->clear();

	QFile file(FileName);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
	QTextStream in(&file);

	while (!in.atEnd())
	{
		line = in.readLine().simplified();
		if (line.contains("offsets"))
		{
			line = in.readLine().simplified();
			connectivities = line.section(" ", -1, -1).toInt();
			break;
		}
	}

	in.seek(0);
	
	while (!in.atEnd())
	{
		line = in.readLine().simplified();
		if (line.contains("<Piece"))
		{
			this->NumberOfPoints = line.section("\"", 1, 1).toInt();
			this->NumberOfCells = line.section("\"", 3, 3).toInt();
		}
		if (line.contains("NumberOfComponents"))
		{
			for (int p = 0; p != this->NumberOfPoints; p++)
			{
				C_Vector3D Po;
				in >> d_value;
				Po.setX(d_value);
				in >> d_value;
				Po.setY(d_value);
				in >> d_value;
				Po.setZ(d_value);
				this->Points.append(Po);
			}
		}
		if (line.contains("connectivity"))
		{
			for (int c = 0; c != connectivities; c++)
			{
				in >> i_value;
				this->connectivity.append(i_value);
			}
		}
		if (line.contains("offsets"))
		{
			for (int c = 0; c != this->NumberOfCells; c++)
			{
				in >> i_value;
				this->offsets.append(i_value);
			}
		}
		if (line.contains("types"))
		{
			for (int c = 0; c != this->NumberOfCells; c++)
			{
				in >> i_value;
				this->types.append(i_value);
			}
		}
		if (line.contains("pointData"))
		{
			for (int p = 0; p != this->NumberOfPoints; p++)
			{
				in >> i_value;
				this->pointData.append(i_value);
			}
		}
		if (line.contains("Object1"))
		{
			for (int c = 0; c != this->NumberOfCells; c++)
			{
				in >> i_value;
				this->Object1.append(i_value);
			}
		}
		if (line.contains("Object2"))
		{
			for (int c = 0; c != this->NumberOfCells; c++)
			{
				in >> i_value;
				this->Object2.append(i_value);
			}
		}
		if (line.contains("matType"))
		{
			for (int c = 0; c != this->NumberOfCells; c++)
			{
				in >> i_value;
				this->matType.append(i_value);
			}
		}
		if (line.contains("\"red\""))
		{
			for (int c = 0; c != this->NumberOfCells; c++)
			{
				in >> i_value;
				this->matR.append(i_value);
			}
		}
		if (line.contains("\"green\""))
		{
			for (int c = 0; c != this->NumberOfCells; c++)
			{
				in >> i_value;
				this->matG.append(i_value);
			}
		}
		if (line.contains("\"blue\""))
		{
			for (int c = 0; c != this->NumberOfCells; c++)
			{
				in >> i_value;
				this->matB.append(i_value);
			}
		}
	}
}

/********** Class C_Line **********/

void
C_Line::AddPosition()
{
/*
	Add a length position to the points of a convex hull.
	The first point of the convex hull has a position of zero. 
	The next point of the convex hull has a position of the previous point plus the length of the segment between these two points.
	Therefore the position increases from one point to the other by the length of the segment and the position is stored in NsPos.
*/
	if (this->Ns.isEmpty()) 
		return;
	
	NsPos.clear();
	double Position = 0;
	NsPos.append(Position);
	
	for (int n = 0; n != Ns.length() - 1; n++)
	{
		double part = (Ns[n + 1].x() - Ns[n].x())*(Ns[n + 1].x() - Ns[n].x()) +
					  (Ns[n + 1].y() - Ns[n].y())*(Ns[n + 1].y() - Ns[n].y()) +
					  (Ns[n + 1].z() - Ns[n].z())*(Ns[n + 1].z() - Ns[n].z());
		Position += std::sqrt(part);
		NsPos.append(Position);
	}
}

void
C_Line::AddPoint(C_Vector3D TP)
{
/*
	Add a point to the belonging line
*/
	C_Vector3D *xProjection = new C_Vector3D;
	for (int n = 0; n != this->Ns.length() - 1; n++)
	{
		projectTo(TP, this->Ns[n], this->Ns[n + 1], xProjection);
		if (lengthSquared(*xProjection) != 0 && lengthSquared(*xProjection - TP) < 1e-24)
		{
			this->Ns.insert(n + 1, TP);
			return;
		}
	}
}

bool
C_Line::isInside(C_Vector3D N)
{
/*
	Point-in-polygon algorithm, created especially for World-Wide Web servers to process image maps with mouse-clickable regions.
	The original source code [INPOLY.C] was developed by Bob Stein and Craig Yap and can be found at (http://www.visibone.com/inpoly/inpoly.c)
*/
	C_Vector3D oldN, newN;
	C_Vector3D N1, N2;
	bool inside = false;

	if (Ns.length() < 3)
		return false;

	oldN = Ns[0];

	for (int n = 1; n != Ns.length(); n++)
	{
		newN = Ns[n];
		if (newN.x() > oldN.x())
		{
			N1 = oldN;
			N2 = newN;
		}
		else
		{
			N1 = newN;
			N2 = oldN;
		}

		if (
			(newN.x() < N.x()) == (N.x() <= oldN.x())	/* edge "open" at left end */
			&&
			((long)N.y() - (long)N1.y())*(long)(N2.x() - N1.x()) < ((long)N2.y() - (long)N1.y())*(long)(N.x() - N1.x())
			)
			inside = !inside;

		oldN = newN;
	}
	return(inside);
}

void
C_Line::calculate_min_max()
{
/*
	Calculate the minimum and maximum extensions to set the Line.
*/
	if (this->Ns.isEmpty())
		return;
	
	this->min = this->Ns.first();
	this->max = this->Ns.first();
	
	for (int n = 0; n != this->Ns.length(); n++)
	{
		if (this->Ns[n].x()<this->min.x())
			this->min.setX(this->Ns[n].x());
		if (this->Ns[n].y()<this->min.y())
			this->min.setY(this->Ns[n].y());
		if (this->Ns[n].z()<this->min.z())
			this->min.setZ(this->Ns[n].z());

		if (this->Ns[n].x()>this->max.x())
			this->max.setX(this->Ns[n].x());
		if (this->Ns[n].y()>this->max.y())
			this->max.setY(this->Ns[n].y());
		if (this->Ns[n].z()>this->max.z())
			this->max.setZ(this->Ns[n].z());
	}
}

void
C_Line::Invert()
{
/*
	Invert the place and position of the points within the list of a Convex Hull and calculate the new positions
*/
	for (int n = 0; n != Ns.length(); n++)
		Ns.move(n, 0);
	this->AddPosition();
}

void
C_Line::EraseSpecialPoints()
{
/*
	Set Type for all points to the default value [DEFAULT].
*/
	for (int n = 0; n != Ns.length(); n++)
		Ns[n].setType("DEFAULT");
}

void
C_Line::MakeCornersSpecial()
{
/*
	Sets the type for a corner point of a Convex Hull.
	It first checks if there is a hard corner (less than 135 degree).
*	If found, it updates the Type to CORNER.
*/
	C_Vector3D diff1, diff2, point;
	
	for (int n = 0; n != Ns.length() - 1; n++)
	{
		if (n == 0)
		{
			diff1 = Ns[Ns.length() - 2] - Ns[0];
			diff2 = Ns[1] - Ns[0];
		}
		else
		{
			diff1 = Ns[n - 1] - Ns[n];
			diff2 = Ns[n + 1] - Ns[n];
		}
		
		normalize(&diff1);
		normalize(&diff2);
		double alpha = dot(diff1, diff2);
		
		if (alpha > (-0.5*SQUAREROOTTWO))
		{
			Ns[n].setType("CORNER");
			if (n == 0)
				Ns.last().setType("CORNER");
		}
	}
}

bool
C_Line::SortByType(QString FLAG)
{
/*
	Sort points of a Convex Hull based on their types [FLAG].
	As long as the first point of the C_Line has a different type than FLAG, it will be pushed to the end of the list.
*/
	this->Ns.removeFirst(); /*the first point is identically with the last and therfore removed*/
	
	for (int n = 0; n != this->Ns.length(); n++)
	{
		if (this->Ns.first().type() != FLAG)
			this->Ns.move(0, this->Ns.length() - 1); /*point is not a special type --> push it to the end*/
		else
		{
			this->Ns.append(this->Ns.first()); /*the first point is duplicated at the end*/
			return true;
		}
	}
	this->Ns.append(this->Ns.first()); /*the first point is duplicated at the end*/
	return false;
}

void
C_Line::CleanIdenticalPoints()
{
/*
	Delete identical points from the line.
*/
	if (this->Ns.length()<2)
		return;

	for (int n = 0; n != Ns.length() - 1; n++)
	{
		if (lengthSquared(this->Ns[n] - this->Ns[n + 1])<1e-24)
		{
			if (Ns[n].type() == "DEFAULT")
				this->Ns.removeAt(n);
			else
				this->Ns.removeAt(n + 1);
			n--;
		}
	}
}

void
C_Line::RefineByLength(double length)
{
/*
	Refine an existing line by length.
	It first deletes all points which has a Type equal to DEFAULT.
	It then subdivided the new parts in intervals of length equals to length inserting new points in the list.
*/
	if (this->Ns.length()<2)
		return;
	
	C_Vector3D newPoint;
	C_Line newLine;
	newLine.Ns.append(Ns.first());
	newLine.NsPos.append(NsPos.first());
	
	for (int n = 1; n != Ns.length() - 1; n++)
	{
		if (Ns[n].type() != "DEFAULT")
		{
			newLine.Ns.append(Ns[n]);
			newLine.NsPos.append(NsPos[n]);
		}
	}
	newLine.Ns.append(Ns.last());
	newLine.NsPos.append(NsPos.last());

	for (int n = 0; n != newLine.Ns.length() - 1; n++)
	{
		double tms1 = (newLine.NsPos[n + 1] - newLine.NsPos[n]) / length;
		double tms2 = std::floor((newLine.NsPos[n + 1] - newLine.NsPos[n]) / length);
		int points;
		if ((tms1 - tms2)<.5)
			points = std::floor((newLine.NsPos[n + 1] - newLine.NsPos[n]) / length);
		else
			points = std::ceil((newLine.NsPos[n + 1] - newLine.NsPos[n]) / length);
		double minLength = (newLine.NsPos[n + 1] - newLine.NsPos[n]) / points;
		double tmpPos = newLine.NsPos[n];
		if (points>1)
		{
			for (int c = 1; c != points; c++)
			{
				tmpPos += minLength;
				newPoint = this->getPointAtPos(tmpPos);
				newPoint.setType("DEFAULT");
				newLine.Ns.insert(n + c, newPoint);
				newLine.NsPos.insert(n + c, tmpPos);
			}
			n += points - 1;
		}
	}
	Ns.clear();

	for (int n = 0; n != newLine.Ns.length(); n++)
		Ns.append(newLine.Ns[n]);

	this->AddPosition();
}

bool
C_Line::IsIdenticallyWith(C_Line Line)
{
/*
	Determine whether an existing line is identical to the given Line
*/
	bool found;
	if (Ns.length() != Line.Ns.length())
		return false;
	if (lengthSquared(Ns.first() - Line.Ns.first())>1e-24)
		Line.Invert();
	for (int n = 0; n != Line.Ns.length(); n++)
	{
		found = false;
		for (int n2 = 0; n2 != Ns.length(); n2++)
			if (lengthSquared(Ns[n2] - Line.Ns[n])<1e-24)
				found = true;
		if (!found)
			return false;
	}
	return true;
}

void
C_Line::appendNonExistingSegment(C_Vector3D N1, C_Vector3D N2)
{
/*
	Append a segment [N1, N2] to an existing C_Line (only if not already listed).
*/
	bool append = true;

	if (lengthSquared(N1 - N2) < 1e-24)	/*points are identical --> not a segment*/
		append = false;

	for (int n = 0; n<Ns.length(); n = n + 2)
	{
		if (lengthSquared(N1 - Ns[n])<1e-24 && lengthSquared(N2 - Ns[n + 1])<1e-24)
			append = false;
		if (lengthSquared(N2 - Ns[n])<1e-24 && lengthSquared(N1 - Ns[n + 1])<1e-24)
			append = false;
	}

	if (append)
	{
		Ns.append(N1);
		Ns.append(N2);
	}
}

void
C_Line::GenerateFirstSplineOfSegments(C_Line *input)
{
/*
	Generate an ordered polyline (no duplicated points of common intersections) defined by the intersection points between two surfaces.
	Update the intersection given the IDs of the two intersecting surfaces.
*/
	this->appendNonExistingSegment(input->Ns[0], input->Ns[1]);
	input->Ns.removeAt(1);
	input->Ns.removeAt(0);

	/*look for the first node of the new list in the old one*/
	bool reloop = true;
	while (reloop)
	{
		reloop = false;
		for (int n = 0; n != input->Ns.length(); n++)
		{
			if (lengthSquared(input->Ns[n] - Ns.first())<1e-24)
			{
				if (n % 2 == 0)
				{
					Ns.prepend(input->Ns[n + 1]);
					input->Ns.removeAt(n + 1);
					input->Ns.removeAt(n);
					reloop = true;
					break;
				}
				else
				{
					Ns.prepend(input->Ns[n - 1]);
					input->Ns.removeAt(n);
					input->Ns.removeAt(n - 1);
					reloop = true;
					break;
				}
			}
		}
	}
	reloop = true;
	while (reloop)
	{
		reloop = false;
		for (int n = 0; n != input->Ns.length(); n++)
		{
			if (lengthSquared(input->Ns[n] - Ns.last())<1e-24)
			{
				if (n % 2 == 0)
				{
					Ns.append(input->Ns[n + 1]);
					input->Ns.removeAt(n + 1);
					input->Ns.removeAt(n);
					reloop = true;
					break;
				}
				else
				{
					Ns.append(input->Ns[n - 1]);
					input->Ns.removeAt(n);
					input->Ns.removeAt(n - 1);
					reloop = true;
					break;
				}
			}
		}
	}

/*	QList<C_Vector3D>::iterator i;
	double x, y, z;
	for (i = Ns.begin(); i != Ns.end(); ++i)
	{
		x = i->x();
		y = i->y();
		z = i->z();
	}
*/
	this->AddPosition();
}

C_Line
C_Line::calculateSkewLineTransversal(const C_Vector3D &P1, const C_Vector3D &P2, const C_Vector3D &Q1, const C_Vector3D &Q2)
{
/*
	Calculates the shortest connector [L0] between two skew lines [L1 : P1+s*(P2-P1)] and [L2 : Q1+t*(Q2-Q1)]
	The shortest connector is defined as the transversal (i.e. common intersection line) of minimum length (i.e. perpendicular to each lines).
	Step 1. Calculate the normal [R21= P21xQ21] vector common to both segments.
	If the two lines are parallel (common normal is equal to zero) or it is not possible to connect --> return empty C_Line.
	Step 2. Calculate the plane equation defined by one of the two lines (i.e. C_Line L1) and the common normal vector [E1 : P1+s*P21+v*R21]
	Step 3. Calculate the intersection between the plane [E1] and the other C_Line [L2].
	Solve the corresponding system of three linear equations for s, t and v.
	Once having determined s and t, the two points defining the Shortest Connector are calculated by inserting those parameters into the equation of the two C_Line L1 and L2.
	return:
	If the two lines are intersecting -> C_Line containing twice the intersection point.
	If the two lines are skew		 --> shortest connector.
*/

	C_Line shortestConnector;

	C_Vector3D P21 = C_Vector3D(P2 - P1);
	C_Vector3D Q21 = C_Vector3D(Q2 - Q1);
	C_Vector3D R21; /* The common normal*/
	cross(P21, Q21, &R21);
	C_Vector3D PQ1 = C_Vector3D(P1 - Q1);

	if (length(R21) == 0)
	{
		/*The 2 lines are parallel --> return empty C_Line*/
		shortestConnector.Ns.clear();
		return shortestConnector;
	}

	double D = P21.x()*Q21.y()*R21.z() + Q21.x()*R21.y()*P21.z() + R21.x()*P21.y()*Q21.z() - R21.x()*Q21.y()*P21.z() - Q21.x()*P21.y()*R21.z() - P21.x()*R21.y()*Q21.z();
	double Ds = R21.x()*Q21.y()*PQ1.z() + Q21.x()*PQ1.y()*R21.z() + PQ1.x()*R21.y()*Q21.z() - PQ1.x()*Q21.y()*R21.z() - Q21.x()*R21.y()*PQ1.z() - R21.x()*PQ1.y()*Q21.z();
	double Dt = P21.x()*PQ1.y()*R21.z() + PQ1.x()*R21.y()*P21.z() + R21.x()*P21.y()*PQ1.z() - R21.x()*PQ1.y()*P21.z() - PQ1.x()*P21.y()*R21.z() - P21.x()*R21.y()*PQ1.z();

	double s = Ds / D;
	double t = Dt / D;

	if (0 <= s && s<1)
		if (0 <= t && t<1)
		{
			shortestConnector.Ns.append(P1 + s*P21);
			shortestConnector.Ns.append(Q1 + t*Q21);
			return shortestConnector;
		}

	return shortestConnector;
}

C_Vector3D
C_Line::getPointAtPos(double Pos)
{
/*
	Return the point at the given postion [Pos] in the list of points describing the C_Line.
	It first determines the two consecutive points within the list which contain the new point.
	Once found, it determines the new coordinates of the point based on the distances of the new point with respect to the existing ones.
*/
	double ratio;
	C_Vector3D newPoint;

	for (int n=0;n!=Ns.length()-1;n++)
		if (NsPos[n]<=Pos && Pos<NsPos[n+1])
		{
			ratio = (Pos - NsPos[n]) / (NsPos[n + 1] - NsPos[n]);
			newPoint = Ns[n + 1] - Ns[n];
			newPoint *= ratio;
			newPoint += Ns[n];
		}

	return newPoint;
}

/********** Class C_Polyline **********/

C_Polyline::C_Polyline()
{
	this->drawScatteredData = false;
	this->drawEdges = false;
	this->drawVertices = false;
	this->drawIntVertices = false;
	this->drawConstraints = false;
}

void
C_Polyline::calculate_segments(bool withConstraints)
{
	this->Path.Ns = this->SDs;
	
	Path.CleanIdenticalPoints();
	Path.AddPosition();
	Path.RefineByLength(this->size);
	
	if (withConstraints)
	{
		for (int c = 0; c != this->Constraints.length(); c++)
		{
			if (this->Constraints[c].Type != "UNDEFINED")
			{
				if (lengthSquared(this->Constraints[c].Ns.first() - Path.Ns.first())<1e-24)
					Path.Ns.first().setType("INTERSECTION_POINT");
				else if (lengthSquared(this->Constraints[c].Ns.first() - Path.Ns.last())<1e-24)
					Path.Ns.last().setType("INTERSECTION_POINT");
				else
					Path.AddPoint(Constraints[c].Ns.first());
			}
		}
		
		while (Path.Ns.length() != 0 && Path.Ns.first().type() != "INTERSECTION_POINT")
			Path.Ns.removeFirst();
		while (Path.Ns.length() != 0 && Path.Ns.last().type() != "INTERSECTION_POINT")
			Path.Ns.removeLast();
		
		Path.CleanIdenticalPoints();
		Path.AddPosition();
		Path.RefineByLength(this->size);
	}
}

void
C_Polyline::calculate_Constraints()
{
	this->Constraints.clear();

	C_Line Constraint;
	int lastPos;
	lastPos = 0;
	Constraint.RGB[0] = 0;
	Constraint.RGB[1] = 0;
	Constraint.RGB[2] = 0;

	Constraint.Type = "UNDEFINED";

	Constraint.Ns.clear();
	Constraint.Ns.append(this->Path.Ns.first());
	if (++Constraint.RGB[0]>255)
	{
		Constraint.RGB[0] = 0; 
		if (++Constraint.RGB[1]>255)
		{ 
			Constraint.RGB[1] = 0; 
			if (++Constraint.RGB[2]>255)
			{}
		}
	}
	this->Constraints.append(Constraint);

	Constraint.Ns.clear();
	Constraint.Ns.append(this->Path.Ns.last());
	if (++Constraint.RGB[0]>255)
	{
		Constraint.RGB[0] = 0;
		if (++Constraint.RGB[1]>255)
		{
			Constraint.RGB[1] = 0;
			if (++Constraint.RGB[2]>255)
			{}
		}
	}
	this->Constraints.append(Constraint);

	for (int i = 0; i != this->Intersections.length(); i++)
	{
		Constraint.Ns = Intersections[i]->Ns;
		if (++Constraint.RGB[0]>255)
		{
			Constraint.RGB[0] = 0;
			if (++Constraint.RGB[1]>255)
			{
				Constraint.RGB[1] = 0;
				if (++Constraint.RGB[2]>255)
				{}
			}
		}
		this->Constraints.append(Constraint);
	}
}

void
C_Polyline::calculate_min_max()
{
	this->min = this->SDs.first();
	this->max = this->SDs.first();

	for (int s = 0; s != this->SDs.length(); s++)
	{
		if (this->SDs[s].x()<this->min.x())
			this->min.setX(this->SDs[s].x());
		if (this->SDs[s].y()<this->min.y())
			this->min.setY(this->SDs[s].y());
		if (this->SDs[s].z()<this->min.z())
			this->min.setZ(this->SDs[s].z());
		if (this->SDs[s].x()>this->max.x())
			this->max.setX(this->SDs[s].x());
		if (this->SDs[s].y()>this->max.y())
			this->max.setY(this->SDs[s].y());
		if (this->SDs[s].z()>this->max.z())
			this->max.setZ(this->SDs[s].z());
	}
	if (this->size == 0)
		this->size = length(max - min) / 16;
}

void
C_Polyline::makeScatteredData()
{
	GLuint list = glGenLists(1);
	glNewList(list, GL_COMPILE);
	glDisable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.White);

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, this->Cols.LightMagenta);

	glBegin(GL_POINTS);
	for (int s = 0; s != SDs.length(); s++)
		glVertex3f(SDs[s].x(), SDs[s].y(), SDs[s].z());
	glEnd();

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.Grey);
	glEnable(GL_LIGHT0);
	glEndList();
	listScatteredData = list;
}

void
C_Polyline::makeEdges()
{
	GLuint list = glGenLists(1);
	glNewList(list, GL_COMPILE);

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, Cols.LightMagenta);

	glBegin(GL_LINE_STRIP);
	for (int n = 0; n != this->Path.Ns.length(); n++)
		glVertex3d(this->Path.Ns[n].x(), this->Path.Ns[n].y(), this->Path.Ns[n].z());
	glEnd();

	glEndList();
	listEdges = list;
}

void
C_Polyline::makeVertices()
{
	GLuint list = glGenLists(1);
	glNewList(list, GL_COMPILE);
	glDisable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.White);

	glBegin(GL_POINTS);
	for (int n = 0; n != this->Path.Ns.length(); n++)
	{
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, this->Cols.LightMagenta);
		glVertex3d(this->Path.Ns[n].x(), this->Path.Ns[n].y(), this->Path.Ns[n].z());
	}
	glEnd();

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.Grey);
	glEnable(GL_LIGHT0);
	glEndList();
	listVertices = list;
}

void
C_Polyline::makeIntVertices()
{
	GLuint list = glGenLists(1);
	glNewList(list, GL_COMPILE);
	glDisable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.White);

	glBegin(GL_POINTS);
	for (int i = 0; i != Intersections.length(); i++)
	{
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, this->Cols.White);
		glVertex3d(Intersections[i]->Ns.first().x(), Intersections[i]->Ns.first().y(), Intersections[i]->Ns.first().z());
	}
	glEnd();

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.Grey);
	glEnable(GL_LIGHT0);
	glEndList();
	listIntVertices = list;
}

void
C_Polyline::makeConstraints()
{
	GLuint list = glGenLists(1);
	glNewList(list, GL_COMPILE);
	glDisable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.White);

	for (int s = 0; s != Constraints.length(); s++)
	{
		if (Constraints[s].Type == "SEGMENTS")
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, this->Cols.Blue);
		if (Constraints[s].Type == "UNDEFINED")
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, this->Cols.White);

		glColor3f(Constraints[s].RGB[0] / 255.0f, Constraints[s].RGB[1] / 255.0f, Constraints[s].RGB[2] / 255.0f);

		glBegin(GL_POINTS);
		glVertex3f(Constraints[s].Ns.first().x(), Constraints[s].Ns.first().y(), Constraints[s].Ns.first().z());
		glEnd();
		glColor3f(0.0f, 0.0f, 0.0f);
	}

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.Grey);
	glEnable(GL_LIGHT0);
	glEndList();
	listConstraints = list;
}

void
C_Polyline::makeVTU_SD()
{
	this->VTU.clear();
	this->VTU.NumberOfPoints = this->SDs.length();
	this->VTU.NumberOfCells = 1;
	this->VTU.Points = this->SDs;
	for (int s = 0; s != this->SDs.length(); s++) 
		this->VTU.connectivity.append(s);
	this->VTU.offsets.append(this->SDs.length());
	this->VTU.types.append(2);
}

void
C_Polyline::makeVTU_SEG()
{
	this->VTU.clear();
	this->VTU.NumberOfPoints = this->Path.Ns.length();
	this->VTU.NumberOfCells = 1;
	this->VTU.Points.append(this->Path.Ns);
	for (int n = 0; n != this->Path.Ns.length(); n++)
		this->VTU.connectivity.append(n);
	this->VTU.offsets.append(this->Path.Ns.length());
	this->VTU.types.append(4);
}

void
C_Polyline::makeVTU_CON()
{
	this->VTU.clear();
	this->VTU.NumberOfPoints = this->Constraints.length();
	this->VTU.NumberOfCells = this->Constraints.length();
	for (int c = 0; c != this->Constraints.length(); c++)
	{
		this->VTU.Points.append(this->Constraints[c].Ns);
		this->VTU.connectivity.append(c);
		this->VTU.matR.append(this->Constraints[c].RGB[0]);
		this->VTU.matG.append(this->Constraints[c].RGB[1]);
		this->VTU.matB.append(this->Constraints[c].RGB[2]);
		if (this->Constraints[c].Type == "SEGMENTS")
			this->VTU.matType.append(0);
		if (this->Constraints[c].Type == "UNDEFINED")
			this->VTU.matType.append(1);
		this->VTU.offsets.append(c + 1);
		this->VTU.types.append(1);
	}
}






















/********** Class C_Triangle **********/

void
C_Triangle::calculate_min_max()
{
/*
	Calculate the minimum and maximum extensions to set the Triangle.
*/
	this->min = *this->Ns[0];
	this->max = *this->Ns[0];

	if (this->Ns[1]->x()<this->min.x())
		this->min.setX(this->Ns[1]->x());
	if (this->Ns[2]->x()<this->min.x())
		this->min.setX(this->Ns[2]->x());
	if (this->Ns[1]->y()<this->min.y())
		this->min.setY(this->Ns[1]->y());
	if (this->Ns[2]->y()<this->min.y())
		this->min.setY(this->Ns[2]->y());
	if (this->Ns[1]->z()<this->min.z())
		this->min.setZ(this->Ns[1]->z());
	if (this->Ns[2]->z()<this->min.z())
		this->min.setZ(this->Ns[2]->z());
	if (this->Ns[1]->x()>this->max.x())
		this->max.setX(this->Ns[1]->x());
	if (this->Ns[2]->x()>this->max.x())
		this->max.setX(this->Ns[2]->x());
	if (this->Ns[1]->y()>this->max.y())
		this->max.setY(this->Ns[1]->y());
	if (this->Ns[2]->y()>this->max.y())
		this->max.setY(this->Ns[2]->y());
	if (this->Ns[1]->z()>this->max.z())
		this->max.setZ(this->Ns[1]->z());
	if (this->Ns[2]->z()>this->max.z())
		this->max.setZ(this->Ns[2]->z());
}

void
C_Triangle::setNormalVector()
{
/*
	Calculate the normal vector to the plane of the Triangle given its three vertices.
*/
	normal(*this->Ns[0], *this->Ns[1], *this->Ns[2], &this->normal_vector);
}

















void C_Surface::Crout_LU_Decomposition_with_Pivoting(double *A, int * pivot){
	int i, j, k;
	double *p_k, *p_row, *p_col;
	double max;
//	For each row and column (k = 0, ..., n-1) ...//
	for(k=0,p_k=A;k!=(this->selectedSDs.length()+1);p_k+=(this->selectedSDs.length()+1),k++){
//	... find the pivot row ... //
		pivot[k]=k;
		max=FABS(*(p_k+k));
		for(j=k+1,p_row=p_k+(this->selectedSDs.length()+1);j!=(this->selectedSDs.length()+1);j++,p_row+=(this->selectedSDs.length()+1)){
			if(max<FABS(*(p_row+k))){
				max=FABS(*(p_row+k));
				pivot[k]=j;
				p_col=p_row;
			}
		}
//	... and if the pivot row differs from the current row, then interchange the two rows ... //
		if(pivot[k]!=k){
			for(j=0;j!=(this->selectedSDs.length()+1);j++){
				max=*(p_k+j);
				*(p_k+j)=*(p_col+j);
				*(p_col+j)=max;
			}
		}
//	... and if the matrix is singular, return error !! ... //
		if(*(p_k+k)==0.0){
			break;
		}
//	... otherwise find the upper triangular matrix elements for row k ... //
		for(j=k+1;j!=(this->selectedSDs.length()+1);j++){
			*(p_k+j)/=*(p_k+k);
		}
//	... and update remaining matrix. //
		for(i=k+1,p_row=p_k+(this->selectedSDs.length()+1);i!=(this->selectedSDs.length()+1);p_row+=(this->selectedSDs.length()+1),i++)
			for(j=k+1;j!=(this->selectedSDs.length()+1);j++)
				*(p_row+j)-=*(p_row+k) * *(p_k+j);
	}
}

void C_Surface::Crout_LU_with_Pivoting_Solve(double * LU, double * B, double * x, int * pivot){
	int i, k ;
	double *p_k;
	double dum;
//	... Solve the linear equation Lx = B for x, where L is a lower triangular matrix. //
	for(k=0,p_k=LU;k!=(this->selectedSDs.length()+1);p_k+=(this->selectedSDs.length()+1),k++){
		if(pivot[k]!=k){
			dum=B[k];
            B[k]=B[pivot[k]];
            B[pivot[k]]=dum;
        }
		x[k]=B[k];
        for(i=0;i<k;i++) x[k]-=x[i] * *(p_k + i);
        x[k]/=*(p_k + k);
	}
//	... Solve the linear equation Ux = y, where y is the solution obtained above of Lx = B and U is an upper triangular matrix. 
//		The diagonal part of the upper triangular part of the matrix is assumed to be 1.0. //
	for (k=(this->selectedSDs.length()+1)-1,p_k=LU+(this->selectedSDs.length()+1)*((this->selectedSDs.length()+1)-1);k>=0;k--,p_k-=(this->selectedSDs.length()+1)){
		if(pivot[k]!=k){
			dum=B[k];
            B[k]=B[pivot[k]];
            B[pivot[k]]=dum;
        }
		for(i=k+1;i!=(this->selectedSDs.length()+1);i++)  x[k]-=x[i]* *(p_k+i);
        if(*(p_k+k)==0.0) break;
	}
}

void C_Surface::fillKriging(){
//	Select the points for interpolation //
	int step=this->SDs.length()/256+1;
	this->selectedSDs.clear();
	for (int s=0;s<this->SDs.length();s+=step){
		this->selectedSDs.append(&this->SDs[s]);
	}
// Allocate memory for the arrays//
	this->MatrixKriging=new double[(this->selectedSDs.length()+1)*(this->selectedSDs.length()+1)];
	this->VectorKriging=new double[(this->selectedSDs.length()+1)];
	this->KrigingWeights=new double[(this->selectedSDs.length()+1)];
	this->KrigingPivot=new int[(this->selectedSDs.length()+1)];
//	Copy data //
	this->KrigingBeta=1.5;
	for(int s=0;s!=this->selectedSDs.length();s++){
		this->VectorKriging[s]=this->selectedSDs[s]->z();
	}
	this->VectorKriging[(this->selectedSDs.length())]=0.0;
//	Fill the Kriging_r2VLU matrix //
	for(int s1=0;s1!=this->selectedSDs.length();s1++){
		for(int s2=s1+1;s2!=this->selectedSDs.length();s2++){
			double r2=(this->selectedSDs[s1]->x()-this->selectedSDs[s2]->x())*(this->selectedSDs[s1]->x()-this->selectedSDs[s2]->x())+
						(this->selectedSDs[s1]->y()-this->selectedSDs[s2]->y())*(this->selectedSDs[s1]->y()-this->selectedSDs[s2]->y());
			this->MatrixKriging[s1*(this->selectedSDs.length()+1)+s2]=this->MatrixKriging[s2*(this->selectedSDs.length()+1)+s1]=r2;
		}
	}
	for(int s1=0;s1!=this->selectedSDs.length();s1++){
		this->MatrixKriging[s1*(this->selectedSDs.length()+2)]=0.0;
		for(int s2=s1+1;s2!=this->selectedSDs.length();s2++){
			this->MatrixKriging[s1*(this->selectedSDs.length()+1)+s2]=this->MatrixKriging[s2*(this->selectedSDs.length()+1)+s1]=
			pow(this->MatrixKriging[s1*(this->selectedSDs.length()+1)+s2],0.5*this->KrigingBeta);
		}
		this->MatrixKriging[s1*(this->selectedSDs.length()+1)+this->selectedSDs.length()]=this->MatrixKriging[this->selectedSDs.length()*(this->selectedSDs.length()+1)+s1] = 1.0;
	}
	this->MatrixKriging[this->selectedSDs.length()*(this->selectedSDs.length()+2)]=0.0;
//	LU decomposition of the Kriging_r2VLU matrix //
	this->Crout_LU_Decomposition_with_Pivoting(this->MatrixKriging, this->KrigingPivot);
//	Solve the Matrix System //
	this->Crout_LU_with_Pivoting_Solve(this->MatrixKriging, this->VectorKriging, this->KrigingWeights, this->KrigingPivot);
}

void C_Surface::fillSpline(){
	double r2;
	int step=this->SDs.length()/256+1;
	this->selectedSDs.clear();
	this->MatrixSpline.clear();
	this->VectorSpline.clear();
	for (int s=0;s<this->SDs.length();s+=step){
		this->selectedSDs.append(&this->SDs[s]);  
		this->VectorSpline.append(this->SDs[s].z());
	}
	for(int n=0;n!=3;n++) this->VectorSpline.append(0);


// Fill Matrix coefficient matrix//
	QList<double> row;
	for(int r=0;r!=this->selectedSDs.length();r++){
		row.clear(); 
		for (int c=0;c!=this->selectedSDs.length();c++){
			if (c==r){
				row.append(0.0);
			}else{
				r2=(this->selectedSDs[r]->x() - this->selectedSDs[c]->x())*(this->selectedSDs[r]->x() - this->selectedSDs[c]->x()) + (this->selectedSDs[r]->y() - this->selectedSDs[c]->y())*(this->selectedSDs[r]->y() - this->selectedSDs[c]->y());
					row.append(0.5*r2*log(r2));
			}
		}
		row.append(1.0);
		row.append(this->selectedSDs[r]->x());
		row.append(this->selectedSDs[r]->y());
		MatrixSpline.append(row); 
	}
	row.clear(); 
	for (int c=0;c!=this->selectedSDs.length();c++) row.append(1.0); 
	for(int n=0;n!=3;n++) row.append(0); 
	MatrixSpline.append(row); 
	row.clear(); 
	for (int c=0;c!=this->selectedSDs.length();c++) row.append(this->selectedSDs[c]->x()); 
	for(int n=0;n!=3;n++) row.append(0); 
	MatrixSpline.append(row); 
	row.clear(); 
	for (int c=0;c!=this->selectedSDs.length();c++) row.append(this->selectedSDs[c]->y()); 
	for(int n=0;n!=3;n++) row.append(0); 
	MatrixSpline.append(row);
	this->SolveSpline();
}

void C_Surface::SolveSpline(){
	double big,dum,pivinv;
	int icol,irow,l;
	QList<int> indxc,indxr,ipiv;

	indxc.clear();
	indxr.clear();
	ipiv.clear();
	for(int i=0;i!=this->VectorSpline.length();i++){
		indxc.append(0);
		indxr.append(0);
		ipiv.append(0);
	}
	icol=irow=0;
	for(int i=0;i!=this->VectorSpline.length();i++){
		big = 0.0;
		for(int j=0;j!=this->VectorSpline.length();j++){
			if((ipiv[j]!=1)){
				for(int k=0;k!=this->VectorSpline.length();k++){
					if((ipiv[k]==0)){
						if((FABS(this->MatrixSpline[j][k])>=big)){
							big=FABS(this->MatrixSpline[j][k]);
							irow=j;
							icol=k;
						}
					}else if(ipiv[k]>1){
						printf("pause 1 in GAUSSJ - singular matrix\n");
						break;
					}
				}
			}
		}
		ipiv[icol]=ipiv[icol]+1;
		if((irow!=icol)){
			for(int l=0;l!=this->VectorSpline.length();l++){
				dum=this->MatrixSpline[irow][l];
				this->MatrixSpline[irow][l]=this->MatrixSpline[icol][l];
				this->MatrixSpline[icol][l]=dum;
			}
			dum=this->VectorSpline[irow];
			this->VectorSpline[irow]=this->VectorSpline[icol];
			this->VectorSpline[icol] = dum;
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if(this->MatrixSpline[icol][icol] == 0.0){
			printf("pause 2 in GAUSSJ - singular matrix\n");
			break;
		}
		pivinv=1.0/(this->MatrixSpline[icol][icol]);
		this->MatrixSpline[icol][icol]=1.0;
		for(l=0;l!=this->VectorSpline.length();l++){
			this->MatrixSpline[icol][l]=this->MatrixSpline[icol][l]*pivinv;
		}
		this->VectorSpline[icol]=this->VectorSpline[icol]*pivinv;
		for(int ll=0;ll!=this->VectorSpline.length();ll++){
			if((ll!=icol)){
				dum=this->MatrixSpline[ll][icol];
				this->MatrixSpline[ll][icol]=0.0;
				for(int l=0;l!=this->VectorSpline.length();l++){
					this->MatrixSpline[ll][l]=this->MatrixSpline[ll][l]-this->MatrixSpline[icol][l]*dum;
				}
				this->VectorSpline[ll]=this->VectorSpline[ll]-this->VectorSpline[icol]*dum;
			}
		}
	}
	for(int j=0;j!=this->VectorSpline.length();j++){
		l=this->VectorSpline.length()-1-j;
		if(indxr[l]!=indxc[l]){
			for(int k=0;k!=this->VectorSpline.length();k++){
				dum=this->MatrixSpline[k][indxr[l]];
				this->MatrixSpline[k][indxr[l]]=this->MatrixSpline[k][indxc[l]];
				this->MatrixSpline[k][indxc[l]]=dum;
			}
		}
	}
}

double C_Surface::KRIGING(double x, double y){
	double z=0.0;
	for(int s=0;s!=this->selectedSDs.length();s++){
		z+=this->KrigingWeights[s]*pow((x-this->selectedSDs[s]->x())*(x-this->selectedSDs[s]->x())+(y-this->selectedSDs[s]->y())*(y-this->selectedSDs[s]->y()),0.5*this->KrigingBeta);
	}
	z+=this->KrigingWeights[this->selectedSDs.length()];
	return z;
}

double C_Surface::SPLINE(double x, double y){
	double distance;
	for (int s=0;s!=this->SDs.length(); s++){
		if (this->SDs[s].x()==x && this->SDs[s].y()==y) return this->SDs[s].z();
	}
	double z=this->VectorSpline[this->selectedSDs.length()]+this->VectorSpline[this->selectedSDs.length()+1]*x+this->VectorSpline[this->selectedSDs.length()+2]*y;
	for (int s=0;s!=this->selectedSDs.length(); s++){
		distance=(this->selectedSDs[s]->x()-x)*(this->selectedSDs[s]->x()-x)+(this->selectedSDs[s]->y()-y)*(this->selectedSDs[s]->y()-y);
		z+=0.5*this->VectorSpline[s]*distance*log(distance);
	}
	return z;
}

double C_Surface::IDW(double x, double y){
	double sum=0;
	double z=0;
	double distance;
	QList<double> lambda;

	for (int s=0;s!=this->SDs.length(); s++){
		if (this->SDs[s].x()==x && this->SDs[s].y()==y) return this->SDs[s].z();
		distance=(this->SDs[s].x()-x)*(this->SDs[s].x()-x)+(this->SDs[s].y()-y)*(this->SDs[s].y()-y);
		lambda.append(1/(distance*distance));
	}
	for (int l=0; l!=lambda.length();l++){
		sum+=lambda[l];
		z+=lambda[l]*this->SDs[l].z();
	}
	lambda.clear();
	return z/sum;
}

void C_Surface::interpolation(QString object, QString method){
	if (object=="Mesh"){
		if(method=="IDW"){
			for (int n=0; n!=this->Ns.length(); n++){
					Ns[n].setZ(IDW(this->Ns[n].x(),this->Ns[n].y()));
				}
		}else if(method=="SPLINE"){
			this->fillSpline();
			for (int n=0; n!=this->Ns.length(); n++){
				Ns[n].setZ(SPLINE(this->Ns[n].x(),this->Ns[n].y()));
			}
		}else if(method=="KRIGING"){
			this->fillKriging();
			for (int n=0; n!=this->Ns.length(); n++){
				Ns[n].setZ(KRIGING(this->Ns[n].x(),this->Ns[n].y()));
			}
		}
	}
	if (object=="ConvexHull"){ 
		if(method=="IDW"){	
			for (int n=0; n!=this->ConvexHull.Ns.length(); n++){
				this->ConvexHull.Ns[n].setZ(IDW(this->ConvexHull.Ns[n].x(),this->ConvexHull.Ns[n].y()));
			}
		}else if(method=="SPLINE"){
			this->fillSpline();
			for (int n=0; n!=this->ConvexHull.Ns.length(); n++){
				this->ConvexHull.Ns[n].setZ(SPLINE(this->ConvexHull.Ns[n].x(),this->ConvexHull.Ns[n].y()));
			}
		}else if(method=="KRIGING"){
			this->fillKriging();
			for (int n=0; n!=this->ConvexHull.Ns.length(); n++){
				this->ConvexHull.Ns[n].setZ(KRIGING(this->ConvexHull.Ns[n].x(),this->ConvexHull.Ns[n].y()));
			}
		}
	}
} 

void C_Surface::rotate(bool onto_z){
	C_Vector3D * axis = new C_Vector3D;
	if (onto_z){
		cross(this->normal_vector, C_Vector3D(0, 0, 1), axis);
	}
	else{
		cross(C_Vector3D(0, 0, 1), this->normal_vector, axis);
		for (int n = 0; n != Ns.length(); n++){
			for (int nc = 0; nc != ConvexHull.Ns.length(); nc++){
				if (ConvexHull.Ns[nc].x() == Ns[n].x() && ConvexHull.Ns[nc].y() == Ns[n].y()){
					Ns[n].setZ(ConvexHull.Ns[nc].z());
				}
			}
			for (int i = 0; i != Intersections.length(); i++){
				for (int ni = 0; ni != Intersections[i]->Ns.length(); ni++){
					if (Intersections[i]->Ns[ni].x() == Ns[n].x() && Intersections[i]->Ns[ni].y() == Ns[n].y()){
						Ns[n].setZ(Intersections[i]->Ns[ni].z());
					}
				}
			}
			for (int s = 0; s != Constraints.length(); s++){
				for (int ns = 0; ns != Constraints[s].Ns.length(); ns++){
					if (Constraints[s].Ns[ns].x() == Ns[n].x() && Constraints[s].Ns[ns].y() == Ns[n].y()){
						Ns[n].setZ(Constraints[s].Ns[ns].z());
					}
				}
			}
		}
	}
	normalize(axis);
	double costheta = dot(C_Vector3D(0, 0, 1), this->normal_vector);
	if (1 - FABS(costheta) < 1e-12) return;
	double c = costheta;
	double s = sqrt(1 - c*c);
	double C = 1 - c;
	C_Vector3D r1 = C_Vector3D(axis->x()*axis->x()*C + c, axis->x()*axis->y()*C - axis->z()*s, axis->x()*axis->z()*C + axis->y()*s);
	C_Vector3D r2 = C_Vector3D(axis->y()*axis->x()*C + axis->z()*s, axis->y()*axis->y()*C + c, axis->y()*axis->z()*C - axis->x()*s);
	C_Vector3D r3 = C_Vector3D(axis->z()*axis->x()*C - axis->y()*s, axis->z()*axis->y()*C + axis->x()*s, axis->z()*axis->z()*C + c);
	C_Vector3D temp;

	for (int n = 0; n != Ns.length(); n++){
		temp = Ns[n];
		Ns[n].setX(dot(r1, temp));
		Ns[n].setY(dot(r2, temp));
		Ns[n].setZ(dot(r3, temp));
	}

	for (int s = 0; s != SDs.length(); s++){
		temp = SDs[s];
		SDs[s].setX(dot(r1, temp));
		SDs[s].setY(dot(r2, temp));
		SDs[s].setZ(dot(r3, temp));
	}
	for (int n = 0; n != ConvexHull.Ns.length(); n++){
		temp = ConvexHull.Ns[n];
		ConvexHull.Ns[n].setX(dot(r1, temp));
		ConvexHull.Ns[n].setY(dot(r2, temp));
		ConvexHull.Ns[n].setZ(dot(r3, temp));
	}
	for (int i = 0; i != Intersections.length(); i++){
		for (int n = 0; n != Intersections[i]->Ns.length(); n++){
			temp = Intersections[i]->Ns[n];
			Intersections[i]->Ns[n].setX(dot(r1, temp));
			Intersections[i]->Ns[n].setY(dot(r2, temp));
			Intersections[i]->Ns[n].setZ(dot(r3, temp));
		}
	}
	for (int s = 0; s != Constraints.length(); s++){
		for (int n = 0; n != Constraints[s].Ns.length(); n++){
			temp = Constraints[s].Ns[n];
			Constraints[s].Ns[n].setX(dot(r1, temp));
			Constraints[s].Ns[n].setY(dot(r2, temp));
			Constraints[s].Ns[n].setZ(dot(r3, temp));
		}
	}
}

void C_Surface::calculate_min_max(){
	this->min=this->SDs.first();  
	this->max=this->SDs.first();  
	for (int s=0;s!=this->SDs.length();s++){
		if (this->SDs[s].x()<this->min.x()) this->min.setX(this->SDs[s].x());
		if (this->SDs[s].y()<this->min.y()) this->min.setY(this->SDs[s].y());
		if (this->SDs[s].z()<this->min.z()) this->min.setZ(this->SDs[s].z());
		if (this->SDs[s].x()>this->max.x()) this->max.setX(this->SDs[s].x());
		if (this->SDs[s].y()>this->max.y()) this->max.setY(this->SDs[s].y());
		if (this->SDs[s].z()>this->max.z()) this->max.setZ(this->SDs[s].z());
	}
	if (this->size == 0) size = length(max - min) / 16;
}

/*! \brief Calculating best fitting plane based on a set of C_Vector3D Points.
	\details 
		Note that this is the "ordinary least squares" fit, which is appropriate only when z is expected to be a linear function of x and y. 
		If you are looking more generally for a "best fit plane" in 3-space, you may want to learn about "geometric" least squares.

		Note also that this will fail if your points are in a plane parallel to z-direction.
		For this case the z-component of the normal vecto is 0 and the x,y- components are calculated by use of the x and y values of the scattered data points.

*/

void C_Surface::calculate_normal_vector(){
//	Determine centre of the plane //
	double centre[3] = {0, 0, 0};
	double weight = 0.0;
	for(int n=0;n!=SDs.length();n++){
		centre[0] += SDs[n].x();
		centre[1] += SDs[n].y();
		centre[2] += SDs[n].z();
		weight ++;
	}
	double recip_w = 1.0/weight;
	centre[0] *= recip_w;
	centre[1] *= recip_w;
	centre[2] *= recip_w;
//	Summing the squares of the differences
	double sum_xx = 0.0, sum_xy = 0.0, sum_xz = 0.0;
	double sum_yy = 0.0, sum_yz = 0.0;
	double sum_zz = 0.0;
	double diff[3];
	for(int n=0;n!=SDs.length();n++){
		diff[0] = (SDs[n].x() - centre[0]);
		diff[1] = (SDs[n].y() - centre[1]);
		diff[2] = (SDs[n].z() - centre[2]);
		sum_xx += diff[0]*diff[0];
		sum_xy += diff[0]*diff[1];
		sum_xz += diff[0]*diff[2];
		sum_yy += diff[1]*diff[1];
		sum_yz += diff[1]*diff[2];
		sum_zz += diff[2]*diff[2];
	}
	sum_xx *= recip_w;
	sum_xy *= recip_w;
	sum_xz *= recip_w;
	sum_yy *= recip_w;
	sum_yz *= recip_w;
	sum_zz *= recip_w;
//	Set up the Eigensolver
	C_Eigenvalue k_es;
	k_es.Element[0][0] = sum_xx;
	k_es.Element[0][1] = sum_xy;
	k_es.Element[0][2] = sum_xz;
	k_es.Element[1][0] = sum_xy;
	k_es.Element[1][1] = sum_yy;
	k_es.Element[1][2] = sum_yz;
	k_es.Element[2][0] = sum_xz;
	k_es.Element[2][1] = sum_yz;
	k_es.Element[2][2] = sum_zz;
//	Compute eigenvalue stuff, smallest eigenvalue is in the last position ...
	k_es.ComputeEigenvalue();
	double k_normal[3];
	double plane[4];
	k_normal[0] = k_es.Element[0][2];
	k_normal[1] = k_es.Element[1][2];
	k_normal[2] = k_es.Element[2][2];
// Here use the minimum energy ... this part for MeshIT is deprecated ... this gives the plane ...
	plane[0] = k_normal[0]; // a
	plane[1] = k_normal[1]; // b 
	plane[2] = k_normal[2]; // c
	double dot = k_normal[0]*centre[0] + k_normal[1]*centre[1] + k_normal[2]*centre[2];
	plane[3] = 0 - dot; // d
//	here we fill the proper normal vector
	normal_vector.setX(k_normal[0]);
	normal_vector.setY(k_normal[1]);
	normal_vector.setZ(k_normal[2]);
}

void C_Surface::calculate_convex_hull(){
	ConvexHull.Ns.clear(); 
	QList<C_Vector3D> mSDs(SDs);
	//int mem_n;
	int mem_m=0;
	for (int m=0;m!=mSDs.length();m++){
		if (mSDs[m].y()<mSDs[mem_m].y() || (mSDs[m].y()==mSDs[mem_m].y() && mSDs[m].x()<mSDs[mem_m].x())) mem_m=m;
	}
	ConvexHull.Ns.append(mSDs[mem_m]);
	mSDs.move(mem_m,mSDs.length()-1);

	while (lengthSquared(ConvexHull.Ns.first()-ConvexHull.Ns.last())!=0 || ConvexHull.Ns.length()<3){
		mem_m=0;
		for (int m=1;m!=mSDs.length();m++){
			double ccw=(mSDs[mem_m].x() - ConvexHull.Ns.last().x())*(mSDs[m].y() - ConvexHull.Ns.last().y()) - (mSDs[mem_m].y() - ConvexHull.Ns.last().y())*(mSDs[m].x() - ConvexHull.Ns.last().x());
			if (ccw>0){
				mem_m=m;
			}
		}
		ConvexHull.Ns.append(mSDs[mem_m]);  
		mSDs.removeAt(mem_m); 
	}
	ConvexHull.EraseSpecialPoints();
	ConvexHull.MakeCornersSpecial();
	ConvexHull.SortByType("CORNER"); 
	ConvexHull.AddPosition();
	ConvexHull.RefineByLength(this->size);
}

/*! \brief
	Merges the first and the last point of intersections to belonging convexhull.

	If the first or last point of the intersection is closer than 1e-12 to one special point of the convexhull,
	then the point gets the coordinates of this convexhull point.

	If the first or last point of the intersection is close to a segment of the convexhull (smaller than 1e-12) and not close to one special point of the convex hull,
	then this point is added to the convex hull as a special point.
*/
void C_Surface::alignIntersectionsToConvexHull()
{
	C_Vector3D * xProjection = new C_Vector3D;
	for (int i = 0; i != Intersections.length(); i++)
	{
		for (int n=0;n!=ConvexHull.Ns.length()-1;n++){
			if (ConvexHull.Ns[n].type() != "DEFAULT" && lengthSquared(ConvexHull.Ns[n] - Intersections[i]->Ns.first())<1e-24){
				Intersections[i]->Ns.first() = ConvexHull.Ns[n];
				break;
			}
			else if (ConvexHull.Ns[n + 1].type() != "DEFAULT" && lengthSquared(ConvexHull.Ns[n + 1] - Intersections[i]->Ns.first())<1e-24){
				Intersections[i]->Ns.first() = ConvexHull.Ns[n + 1];
				break;
			}
			else{
				projectTo(Intersections[i]->Ns.first(), ConvexHull.Ns[n], ConvexHull.Ns[n + 1], xProjection);
				if (lengthSquared(*xProjection) != 0 && lengthSquared(*xProjection - Intersections[i]->Ns.first())<1e-24){
					Intersections[i]->Ns.first().setType("COMMON_INTERSECTION_CONVEXHULL_POINT");
					ConvexHull.Ns.insert(n + 1, Intersections[i]->Ns.first());
					break;
				}
			}
		}
		for (int n = 0; n != ConvexHull.Ns.length() - 1; n++){
			if (ConvexHull.Ns[n].type() != "DEFAULT" && lengthSquared(ConvexHull.Ns[n] - Intersections[i]->Ns.last())<1e-24){
				Intersections[i]->Ns.last() = ConvexHull.Ns[n];
				break;
			}
			else if (ConvexHull.Ns[n + 1].type() != "DEFAULT" && lengthSquared(ConvexHull.Ns[n + 1] - Intersections[i]->Ns.last())<1e-24){
				Intersections[i]->Ns.last() = ConvexHull.Ns[n + 1];
				break;
			}
			else{
				projectTo(Intersections[i]->Ns.last(), ConvexHull.Ns[n], ConvexHull.Ns[n + 1], xProjection);
				if (lengthSquared(*xProjection) != 0 && lengthSquared(*xProjection - Intersections[i]->Ns.last())<1e-24){
					Intersections[i]->Ns.last().setType("COMMON_INTERSECTION_CONVEXHULL_POINT");
					ConvexHull.Ns.insert(n + 1, Intersections[i]->Ns.last());
					break;
				}
			}
		}
	}
	ConvexHull.CleanIdenticalPoints();
	ConvexHull.AddPosition();
	ConvexHull.RefineByLength(this->size);
}

void C_Surface::calculate_triangles(bool withConstraints, double gradient){
	struct triangulateio in, out;
	int s_start,s_end;
	bool found_1,found_2;
	QList<C_Vector3D> points;	// Pointlist storing points of ConvexHull and Intersections //
	QList<int> segments;		// Segmentlist storing segments of ConvexHull and Intersections //
	QList<double> refValues;

	//start_index
	// adding segments of the convex hull and intersections as constraints
	if (withConstraints){
		for (int s=0;s!=this->Constraints.length();s++){
			if (Constraints[s].Ns.length()==1 && Constraints[s].Type!="UNDEFINED"){
				points.append(Constraints[s].Ns);
				refValues.append(Constraints[s].size);
			}
			if (Constraints[s].Ns.length()>1 && Constraints[s].Type!="UNDEFINED"){
				for (int n=0;n!=Constraints[s].Ns.length()-1;n++){
					found_1=false;
					for (int p=0;p!=points.length();p++){
						if (lengthSquared(Constraints[s].Ns[n] - points[p])<1e-24){
							found_1=true;
							s_start=p;
						}
					}
					if (!found_1){
						s_start=points.length();
						points.append(Constraints[s].Ns[n]);
						refValues.append(Constraints[s].size);
					}

					found_2=false;
					for (int p=0;p!=points.length();p++){
						if (lengthSquared(Constraints[s].Ns[n + 1] - points[p])<1e-24){
							found_2=true;
							s_end=p;
						}
					}
					if (!found_2){
						s_end=points.length();
						points.append(Constraints[s].Ns[n+1]);
						refValues.append(Constraints[s].size);
					}
					segments.append(s_start); 
					segments.append(s_end); 
				}
			}
		}
	}else{
		for (int n=0;n!=ConvexHull.Ns.length()-1;n++){
			points.append(ConvexHull.Ns[n]);
			refValues.append(this->size);
			segments.append(n);
			if (n!=ConvexHull.Ns.length()-2){
				segments.append(n+1);
			}else{
				segments.append(0);
			}
		}
	}
	//end_index

	in.numberofpoints = points.length();
	in.pointlist = new REAL[in.numberofpoints * 2];
	in.pointmarkerlist = new int[in.numberofpoints];
	for (int p=0;p<points.length();p++){
		in.pointlist[2*p] = points[p].x();
		in.pointlist[2*p+1] = points[p].y();
		in.pointmarkerlist[p]=0;
	}
	in.numberofpointattributes = 0;
	in.pointattributelist = (REAL *) NULL;
/*
 *	Fill in the segments in the struct in
 */	
	in.numberofsegments = segments.length()/2;
	in.segmentlist = new int[in.numberofsegments * 2];
	in.segmentmarkerlist = new int[in.numberofsegments];
	for (int s=0;s<segments.length()/2;s++){ 
		in.segmentlist[2*s]=segments[2*s];
		in.segmentlist[2*s+1]=segments[2*s+1];
		in.segmentmarkerlist[s]=0;
	}
/*
 *	Fill in the Holes in the struct in
 */
	in.numberofholes = this->HoleCoords.length();
	in.holelist = new REAL[in.numberofholes * 2];
	for (int h=0;h<HoleCoords.length();h++){  
		in.holelist[2*h] = HoleCoords[h].x();
		in.holelist[2*h+1] = HoleCoords[h].y();
	}
/*
 *	Fill other parameters (not need to be specified) in the struct in
 */
	in.numberofregions = 0;
	in.regionlist = (REAL *) NULL;
	in.triangleattributelist = (REAL *) NULL;
	in.trianglearealist = (REAL *) NULL;

	in.meshsize = this->size;
	in.refinesize = new REAL[in.numberofpoints];
	for (int p = 0; p<points.length(); p++){
		in.refinesize[p] = refValues[p];
	}

	

		/*
 *	out struct 
 */	  
	out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
  /* Not needed if -N switch used or number of point attributes is zero: */
	out.pointattributelist = (REAL *) NULL;
	out.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
	out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
  /* Not needed if -E switch used or number of triangle attributes is zero: */
	out.triangleattributelist = (REAL *) NULL;
	out.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
/* Needed only if segments are output (-p or -c) and -P not used: */
	out.segmentlist = (int *) NULL;
  /* Needed only if segments are output (-p or -c) and -P and -B not used: */
	out.segmentmarkerlist = (int *) NULL;
	out.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
	out.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */

	this->Ns.clear(); 
	this->Ts.clear(); 

	GradientControl& gc = GradientControl::getInstance();
	gc.update(gradient, in.meshsize, in.numberofpoints, in.pointlist,
	          in.refinesize);

	if (in.numberofsegments>2){
		triangulate(const_cast<char*>("pzYYu"), &in, &out, 0);

		if (out.numberofpoints>0){
			for (int i = 0; i < out.numberofpoints; i++){
				C_Vector3D No;
				No.triID=i; 
				No.setX(out.pointlist[2*i]);
				No.setY(out.pointlist[2*i+1]);
				Ns.append(No); 
			}
		}
		if (out.numberoftriangles>0){
			for (int t=0;t!=out.numberoftriangles;t++){
				C_Triangle Tr;
				Tr.Ns[0]=&Ns[out.trianglelist[t*3]];
				Tr.Ns[1]=&Ns[out.trianglelist[t*3+1]];
				Tr.Ns[2]=&Ns[out.trianglelist[t*3+2]];
				Ts.append(Tr); 
			}
		}
	}
	this->duplicates = points.length();
	points.clear();
	segments.clear();
	free(in.pointlist);
	free(in.pointattributelist);
	free(in.pointmarkerlist);
	free(in.segmentlist);
	free(in.segmentmarkerlist);
	free(in.regionlist);
	free(out.pointlist);
	free(out.pointattributelist);
	free(out.trianglelist);
	free(out.triangleattributelist);
}

C_Surface::C_Surface(){
	this->drawScatteredData = false;
	this->drawConvexHull = false;
	this->drawFaces = false;
	this->drawEdges = false;
	this->drawIntEdges = false;
	this->drawIntVertices = false;
	this->drawMatFaces = false;
	this->drawMatEdges = false;
	this->drawConstraints=false;
	this->size = 0;
}

C_Surface::~C_Surface(){
}

void quickSort(double x[], double y[], double z[], int left, int right) {
int i = left, j = right;
double tmpX,tmpY,tmpZ;
double pivot = x[(left + right) / 2];
/* partition */
	while (i <= j){
		while (x[i] < pivot)
			i++;
			while (x[j] > pivot)
				j--;
			if (i <= j){
				tmpX = x[i];
				tmpY = y[i];
				tmpZ = z[i];
				x[i]=x[j];
				y[i]=y[j];
				z[i]=z[j];
				x[j]=tmpX;
				y[j]=tmpY;
				z[j]=tmpZ;
				i++;
				j--;
			}
		};
/* recursion */
		if (left < j) quickSort(x, y, z, left, j);
		if (i < right) quickSort(x, y, z, i, right);
}

void C_Surface::clearScatteredData(){
	int dim = this->SDs.length();

	double * x = new double[dim];
	double * y = new double[dim];
	double * z = new double[dim];

	for (int s=0;s!=dim;s++){
		x[s]=this->SDs[s].x();
		y[s]=this->SDs[s].y();
		z[s]=this->SDs[s].z();
	}
	quickSort(x,y,z,0,dim-1);
	for (int s1=0;s1!=dim-1;s1++){
		if (x[s1]==x[s1+1]){
			if (y[s1]==y[s1+1]){
				if (z[s1]==z[s1+1]){
					for (int s2=0;s2!=this->SDs.length();s2++){
						if (x[s1]==this->SDs[s2].x()){
							if (y[s1]==this->SDs[s2].y()){
								if (z[s1]==this->SDs[s2].z()){
									this->SDs.removeAt(s2);
									break;
								}
							}
						}
					}
				}
			}
		}
	}

	delete[] x;
	delete[] y;
	delete[] z;
}

void C_Surface::makeScatteredData(){
    GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);
	glDisable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.White);

	if (Type=="UNIT") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.LightBlue);
	if (Type=="FAULT") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.LightRed);
	if (Type=="BORDER") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.LightGreen);
	glBegin(GL_POINTS);
    for (int s = 0; s!=SDs.length();s++){
		glVertex3f(SDs[s].x(),SDs[s].y(),SDs[s].z());
	}
    glEnd();

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.Grey);
	glEnable(GL_LIGHT0);
    glEndList();
	listScatteredData = list;
}

void C_Surface::makeConvexHull(){
    GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);
	glDisable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.White);

	glBegin(GL_LINE_STRIP);
	for (int n=0;n!=ConvexHull.Ns.length();n++){
		if (ConvexHull.Ns[n].type()=="DEFAULT"){
			if (Type=="UNIT") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Blue);
			if (Type=="FAULT") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Red);
			if (Type=="BORDER") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Green);
		}
		glVertex3f(ConvexHull.Ns[n].x(),ConvexHull.Ns[n].y(),ConvexHull.Ns[n].z());
	}
	glEnd();
	glBegin(GL_POINTS);
	for (int n=0;n!=ConvexHull.Ns.length();n++){
		if (ConvexHull.Ns[n].type()=="DEFAULT"){
			if (Type=="UNIT") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Blue);
			if (Type=="FAULT") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Red);
			if (Type=="BORDER") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Green);
		}else{
			if (Type=="UNIT") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.LightBlue);
			if (Type=="FAULT") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.LightRed);
			if (Type=="BORDER") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.LightGreen);
		}
		glVertex3f(ConvexHull.Ns[n].x(),ConvexHull.Ns[n].y(),ConvexHull.Ns[n].z());
	}
	glEnd();

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.Grey);
	glEnable(GL_LIGHT0);
	glEndList();
	listConvexHull=list;
}

void C_Surface::makeFaces(){
	GLuint list = glGenLists(1);
	glNewList(list, GL_COMPILE);

	if (Type=="UNIT") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.LightBlue);
	if (Type=="FAULT") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.LightRed);
	if (Type=="BORDER") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.LightGreen);
	glBegin(GL_TRIANGLES);
    for (int t = 0; t!=Ts.length();t++){
		glNormal3d(Ts[t].normal_vector.x(),Ts[t].normal_vector.y(),Ts[t].normal_vector.z());
		glVertex3d(Ts[t].Ns[0]->x(),Ts[t].Ns[0]->y(),Ts[t].Ns[0]->z());
		glVertex3d(Ts[t].Ns[1]->x(),Ts[t].Ns[1]->y(),Ts[t].Ns[1]->z());
		glVertex3d(Ts[t].Ns[2]->x(),Ts[t].Ns[2]->y(),Ts[t].Ns[2]->z());
	}
    glEnd();

    glEndList();
	listFaces = list;
}

void C_Surface::makeEdges(){
    GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);
	glDisable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.White);

	if (Type=="UNIT") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.LightBlue);
	if (Type=="FAULT") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.LightRed);
	if (Type=="BORDER") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.LightGreen);
    for (int t = 0; t!=Ts.length();t++){
		glBegin(GL_LINE_LOOP);
		glVertex3d(Ts[t].Ns[0]->x(),Ts[t].Ns[0]->y(),Ts[t].Ns[0]->z());
		glVertex3d(Ts[t].Ns[1]->x(),Ts[t].Ns[1]->y(),Ts[t].Ns[1]->z());
		glVertex3d(Ts[t].Ns[2]->x(),Ts[t].Ns[2]->y(),Ts[t].Ns[2]->z());
	    glEnd();
	}

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.Grey);
	glEnable(GL_LIGHT0);
    glEndList();
	listEdges = list;
}

void C_Surface::makeIntEdges(){
    GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);
	glDisable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.White);

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,this->Cols.White);
    for (int i = 0; i!=Intersections.length();i++){
		glBegin(GL_LINE_STRIP);
	    for (int n = 0; n!=Intersections[i]->Ns.length();n++){
			glVertex3d(Intersections[i]->Ns[n].x(),Intersections[i]->Ns[n].y(),Intersections[i]->Ns[n].z());
		}
	    glEnd();
	}

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.Grey);
	glEnable(GL_LIGHT0);
    glEndList();
	listIntEdges = list;
}

void C_Surface::makeIntVertices(){
    GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);
	glDisable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.White);

	glBegin(GL_POINTS);
    for (int i = 0; i!=Intersections.length();i++){
		if (Intersections[i]->Ns.length()==1){
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,this->Cols.White);
			glVertex3d(Intersections[i]->Ns[0].x(),Intersections[i]->Ns[0].y(),Intersections[i]->Ns[0].z());
		}else{
			for (int n = 0; n!=Intersections[i]->Ns.length();n++){
				if (n==0) {
					glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,this->Cols.Red);
				}else if (n==(Intersections[i]->Ns.length()-1)){
					glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,this->Cols.Blue);
				}else if(Intersections[i]->Ns[n].type()=="TRIPLE_POINT"){
					glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,this->Cols.Green);
				}else{
					glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,this->Cols.White);
				}
				glVertex3d(Intersections[i]->Ns[n].x(),Intersections[i]->Ns[n].y(),Intersections[i]->Ns[n].z());
			}
		}
	}
    glEnd();

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.Grey);
	glEnable(GL_LIGHT0);
    glEndList();
	listIntVertices = list;
}

void C_Surface::makeConstraints(bool selectionMode){
    GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);
	glDisable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.White);

	/* In selection mode, set a uniform width to avoid overlapping of points
	 * and lines. */
	GLfloat previousPointWidth;
	GLfloat previousLineWidth;
	GLfloat desiredWidth = 2.0f;
	if( selectionMode ) {
		glGetFloatv(GL_POINT_SIZE, &previousPointWidth);
		glGetFloatv(GL_LINE_WIDTH, &previousLineWidth);
		glPointSize(desiredWidth);
		glLineWidth(desiredWidth);
	}

	for (int s = 0; s!=Constraints.length();s++){
		if (Constraints[s].Type=="SEGMENTS") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,this->Cols.Blue);
		if (Constraints[s].Type=="HOLES") glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,this->Cols.Red);
		if (Constraints[s].Type=="UNDEFINED"){
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,this->Cols.White);

			/* In selection mode, draw constraints with a solid line to allow
			 * the paint fill bucket tool to work correctly. */
			if( ! selectionMode )
				glEnable(GL_LINE_STIPPLE);
		}
		glColor3f(Constraints[s].RGB[0]/255.0f, Constraints[s].RGB[1]/255.0f, Constraints[s].RGB[2]/255.0f);

		glBegin(GL_LINE_STRIP);
			for (int n=0;n!=Constraints[s].Ns.length();n++){
				glVertex3f(Constraints[s].Ns[n].x(),Constraints[s].Ns[n].y(),Constraints[s].Ns[n].z());
			}
		glEnd();

		glBegin(GL_POINTS);
		if (Constraints[s].Ns.length()==1){
			glVertex3f(Constraints[s].Ns.first().x(),Constraints[s].Ns.first().y(),Constraints[s].Ns.first().z());
		}
		else{
			glVertex3f(Constraints[s].Ns.first().x(),Constraints[s].Ns.first().y(),Constraints[s].Ns.first().z());
			glVertex3f(Constraints[s].Ns.last().x(),Constraints[s].Ns.last().y(),Constraints[s].Ns.last().z());
		}
		glEnd();

		glColor3f(0.0f, 0.0f, 0.0f);

		/* Take selection mode into account. */
		if( ! selectionMode )
			glDisable (GL_LINE_STIPPLE);
	}

	/* Restore previous point and line width. */
	if( selectionMode ) {
		glPointSize(previousPointWidth);
		glLineWidth(previousLineWidth);
	}

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.Grey);
	glEnable(GL_LIGHT0);
	glEndList();
	listConstraints = list;
}

void C_Surface::makeVTU_SD(){
	this->VTU.clear();  
	this->VTU.NumberOfPoints=this->SDs.length();
	this->VTU.NumberOfCells=1;
	this->VTU.Points=this->SDs;
	for (int s=0;s!=this->SDs.length();s++) this->VTU.connectivity.append(s);
	this->VTU.offsets.append(this->SDs.length());
	this->VTU.types.append(2);  
}

void C_Surface::makeVTU_CH(){
	this->VTU.clear();  
	this->VTU.NumberOfPoints=this->ConvexHull.Ns.length();  
	this->VTU.NumberOfCells=1;
	this->VTU.Points.append(this->ConvexHull.Ns);
	for (int n=0;n!=this->ConvexHull.Ns.length();n++) this->VTU.connectivity.append(n);
	this->VTU.offsets.append(this->ConvexHull.Ns.length());
	this->VTU.types.append(4);  
	for (int n=0;n!=this->ConvexHull.Ns.length();n++){
		if (this->ConvexHull.Ns[n].type()=="DEFAULT") this->VTU.pointData.append(0);  
		if (this->ConvexHull.Ns[n].type()=="CORNER") this->VTU.pointData.append(1); 
		if (this->ConvexHull.Ns[n].type()=="MERGED_CONVEXHULL_POINT") this->VTU.pointData.append(2); 
		if (this->ConvexHull.Ns[n].type()=="COMMON_INTERSECTION_CONVEXHULL_POINT") this->VTU.pointData.append(3); 
	}
}

void C_Surface::makeVTU_TRI(){
	this->VTU.clear();  
	this->VTU.NumberOfPoints=this->Ns.length();  
	this->VTU.NumberOfCells=this->Ts.length();
	this->VTU.Points.append(this->Ns);
	for (int t=0;t!=this->Ts.length();t++){
		this->VTU.connectivity.append(this->Ts[t].Ns[0]->triID);  
		this->VTU.connectivity.append(this->Ts[t].Ns[1]->triID);  
		this->VTU.connectivity.append(this->Ts[t].Ns[2]->triID);  
		this->VTU.offsets.append(t*3+3);
		this->VTU.types.append(5);  
	}
}

void C_Surface::makeVTU_CON(){
	this->VTU.clear();
	for (int s=0; s!=this->Constraints.length();s++){
		this->VTU.NumberOfPoints+=this->Constraints[s].Ns.length(); 
		this->VTU.NumberOfCells++;
		this->VTU.Points.append(Constraints[s].Ns);
		for (int n=0;n!=this->Constraints[s].Ns.length();n++) this->VTU.connectivity.append(this->VTU.connectivity.length());
		if (this->VTU.offsets.length()==0){
			this->VTU.offsets.append(this->Constraints[s].Ns.length());
		}else{
			this->VTU.offsets.append(this->Constraints[s].Ns.length()+this->VTU.offsets.last());
		}
		this->VTU.types.append(4);
		if (this->Constraints[s].Type=="SEGMENTS") this->VTU.matType.append(0);
		if (this->Constraints[s].Type=="UNDEFINED") this->VTU.matType.append(1);
		if (this->Constraints[s].Type=="HOLES") this->VTU.matType.append(2);
		this->VTU.matR.append(this->Constraints[s].RGB[0]);     
		this->VTU.matG.append(this->Constraints[s].RGB[1]);     
		this->VTU.matB.append(this->Constraints[s].RGB[2]);     
	}
}

C_Polyline *C_Model::findPolyline(const QString &name)
{
	for (int i = 0; i != Polylines.length(); i++)
		if( Polylines[i].Name == name )
			return &Polylines[i];
	return NULL;
}

C_Surface *C_Model::findSurface(const QString &name)
{
	for (int i = 0; i != Surfaces.length(); i++)
		if( Surfaces[i].Name == name )
			return &Surfaces[i];
	return NULL;
}

void C_Model::makeVTU_INT(){
	this->VTU.clear();  
	for (int i=0; i!=this->Intersections.length();i++){
		this->VTU.NumberOfPoints+=this->Intersections[i].Ns.length(); 
		this->VTU.NumberOfCells++;
		this->VTU.Points.append(this->Intersections[i].Ns);
		for (int n=0;n!=this->Intersections[i].Ns.length();n++) this->VTU.connectivity.append(this->VTU.connectivity.length());
		if (this->VTU.offsets.length()==0){
			this->VTU.offsets.append(this->Intersections[i].Ns.length());
		}else{
			this->VTU.offsets.append(this->Intersections[i].Ns.length()+this->VTU.offsets.last());
		}
		if (this->Intersections[i].Ns.length()==1){
			VTU.types.append(1);
		}else{
			VTU.types.append(4);
		}
		this->VTU.Object1.append(this->Intersections[i].Object[0]);     
		this->VTU.Object2.append(this->Intersections[i].Object[1]);     
	}
	for (int p=0;p!=this->VTU.Points.length();p++){
		if (this->VTU.Points[p].type()=="DEFAULT") this->VTU.pointData.append(0);  
		if (this->VTU.Points[p].type()=="TRIPLE_POINT") this->VTU.pointData.append(1); //intersction line cuts intersection line
		if (this->VTU.Points[p].type()=="INTERSECTION_POINT") this->VTU.pointData.append(2); //well cuts surface  
	}
}

void C_Model::makeVTU_MAT(){
	this->VTU.clear();  
	for (int m=0;m!=this->Mats.length();m++){ 
		this->VTU.NumberOfPoints+=this->Mats[m].Locations.length(); 
		this->VTU.NumberOfCells+=this->Mats[m].Locations.length();
		this->VTU.Points.append(this->Mats[m].Locations);
	}
	int index=0;
	for (int m=0;m!=this->Mats.length();m++){ 
		for (int l=0;l!=this->Mats[m].Locations.length();l++){ 
			this->VTU.connectivity.append(index);
			this->VTU.offsets.append(++index);
			this->VTU.types.append(1);
			this->VTU.matType.append(m);
		}
	}
}

void C_Model::makeVTU_TET(){
	this->VTU.clear();  
	this->VTU.NumberOfPoints=this->Mesh->numberofpoints;  
	this->VTU.NumberOfCells=this->Mesh->numberofedges; 
	this->VTU.NumberOfCells+=this->Mesh->numberoftriangles;
	this->VTU.NumberOfCells+=this->Mesh->numberoftetrahedra;
	
	for (int p=0;p!=this->Mesh->numberofpoints;p++){
		this->VTU.Points.append(C_Vector3D(this->Mesh->pointlist[p*3+0],this->Mesh->pointlist[p*3+1],this->Mesh->pointlist[p*3+2]));
	}
	for (int t=0;t!=this->Mesh->numberoftetrahedra;t++){   
		this->VTU.connectivity.append(this->Mesh->tetrahedronlist[t*4+0]);
		this->VTU.connectivity.append(this->Mesh->tetrahedronlist[t*4+1]);
		this->VTU.connectivity.append(this->Mesh->tetrahedronlist[t*4+2]);
		this->VTU.connectivity.append(this->Mesh->tetrahedronlist[t*4+3]);
		this->VTU.offsets.append(t*4+4);
		this->VTU.types.append(10);  
		this->VTU.matType.append(this->Mesh->tetrahedronmarkerlist[t]);
	}
	for (int f=0;f!=this->Mesh->numberoftriangles;f++){  
		this->VTU.connectivity.append(this->Mesh->trianglelist[f*3+0]); 
		this->VTU.connectivity.append(this->Mesh->trianglelist[f*3+1]); 
		this->VTU.connectivity.append(this->Mesh->trianglelist[f*3+2]); 
		this->VTU.offsets.append(this->Mesh->numberoftetrahedra*4+f*3+3);  
		this->VTU.types.append(5);  
		this->VTU.matType.append(this->Mesh->trianglemarkerlist[f]); 
	}
	for (int e=0;e!=this->Mesh->numberofedges;e++){ 
		this->VTU.connectivity.append(this->Mesh->edgelist[e*2+0]);  
		this->VTU.connectivity.append(this->Mesh->edgelist[e*2+1]);  
		this->VTU.offsets.append(this->Mesh->numberoftetrahedra*4+this->Mesh->numberoftriangles*3+e*2+2);  
		this->VTU.types.append(3);  
		this->VTU.matType.append(this->Mesh->edgemarkerlist[e]);
	}
}

void C_Model::VTU_SD_to_Polyline(QString name, QString type, int mat, double size){
	C_Polyline tmpPolyline;
	tmpPolyline.Type=type; 
	tmpPolyline.Name=name.section("_SD",0,0);
	tmpPolyline.MaterialID=mat; 
	tmpPolyline.size = size;
	tmpPolyline.SDs=this->VTU.Points;
	this->Polylines.append(tmpPolyline); 
	this->VTU.clear();
}

void C_Model::VTU_SD_to_Surface(QString name, QString type, int mat, double size){
	C_Surface tmpSurface;
	tmpSurface.Type=type; 
	tmpSurface.Name=name.section("_SD",0,0);
	tmpSurface.MaterialID = mat;
	tmpSurface.size = size;
	tmpSurface.SDs=this->VTU.Points;
	this->Surfaces.append(tmpSurface);  
	this->VTU.clear();
}

void C_Model::VTU_CH_to_Surface(QString name, QString type){
	for (int p=0;p!=this->VTU.pointData.length();p++){
		if (this->VTU.pointData[p]==0) this->VTU.Points[p].setType("DEFAULT");
		if (this->VTU.pointData[p]==1) this->VTU.Points[p].setType("CORNER");
		if (this->VTU.pointData[p]==2) this->VTU.Points[p].setType("MERGED_CONVEXHULL_POINT");
		if (this->VTU.pointData[p]==3) this->VTU.Points[p].setType("COMMON_INTERSECTION_CONVEXHULL_POINT");
	}
	for (int s=0;s!=this->Surfaces.length();s++){
		if (this->Surfaces[s].Name==name.section("_CH",0,0) && this->Surfaces[s].Type==type){
			Surfaces[s].ConvexHull.Ns=this->VTU.Points;
		}
	}
	this->VTU.clear();
}

void C_Model::VTU_SEG_to_Polyline(QString name, QString type){
	for (int p=0;p!=this->Polylines.length();p++){
		if (this->Polylines[p].Name==name.section("_SEG",0,0) && this->Polylines[p].Type==type){
			this->Polylines[p].Path.Ns=this->VTU.Points;
		}
	}
	this->VTU.clear();
}

void C_Model::VTU_TRI_to_Surface(QString name, QString type){
	for (int s=0;s!=this->Surfaces.length();s++){
		if (this->Surfaces[s].Name==name.section("_TRI",0,0) && this->Surfaces[s].Type==type){
			this->Surfaces[s].Ns=this->VTU.Points;
			for (int t=0;t!=this->VTU.NumberOfCells;t++){
				C_Triangle To;
				To.Ns[0]=&this->Surfaces[s].Ns[this->VTU.connectivity[t*3+0]];
				To.Ns[0]->triID=this->VTU.connectivity[t*3+0];
				To.Ns[1]=&this->Surfaces[s].Ns[this->VTU.connectivity[t*3+1]];
				To.Ns[1]->triID=this->VTU.connectivity[t*3+1];
				To.Ns[2]=&this->Surfaces[s].Ns[this->VTU.connectivity[t*3+2]];
				To.Ns[2]->triID=this->VTU.connectivity[t*3+2];
				To.setNormalVector(); 
				this->Surfaces[s].Ts.append(To);   
			}
		}
	}
	this->VTU.clear();
}

void C_Model::VTU_CON_to_Polyline(QString name, QString type){
	for (int p=0;p!=this->Polylines.length();p++){
		if (this->Polylines[p].Name==name.section("_CON",0,0) && this->Polylines[p].Type==type){
			for (int n=0;n!=this->VTU.NumberOfPoints;n++){
				C_Line Constraint;
				Constraint.RGB[0]=this->VTU.matR[n]; 
				Constraint.RGB[1]=this->VTU.matG[n]; 
				Constraint.RGB[2]=this->VTU.matB[n];
				if (this->VTU.matType[n]==0){
					Constraint.Type="SEGMENTS";
					this->VTU.Points[n].setType("INTERSECTION_POINT");  
				}
				if (this->VTU.matType[n]==1){
					Constraint.Type="UNDEFINED";
					this->VTU.Points[n].setType("DEFAULT");  
				}
				Constraint.Ns.clear(); 
				Constraint.Ns.append(this->VTU.Points[n]);
				this->Polylines[p].Constraints.append(Constraint);   
			}
		}
	}
	this->VTU.clear();
}

void C_Model::VTU_CON_to_Surface(QString name, QString type){
	int pos=0;
	for (int s=0;s!=this->Surfaces.length();s++){
		if (this->Surfaces[s].Name==name.section("_CON",0,0) && this->Surfaces[s].Type==type){
			for (int c=0;c!=this->VTU.NumberOfCells;c++){
				C_Line Constraint;
				Constraint.RGB[0]=this->VTU.matR[c]; 
				Constraint.RGB[1]=this->VTU.matG[c]; 
				Constraint.RGB[2]=this->VTU.matB[c]; 
				if (this->VTU.matType[c] == 0) Constraint.Type = "SEGMENTS";
				if (this->VTU.matType[c]==1) Constraint.Type="UNDEFINED";
				if (this->VTU.matType[c]==2) Constraint.Type="HOLES";
				Constraint.Ns=this->VTU.Points.mid(pos,this->VTU.offsets[c]-pos);
				pos=this->VTU.offsets[c];
				this->Surfaces[s].Constraints.append(Constraint);   
			}
		}
	}
	this->VTU.clear();
}

void C_Model::VTU_INT_to_Model(){
	int pos=0;
	for (int p=0;p!=this->VTU.Points.length();p++){
		if (this->VTU.pointData[p]==0) this->VTU.Points[p].setType("DEFAULT");
		if (this->VTU.pointData[p]==1) this->VTU.Points[p].setType("TRIPLE_POINT");
		if (this->VTU.pointData[p]==2) this->VTU.Points[p].setType("INTERSECTION_POINT");
	}
	for (int c=0;c!=this->VTU.NumberOfCells;c++){
		C_Line Intersection;
		Intersection.Ns=this->VTU.Points.mid(pos,this->VTU.offsets[c]-pos);
		Intersection.Object[0]=this->VTU.Object1[c];
		Intersection.Object[1]=this->VTU.Object2[c];
		pos=this->VTU.offsets[c];
		this->Intersections.append(Intersection);
		this->Surfaces[this->VTU.Object1[c]].Intersections.append(&this->Intersections.last());   
		if (this->VTU.types[c]==1){
			this->Polylines[this->VTU.Object2[c]].Intersections.append(&this->Intersections.last());   
		}else{
			this->Surfaces[this->VTU.Object2[c]].Intersections.append(&this->Intersections.last());   
		}
	}
	this->VTU.clear();
}

void C_Model::VTU_MAT_to_Model(){
	this->Mats.clear();
	C_Material * Ma;
	for (int m=0;m!=this->VTU.matType.length();m++){
		while (this->VTU.matType[m]>=this->Mats.length()){
			Ma = new C_Material;
			Mats.append(*Ma);
		}
		Mats[VTU.matType[m]].Locations.append(this->VTU.Points[m]);
	}
	this->VTU.clear();
}

void C_Model::VTU_TET_to_Model(){
	if (!this->Mesh) this->Mesh = new C_Mesh3D;
	for (int p=0;p!=this->VTU.NumberOfPoints;p++){
		this->Mesh->pointlist.append(VTU.Points[p].x());   
		this->Mesh->pointlist.append(VTU.Points[p].y());   
		this->Mesh->pointlist.append(VTU.Points[p].z());
		this->Mesh->numberofpoints++; 
	}

	for (int c=0;c!=this->VTU.NumberOfCells;c++){
		if (this->VTU.types[c]==3){
			this->Mesh->edgelist.append(this->VTU.connectivity[this->VTU.offsets[c]-2]);
			this->Mesh->edgelist.append(this->VTU.connectivity[this->VTU.offsets[c]-1]);
			this->Mesh->edgemarkerlist.append(this->VTU.matType[c]); 
			this->Mesh->numberofedges++;
		}
		if (this->VTU.types[c]==5){
			this->Mesh->trianglelist.append(this->VTU.connectivity[this->VTU.offsets[c]-3]);
			this->Mesh->trianglelist.append(this->VTU.connectivity[this->VTU.offsets[c]-2]);
			this->Mesh->trianglelist.append(this->VTU.connectivity[this->VTU.offsets[c]-1]);
			this->Mesh->trianglemarkerlist.append(this->VTU.matType[c]); 
			this->Mesh->numberoftriangles++;
		}
		if (this->VTU.types[c]==10){
			this->Mesh->tetrahedronlist.append(this->VTU.connectivity[this->VTU.offsets[c]-4]);
			this->Mesh->tetrahedronlist.append(this->VTU.connectivity[this->VTU.offsets[c]-3]);
			this->Mesh->tetrahedronlist.append(this->VTU.connectivity[this->VTU.offsets[c]-2]);
			this->Mesh->tetrahedronlist.append(this->VTU.connectivity[this->VTU.offsets[c]-1]);
			this->Mesh->tetrahedronmarkerlist.append(this->VTU.matType[c]); 
			this->Mesh->numberoftetrahedra++;
		}
	}
	this->VTU.clear();
}

void C_Surface::calculate_Constraints(){
	this->Constraints.clear(); 
	C_Line Constraint;
	int lastPos;
	lastPos=0;

	Constraint.RGB[0]=0; 
	Constraint.RGB[1]=0; 
	Constraint.RGB[2]=0; 
	Constraint.Type="UNDEFINED";
//	Constraint.Type = "SEGMENTS";
	for (int n = 1; n != this->ConvexHull.Ns.length(); n++){
		if (ConvexHull.Ns[n].type()!="DEFAULT" || n==ConvexHull.Ns.length()-1){
			Constraint.Ns=ConvexHull.Ns.mid(lastPos,n-lastPos+1);
			if (++Constraint.RGB[0]>255){Constraint.RGB[0]=0;if (++Constraint.RGB[1]>255){Constraint.RGB[1]=0;if (++Constraint.RGB[2]>255){}}}
			this->Constraints.append(Constraint);
			lastPos=n;
		}
	}
	for (int i=0;i!=this->Intersections.length();i++){
		if (this->Intersections[i]->Ns.length() == 1){
			Constraint.Ns=Intersections[i]->Ns;
			if (++Constraint.RGB[0]>255){Constraint.RGB[0]=0;if (++Constraint.RGB[1]>255){Constraint.RGB[1]=0;if (++Constraint.RGB[2]>255){}}}
			this->Constraints.append(Constraint);
		}else{
			lastPos=0;
			for (int n=1;n!=this->Intersections[i]->Ns.length();n++){
				if (Intersections[i]->Ns[n].type()!="DEFAULT" || n==Intersections[i]->Ns.length()-1){
					Constraint.Ns=Intersections[i]->Ns.mid(lastPos,n-lastPos+1);
					if (++Constraint.RGB[0]>255){Constraint.RGB[0]=0;if (++Constraint.RGB[1]>255){Constraint.RGB[1]=0;if (++Constraint.RGB[2]>255){}}}
					this->Constraints.append(Constraint);
					lastPos=n;
				}
			}
		}
	}
}

void C_Surface::separate_Constraints(){
	double x_sum,y_sum, area, area_sum;
	QList<C_Line> tmpHoles;
	C_Vector3D HoleCoord,Middle,Offset;
	bool extendable;
	for (int s=0;s!=this->Constraints.length();s++){
		if (this->Constraints[s].Type=="HOLES") tmpHoles.append(this->Constraints[s]);    
	}
	for (int h1=0;h1!=tmpHoles.length();h1++){
		extendable=true;
		while (extendable){
			extendable=false;
			for (int h2=h1+1;h2!=tmpHoles.length();h2++){
				if (lengthSquared(tmpHoles[h1].Ns.last() - tmpHoles[h2].Ns.first())<1e-24){
					for (int n=1;n!=tmpHoles[h2].Ns.length();n++) tmpHoles[h1].Ns.append(tmpHoles[h2].Ns[n]);
					tmpHoles.removeAt(h2--); 
					extendable=true;
					break;
				}
				if (lengthSquared(tmpHoles[h1].Ns.last() - tmpHoles[h2].Ns.last())<1e-24){
					tmpHoles[h2].Invert(); 
					for (int n=1;n!=tmpHoles[h2].Ns.length();n++) tmpHoles[h1].Ns.append(tmpHoles[h2].Ns[n]);
					tmpHoles.removeAt(h2--); 
					extendable=true;
					break;
				}
			}
		}
		if (lengthSquared(tmpHoles[h1].Ns.first() - tmpHoles[h1].Ns.last())>1e-24) tmpHoles.removeAt(h1--);
	}
    this->HoleCoords.clear();
    for (int t=0;t!=tmpHoles.length();t++){
        x_sum = 0.0, y_sum = 0.0, area = 0.0;
        for(int s=0;s!=tmpHoles[t].Ns.length()-1;s++){
            area_sum = (tmpHoles[t].Ns[s].x() * tmpHoles[t].Ns[s + 1].y()) - (tmpHoles[t].Ns[s+1].x() * tmpHoles[t].Ns[s].y());
            x_sum += area_sum * (tmpHoles[t].Ns[s].x() + tmpHoles[t].Ns[s+1].x());
            y_sum += area_sum * (tmpHoles[t].Ns[s].y() + tmpHoles[t].Ns[s+1].y());
            area += 0.5*area_sum;
        }
		HoleCoords.append(C_Vector3D(x_sum/(area*6),y_sum/(area*6),0.0));
    }
    tmpHoles.clear();
}

C_Model::C_Model(){
	this->shift = C_Vector3D(0,0,0);
	this->scale=1;
	this->ExportRotationAngle = 0.0;
	this->preMeshGradient = 2.0;
	this->meshGradient = 2.0;
}

C_Model::~C_Model(){
}

void C_Model::ReadGocadFile(){
	bool surfaceExists;
    QFile file(FileNameTmp);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
		return;
	QTextStream in(&file);
    QString line = in.readLine();
    while (!line.isNull()) {
		if (line.left(6).compare("*name:")==0){
			C_Surface tmpSurface;
			line = line.section(":",1,1);
			if (line.contains("/")){  
				line = line.section("/",1,1);
				tmpSurface.Type="UNIT";
			 }else{
				tmpSurface.Type="FAULT";
			 }
			 tmpSurface.Name=line; 
			 while (line.left(3).compare("END")!=0){
		         line = in.readLine();
				 if (line.left(4).compare("VRTX")==0){
					 C_Vector3D N;
					 N.setX(line.section(" ",2,2).toDouble());
					 N.setY(line.section(" ",3,3).toDouble());
					 N.setZ(line.section(" ",4,4).toDouble());
					 tmpSurface.SDs.append(N);  
				 }
			 }
			surfaceExists=false;
			for (int s=0;s!=this->Surfaces.length();s++){
				if (Surfaces[s].Name==tmpSurface.Name){
					Surfaces[s].SDs.append(tmpSurface.SDs);    
					surfaceExists=true;
				}
			}
			if (!surfaceExists){
				tmpSurface.MaterialID=-1; 
				this->Surfaces.append(tmpSurface);
			}
		 }
         line = in.readLine();
    }
	this->calculate_min_max();
	this->tranformForward();
}

void C_Model::AddSurface(QString type)
{
	this->tranformBackward();
	C_Surface tmpSurface;
	tmpSurface.Type=type;
	int NumberOfPoints;
	double d_value;
	QRegExp sep("(\\t+|\\s+|\\,|\\;)");
	QString line;
	QFile file(FileNameTmp);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
	QTextStream in(&file);
	QString ext=FileNameTmp.section(".",-1,-1);
	if (ext=="vtu")
	{
		while (!in.atEnd()){
			line=in.readLine().simplified();
			if (line.contains("<Piece"))
			{
				NumberOfPoints=line.section("\"",1,1).toInt();
			}
			if (line.contains("NumberOfComponents"))
			{
				for (int p=0;p!=NumberOfPoints;p++)
				{
					C_Vector3D SD;
					in>>d_value;
					SD.setX(d_value);
					in>>d_value;
					SD.setY(d_value);
					in>>d_value;
					SD.setZ(d_value);
					tmpSurface.SDs.append(SD);
				}
			}
		}
	}
	else
	{
		while (!in.atEnd())
		{
			line = in.readLine();
			line = line.trimmed(); // ignore spaces at the beginning and end
			if (line.isEmpty()) continue; // skip empty lines
			C_Vector3D SD;
			SD.setX(line.section(sep,0,0).toDouble());
			SD.setY(line.section(sep,1,1).toDouble());
			SD.setZ(line.section(sep,2,2).toDouble());
			tmpSurface.SDs.append(SD);
		}
	}
	FileNameTmp = FileNameTmp.section("/",-1,-1);
	FileNameTmp = FileNameTmp.section(".",0,0);
	tmpSurface.Name=FileNameTmp.section("_SD",0,0);
	tmpSurface.MaterialID=-1; 
	this->Surfaces.append(tmpSurface);
	this->calculate_min_max();
	this->tranformForward();
}

void C_Model::AddPolyline(QString type){
	this->tranformBackward();
	C_Polyline tmpPolyline;
	tmpPolyline.Type=type;
	int NumberOfPoints;
	double d_value;
	QRegExp sep("(\\t+|\\s+|\\,|\\;)");
	QString line;
	QFile file(FileNameTmp);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
	QTextStream in(&file);
	QString ext=FileNameTmp.section(".",-1,-1);
	if (ext=="vtu"){
		while (!in.atEnd()){
			line=in.readLine().simplified();
			if (line.contains("<Piece")){
				NumberOfPoints=line.section("\"",1,1).toInt();  
			}
			if (line.contains("NumberOfComponents")){
				for (int p=0;p!=NumberOfPoints;p++){
					C_Vector3D N;
					in>>d_value;
					N.setX(d_value);   
					in>>d_value;
					N.setY(d_value);   
					in>>d_value;
					N.setZ(d_value);   
					tmpPolyline.SDs.append(N); 
				}
			}
		}
	}else{
		while (!in.atEnd()) {
			line = in.readLine();
			C_Vector3D N;
			N.setX(line.section(sep,0,0).toDouble());
			N.setY(line.section(sep,1,1).toDouble());
			N.setZ(line.section(sep,2,2).toDouble());
			tmpPolyline.SDs.append(N);
		}
	}
	FileNameTmp = FileNameTmp.section("/",-1,-1);
	FileNameTmp = FileNameTmp.section(".",0,0);
	tmpPolyline.Name=FileNameTmp.section("_SD",0,0);
	tmpPolyline.MaterialID=-1; 
	Polylines.append(tmpPolyline);
	this->calculate_min_max();
	this->tranformForward(); 
}

void C_Model::DeleteSurface(QString surfaceName)
{
	for (int s = 0; s != this->Surfaces.length(); s++)
	{
		if (surfaceName == Surfaces[s].Name)
		{
			this->Surfaces.removeAt(s);
			return;
		}
	}
	for (int p = 0; p != this->Polylines.length(); p++)
	{
		if (surfaceName == Polylines[p].Name)
		{
			this->Polylines.removeAt(p);
			return;
		}
	}
}

void C_Model::calculateNumberWithMaterials(){
	this->Mesh->numberoftetrahedrawithmaterial = 0;
	for (int t = 0; t != this->Mesh->numberoftetrahedra; t++){
		if (this->getMaterial(10,t)!=-1) this->Mesh->numberoftetrahedrawithmaterial++;

	}
	this->Mesh->numberoftriangleswithmaterial = 0;
	for (int t = 0; t != this->Mesh->numberoftriangles; t++){
		if (this->getMaterial(5, t) != -1) this->Mesh->numberoftriangleswithmaterial++;

	}
	this->Mesh->numberofedgeswithmaterial = 0;
	for (int t = 0; t != this->Mesh->numberofedges; t++){
		if (this->getMaterial(3, t) != -1) this->Mesh->numberofedgeswithmaterial++;

	}

}

void C_Model::ExportLists_clean(){
	this->Mesh->numberofpoints_export = this->Mesh->numberofedges_export = this->Mesh->numberoftriangles_export = this->Mesh->numberoftetrahedra_export = 0;
	this->Mesh->pointlist_export.clear();
	this->Mesh->edgelist_export.clear();
	this->Mesh->edgemarkerlist_export.clear();
	this->Mesh->trianglelist_export.clear();
	this->Mesh->trianglemarkerlist_export.clear();
	this->Mesh->tetrahedronlist_export.clear();
	this->Mesh->tetrahedronmarkerlist_export.clear();
}

void C_Model::ExportLists_make(){
	this->ExportLists_clean();
	for (int m = 0; m != this->Mats.length(); m++){
		for (int t = 0; t != this->Mesh->tetrahedronmarkerlist.length(); t++){
			if (m == this->Mesh->tetrahedronmarkerlist[t]){
				this->Mesh->tetrahedronmarkerlist_export.append(m);
				this->Mesh->tetrahedronlist_export.append(this->Mesh->tetrahedronlist[t * 4 + 0]);
				this->Mesh->tetrahedronlist_export.append(this->Mesh->tetrahedronlist[t * 4 + 1]);
				this->Mesh->tetrahedronlist_export.append(this->Mesh->tetrahedronlist[t * 4 + 2]);
				this->Mesh->tetrahedronlist_export.append(this->Mesh->tetrahedronlist[t * 4 + 3]);
				this->Mesh->numberoftetrahedra_export++;
			}
		}
	}
	for (int s = 0; s != this->Surfaces.length(); s++){
		if (this->Surfaces[s].MaterialID != -1){
			for (int t = 0; t != this->Mesh->trianglemarkerlist.length(); t++){
				if (s==this->Mesh->trianglemarkerlist[t]){
					this->Mesh->trianglemarkerlist_export.append(this->Surfaces[s].MaterialID);
					this->Mesh->trianglelist_export.append(this->Mesh->trianglelist[t * 3 + 0]);
					this->Mesh->trianglelist_export.append(this->Mesh->trianglelist[t * 3 + 1]);
					this->Mesh->trianglelist_export.append(this->Mesh->trianglelist[t * 3 + 2]);
					this->Mesh->numberoftriangles_export++;
				}
			}
		}
	}
	for (int p = 0; p != this->Polylines.length(); p++){
		if (this->Polylines[p].MaterialID != -1){
			for (int e = 0; e != this->Mesh->edgemarkerlist.length(); e++){
				if (p == this->Mesh->edgemarkerlist[e]){
					this->Mesh->edgemarkerlist_export.append(this->Polylines[p].MaterialID);
					this->Mesh->edgelist_export.append(this->Mesh->edgelist[e * 2 + 0]);
					this->Mesh->edgelist_export.append(this->Mesh->edgelist[e * 2 + 1]);
					this->Mesh->numberofedges_export++;
				}
			}
		}
	}
	this->Mesh->numberofpoints_export = this->Mesh->numberofpoints;
	this->Mesh->pointlist_export = this->Mesh->pointlist;
}

void C_Model::ExportFeFlow(){
	//initialize FeFlow without having a individual constructur
	C_FeFlow * FeFlow = new C_FeFlow;
	int havewritten = 0;
	int minMat, maxMat;
	this->tranformBackward();
	QFile file(FileNameTmp);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;
	QTextStream out(&file);
	out.setRealNumberPrecision(24);
	out << "PROBLEM:\n";
	out << "CLASS (v.7)\n";
	out << "   2    1    0    3    0    0    8    8    0    0\n";
	out << "DIMENS\n";
	out << "   " << this->Mesh->numberofpoints << "     " << this->Mesh->numberoftetrahedra <<"     0      1      0      0      0      0      0      2     0      0      1      0      0      0      0\n";
	out << "SCALE\n";
	out << "   1.0, 1.0, 1.0, 1.0, 0.0, 0.0\n";
	out << "VARNODE\n";
	out << "   " << this->Mesh->numberoftetrahedra << "     4     4\n";
	for (long t = 0; t != this->Mesh->numberoftetrahedra; t++)
		out << "   6     " << this->Mesh->tetrahedronlist[t * 4 + 0]+1 << "     " << this->Mesh->tetrahedronlist[t * 4 + 1]+1 << "     " << this->Mesh->tetrahedronlist[t * 4 + 2]+1 << "     " << this->Mesh->tetrahedronlist[t * 4 + 3]+1 << "\n";
	out << "XYZCOOR\n";
	for (long n = 0; n != this->Mesh->numberofpoints; n++)
		out << "     " << this->Mesh->pointlist[n * 3 + 0] << ", " << this->Mesh->pointlist[n * 3 + 1] << ", " << this->Mesh->pointlist[n * 3 + 2] << "\n";
	if (this->Mesh->numberoftetrahedra>0)
	{
		minMat = this->Mesh->tetrahedronmarkerlist[0];
		maxMat = this->Mesh->tetrahedronmarkerlist[0];
		for (long t = 0; t != this->Mesh->numberoftetrahedra; t++)
		{
			if (this->Mesh->tetrahedronmarkerlist[t] < minMat) minMat = this->Mesh->tetrahedronmarkerlist[t];
			if (this->Mesh->tetrahedronmarkerlist[t] > maxMat) maxMat = this->Mesh->tetrahedronmarkerlist[t];
		}
		out << "ELEMENTALSETS\n";
		for (int m = minMat; m != maxMat+1; m++)
		{
			out << "     \"Region: Name: R" << m << "\"";
			havewritten = 0;
			for (long t = 0; t != this->Mesh->numberoftetrahedra; t++)
			{
				if (this->Mesh->tetrahedronmarkerlist[t] == m)
				{
					if ((havewritten % 10) == 0)
					{
						out << "\n";
						out << "\t\t";
					}
					out << t + 1 << " ";
					havewritten++;
				}
			}
			out << "\n";
		}
	}
	//based on the tetrahedronlist all possible triangles are generated (4 for each tet). Common triangles will deleted and each triangle just exists once.
	FeFlow->generateAllTriangles(this->Mesh->numberoftetrahedra, this->Mesh->tetrahedronlist);
	if (this->Mesh->numberoftriangles > 0)
	{
		minMat = this->Mesh->trianglemarkerlist[0];
		maxMat = this->Mesh->trianglemarkerlist[0];
		for (long t = 0; t != this->Mesh->numberoftriangles; t++)
		{
			if (this->Mesh->trianglemarkerlist[t] < minMat) minMat = this->Mesh->trianglemarkerlist[t];
			if (this->Mesh->trianglemarkerlist[t] > maxMat) maxMat = this->Mesh->trianglemarkerlist[t];
		}
		out << "FACESETS\n";
		for (int m = minMat; m != maxMat+1; m++)
		{
			//all tetgen triangles with the marker m are stored in a list having not the right index
			FeFlow->generateUndefinedTriangles(m, this->Mesh->numberoftriangles, this->Mesh->trianglemarkerlist, this->Mesh->trianglelist);
			//all tetgen triangles with the marker m will get a feflow index
			FeFlow->generateDefinedTriangles();
			out << "     \"Surface: Name: S" << m << "\"";
			havewritten = 0;
			for (long t = 0; t != FeFlow->definedTriangles.length(); t++){
				if ((havewritten % 10) == 0){
					out << "\n";
					out << "\t\t";
				}
				out << FeFlow->definedTriangles[t].index+1 << " ";
				havewritten++;
			}
			out << "\n";
		}
	}

	//based on the tetrahedronlist all possible edges are generated (6 for each tet). Common edges will deleted and each edge just exists once.  
	FeFlow->generateAllEdges(this->Mesh->numberoftetrahedra, this->Mesh->tetrahedronlist);

	if (this->Mesh->numberofedges > 0){
		minMat = this->Mesh->edgemarkerlist[0];
		maxMat = this->Mesh->edgemarkerlist[0];
		for (long e = 0; e != this->Mesh->numberofedges; e++){
			if (this->Mesh->edgemarkerlist[e] < minMat) minMat = this->Mesh->edgemarkerlist[e];
			if (this->Mesh->edgemarkerlist[e] > maxMat) maxMat = this->Mesh->edgemarkerlist[e];
		}
		out << "EDGESETS\n";
		for (int m = minMat; m != maxMat + 1; m++){
			//all tetgen edges with the marker m are stored in a list having not the right index
			FeFlow->generateUndefinedEdges(m, this->Mesh->numberofedges, this->Mesh->edgemarkerlist, this->Mesh->edgelist);
			//all tetgen edges with the marker m will get a feflow index
			FeFlow->generateDefinedEdges();
			out << "     \"Polyline: Name: P" << m << "\"";
			havewritten = 0;
			for (long e = 0; e != FeFlow->definedEdges.length(); e++){
				if ((havewritten % 10) == 0){
					out << "\n";
					out << "\t\t";
				}
				out << FeFlow->definedEdges[e].index + 1 << " ";
				havewritten++;
			}
			out << "\n";
		}
	}

	out << "END\n";
	file.close();
	this->tranformForward();

	//feflow is deleted and all lists will be cleared.
	delete FeFlow;
}

void C_Model::ExportOGS(){
	long elementnumber = 0;
	this->tranformBackward();
	this->ExportLists_make();

	QFile file(FileNameTmp);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;
	QTextStream out(&file);
	out.setRealNumberPrecision(24);
	out << "#FEM_MSH" << "\n";
	out << " $PCS_TYPE" << "\n";
	out << "  NO_PCS" << "\n";
	out << " $NODES" << "\n";
	out << "  " << this->Mesh->numberofpoints_export << "\n";
	for (int n = 0; n != this->Mesh->numberofpoints_export; n++){
		out << "  " << n << "  " << this->Mesh->pointlist_export[n * 3 + 0] << " " << this->Mesh->pointlist_export[n * 3 + 1] << " " << this->Mesh->pointlist_export[n * 3 + 2] << "\n";
	}
	out << " $ELEMENTS" << "\n";
	out << "  " << this->Mesh->numberoftetrahedra_export + this->Mesh->numberoftriangles_export + this->Mesh->numberofedges_export << "\n";
	for (int t = 0; t != this->Mesh->numberoftetrahedra_export; t++){
		out << "  " << elementnumber++ << " " << this->Mesh->tetrahedronmarkerlist_export[t] << " tet " << this->Mesh->tetrahedronlist_export[t * 4 + 0] << " " << this->Mesh->tetrahedronlist_export[t * 4 + 1] << " " << this->Mesh->tetrahedronlist_export[t * 4 + 2] << " " << this->Mesh->tetrahedronlist_export[t * 4 + 3] << "\n";
	}
	for (int f = 0; f != this->Mesh->numberoftriangles_export; f++){
		out << "  " << elementnumber++ << " " << this->Mesh->trianglemarkerlist_export[f] << " tri " << this->Mesh->trianglelist_export[f * 3 + 0] << " " << this->Mesh->trianglelist_export[f * 3 + 1] << " " << this->Mesh->trianglelist_export[f * 3 + 2] << " " << "\n";
	}
	for (int e = 0; e != this->Mesh->numberofedges_export; e++){
		out << "  " << elementnumber++ << " " << this->Mesh->edgemarkerlist_export[e] << " line " << this->Mesh->edgelist_export[e * 2 + 0] << " " << this->Mesh->edgelist_export[e * 2 + 1] << " " << "\n";
	}
	out << "#STOP" << "\n";
	file.close();
	this->ExportLists_clean();
	this->tranformForward();
}

void C_Model::ExportABAQUS(QString borderIDs)
{
	int currentMaterialID = -1;
	long elementCount = 0;
	this->tranformBackward();
	this->ExportLists_make();
	QFile file(FileNameTmp);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;
	QTextStream out(&file);
	out.setRealNumberPrecision(24);
	out << "**\n**ABAQUSInputDeckGeneratedbyMeshIT\n**GeneratedusingMeshIT\n";
	out << "**\n**Template:ABAQUS/STANDARD3D\n**\n";
	// write NODE
	out << "*NODE\n";
	for (int n = 0; n != this->Mesh->numberofpoints_export; n++)
		out << n + 1 << "," << this->Mesh->pointlist_export[n * 3 + 0] << "," << this->Mesh->pointlist_export[n * 3 + 1] << "," << this->Mesh->pointlist_export[n * 3 + 2] << "\n";
	// write ELEMENT
	// 3D tetrahedrons --> C3D4 - TET4
	
	for (int t = 0; t != this->Mesh->numberoftetrahedra_export; t++)
	{
		if (currentMaterialID != this->Mesh->tetrahedronmarkerlist_export[t]){
			currentMaterialID = this->Mesh->tetrahedronmarkerlist_export[t];
			out << "*ELEMENT,TYPE=C3D4,ELSET=UNIT_" << currentMaterialID << "\n";
		}
		out << ++elementCount << "," << this->Mesh->tetrahedronlist_export[t * 4 + 0] + 1 << "," << this->Mesh->tetrahedronlist_export[t * 4 + 1] + 1 << "," << this->Mesh->tetrahedronlist_export[t * 4 + 2] + 1 << "," << this->Mesh->tetrahedronlist_export[t * 4 + 3] + 1 << "\n";
	}
	// 2D triangles --> S3R - S3 - CPS3
	// but first organize the lists
	for (int f = 0; f != this->Mesh->numberoftriangles_export; f++)
	{
		if (currentMaterialID != this->Mesh->trianglemarkerlist_export[f]){
			currentMaterialID = this->Mesh->trianglemarkerlist_export[f];
			out << "*ELEMENT,TYPE=S3,ELSET=FAULT_" << currentMaterialID << "\n";
		}
		out << ++elementCount << "," << this->Mesh->trianglelist_export[f * 3 + 0] + 1 << "," << this->Mesh->trianglelist_export[f * 3 + 1] + 1 << "," << this->Mesh->trianglelist_export[f * 3 + 2] + 1 << "\n";
	}
	// 1D edges --> S3R - S3 - CPS3
	// but first organize the lists
	for (int e = 0; e != this->Mesh->numberofedges_export; e++)
	{
		if (currentMaterialID != this->Mesh->edgemarkerlist_export[e]){
			currentMaterialID = this->Mesh->edgemarkerlist_export[e];
			out << "*ELEMENT,TYPE=T3D2,ELSET=WELL_" << currentMaterialID << "\n";
		}
		out << ++elementCount << "," << this->Mesh->edgelist_export[e * 2 + 0] + 1 << "," << this->Mesh->edgelist_export[e * 2 + 1] + 1 << "\n";
	}
	int pos = 0;
	int marker;
	while (borderIDs.section(",", pos, pos).length()>0){
		marker = borderIDs.section(",", pos, pos).toInt();
		out << "*ELEMENT,TYPE=S3,ELSET=BORDER_" << marker << "\n";
		for (int f = 0; f != this->Mesh->trianglemarkerlist.length(); f++){
			if (this->Mesh->trianglemarkerlist[f] == marker){
				out << ++elementCount << "," << this->Mesh->trianglelist[f * 3 + 0] + 1 << "," << this->Mesh->trianglelist[f * 3 + 1] + 1 << "," << this->Mesh->trianglelist[f * 3 + 2] + 1 << "\n";
			}
		}
		pos++;
	}
	out << "*****";
	file.close();
	this->ExportLists_clean();
	this->tranformForward();
}

void C_Model::ExportVTU3D(){
	long offset=0;

	this->tranformBackward();
	this->ExportLists_make();

	this->VTU.clear();  
	this->VTU.NumberOfPoints=this->Mesh->numberofpoints_export;  
	this->VTU.NumberOfCells = this->Mesh->numberoftetrahedra_export + this->Mesh->numberoftriangles_export + this->Mesh->numberofedges_export;

	
	for (int p = 0; p != this->Mesh->numberofpoints_export; p++){
		this->VTU.Points.append(C_Vector3D(this->Mesh->pointlist_export[p * 3 + 0], this->Mesh->pointlist_export[p * 3 + 1], this->Mesh->pointlist_export[p * 3 + 2]));
		this->VTU.pointData.append(0);
	}
	for (int t = 0; t != this->Mesh->numberoftetrahedra_export; t++){
		this->VTU.connectivity.append(this->Mesh->tetrahedronlist_export[t * 4 + 0]);
		this->VTU.connectivity.append(this->Mesh->tetrahedronlist_export[t * 4 + 1]);
		this->VTU.connectivity.append(this->Mesh->tetrahedronlist_export[t * 4 + 2]);
		this->VTU.connectivity.append(this->Mesh->tetrahedronlist_export[t * 4 + 3]);
		offset+=4;
		this->VTU.offsets.append(offset);
		this->VTU.types.append(10);  
		this->VTU.matType.append(this->Mesh->tetrahedronmarkerlist_export[t]);
	}
	for (int f = 0; f != this->Mesh->numberoftriangles_export; f++){
		this->VTU.connectivity.append(this->Mesh->trianglelist_export[f * 3 + 0]);
		this->VTU.connectivity.append(this->Mesh->trianglelist_export[f * 3 + 1]);
		this->VTU.connectivity.append(this->Mesh->trianglelist_export[f * 3 + 2]);
		offset+=3;
		this->VTU.offsets.append(offset);  
		this->VTU.types.append(5);  
		this->VTU.matType.append(this->Mesh->trianglemarkerlist_export[f]);
	}
	for (int e = 0; e != this->Mesh->numberofedges_export; e++){
		this->VTU.connectivity.append(this->Mesh->edgelist_export[e * 2 + 0]);
		this->VTU.connectivity.append(this->Mesh->edgelist_export[e * 2 + 1]);
		offset+=2;
		this->VTU.offsets.append(offset);  
		this->VTU.types.append(3);  
		this->VTU.matType.append(this->Mesh->edgemarkerlist_export[e]);
	}
	this->VTU.write(this->FileNameTmp);
	this->ExportLists_clean();
	this->tranformForward();
}

#ifndef NOEXODUS

void C_Model::ExportEXODUS(QString borderIDs)
{
    int CPU_word_size = 0, IO_word_size = 0, exoid = 0, error = 0;
	tranformBackward();
	calculateNumberWithMaterials();
	int number = Mats.length();
	for (int s = 0; s != Surfaces.length(); s++)
	{
		if (Surfaces[s].MaterialID >= 0)
			number++;
	}
	for (int p = 0; p != Polylines.length(); p++)
	{
		if (Polylines[p].MaterialID >= 0)
			number++;
	}
	C_Exodus ex(3, Mesh->numberofpoints, Mesh->numberoftetrahedrawithmaterial + Mesh->numberoftriangleswithmaterial + Mesh->numberofedgeswithmaterial, Surfaces.length(), number, Mesh->numberoftetrahedra);
	EXODUS_sides(&ex, borderIDs);
	EXODUS_element(&ex);
	EXODUS_nodes(&ex);
    exoid = ex_create(FileNameTmp.toLatin1().data(), EX_CLOBBER | EX_LARGE_MODEL, &CPU_word_size, &IO_word_size);
    error = ex_put_init(exoid, "MeshIt export", ex.num_dim, ex.num_nodes, ex.num_elem, ex.num_elem_blk, ex.num_node_sets, ex.num_side_sets);
	error = ex_put_coord(exoid, ex.x, ex.y, ex.z);
	error = ex_put_coord_names(exoid, const_cast<char**>(ex.coord_names));
	error = ex_put_node_num_map(exoid, ex.node_num_map);
	int *connect;
	int count = 0;
	QString elem_type;
	for (int b = 0; b != ex.num_elem_blk; b++)
	{
		if (ex.num_nodes_per_elem[b] == 4)
			elem_type = "TETRA";
		else if (ex.num_nodes_per_elem[b] == 3)
			elem_type = "TRIANGLE";
		else if (ex.num_nodes_per_elem[b] == 2)
			elem_type = "TRUSS"; // BEAM
		error = ex_put_elem_block(exoid, ex.ebids[b], elem_type.toLatin1().data(), ex.num_elem_in_blk[b], ex.num_nodes_per_elem[b], 0);
		connect = new int[ex.num_elem_in_blk[b] * ex.num_nodes_per_elem[b]];
		for (int c = 0; c != ex.connections[b].length(); c++)
			connect[c] = ex.connections[b][c];
		error = ex_put_elem_conn(exoid, ex.ebids[b], connect);
		delete[] connect;
	}
	for (int e = 0; e != ex.num_elem; e++)
		ex.elem_num_map[e] = e + 1;
	error = ex_put_map(exoid, ex.elem_num_map);
	int pos = 0;
	int marker = 0;
	int nborder = 0;
	int num_side = 0;
	int tri_pos = 0;
	int mem_tri = 0;
	ex.generateAllTriangles(Mesh->numberoftetrahedra, Mesh->tetrahedronlist);
	while (borderIDs.section(",", pos, pos).length() > 0)
	{
		marker = borderIDs.section(",", pos, pos).toInt();
		for (int s = 0; s != Surfaces.length(); s++)
		{
			if (s == marker)
			{
				++nborder;
				num_side = ex.generateMarkedTriangles(marker, Mesh->numberoftriangles, Mesh->trianglemarkerlist, Mesh->trianglelist);
				ex.elem_list = new int[num_side];
				ex.side_list = new int[num_side];
				mem_tri = 0;
				for (int n = 0; n != num_side; n++)
				{
					tri_pos = ex.findTriangle(ex.markedTriangles[n], mem_tri);
					ex.side_list[n] = ex.allTriangles[tri_pos].index_side;
					ex.elem_list[n] = ex.allTriangles[tri_pos].index_blk;
					mem_tri = tri_pos;
				}
				std::cout << "side [ " << pos + 1 << " out of " << ex.num_side_sets << " ] ... done\n";
				error = ex_put_side_set_param(exoid, nborder, num_side, 0);
				error = ex_put_side_set(exoid, nborder, ex.elem_list, ex.side_list);
				delete[] ex.elem_list;
				delete[] ex.side_list;
				break;
			}
		}
		pos++;
	}
	ex.deallocate();
	error = ex_close(exoid);
    this->tranformForward();
}

void C_Model::EXODUS_sides(C_Exodus *exo, QString borderIDs)
{
	int n_side = 0;
	int marker;
	int pos = 0;
	for (int s = 0; s != exo->num_surfaces; s++)
		exo->num_side_ss[s] = 0;
	while (borderIDs.section(",", pos, pos).length()>0)
	{
		marker = borderIDs.section(",", pos, pos).toInt();
		for (int s = 0; s != exo->num_surfaces; s++)
			if (s == marker) n_side++;
		pos++;
	}
	exo->num_side_sets = n_side;
	exo->num_node_sets = 0;
}

void C_Model::EXODUS_element(C_Exodus *exo)
{
	bool found = false;
	QString elem_type;
	int pos = 0;
	int blk_id_length = 0;
	int elem_count = 1;
	int mat = 0;
	for (int t = 0; t != Mesh->numberoftetrahedra; t++)
	{
		found = true;
		mat = getMaterial(10, t);
		for (int b = 0; b != blk_id_length; b++)
		{
			if (mat == exo->ebids[b])
			{
				exo->num_elem_in_blk[b]++;
				exo->connections[b].append((Mesh->tetrahedronlist[t * 4 + 0]) + 1);
				exo->connections[b].append((Mesh->tetrahedronlist[t * 4 + 1]) + 1);
				exo->connections[b].append((Mesh->tetrahedronlist[t * 4 + 2]) + 1);
				exo->connections[b].append((Mesh->tetrahedronlist[t * 4 + 3]) + 1);
				found = false;
				break;
			}
		}
		if (found)
		{
			if (mat != -1)
			{
				exo->ebids[blk_id_length] = mat;
				exo->num_nodes_per_elem[pos++] = 4;
				exo->num_elem_in_blk[blk_id_length]++;
				exo->connections[blk_id_length].append((Mesh->tetrahedronlist[t * 4 + 0]) + 1);
				exo->connections[blk_id_length].append((Mesh->tetrahedronlist[t * 4 + 1]) + 1);
				exo->connections[blk_id_length].append((Mesh->tetrahedronlist[t * 4 + 2]) + 1);
				exo->connections[blk_id_length].append((Mesh->tetrahedronlist[t * 4 + 3]) + 1);
				blk_id_length++;
			}
		}
	}
	found = false;
	for (int t = 0; t != Mesh->numberoftriangles; t++)
	{
		found = true;
		mat = getMaterial(5, t);
		for (int b = 0; b != blk_id_length; b++)
		{
			if (mat == exo->ebids[b])
			{
				exo->num_elem_in_blk[b]++;
				exo->connections[b].append((Mesh->trianglelist[t * 3 + 0]) + 1);
				exo->connections[b].append((Mesh->trianglelist[t * 3 + 1]) + 1);
				exo->connections[b].append((Mesh->trianglelist[t * 3 + 2]) + 1);
				found = false;
				break;
			}
		}
		if (found)
		{
			if (mat != -1)
			{
				exo->ebids[blk_id_length] = mat;
				exo->num_nodes_per_elem[pos++] = 3;
				exo->num_elem_in_blk[blk_id_length]++;
				exo->connections[blk_id_length].append((Mesh->trianglelist[t * 3 + 0]) + 1);
				exo->connections[blk_id_length].append((Mesh->trianglelist[t * 3 + 1]) + 1);
				exo->connections[blk_id_length].append((Mesh->trianglelist[t * 3 + 2]) + 1);
				blk_id_length++;
			}
		}
	}
	found = false;
	for (int e = 0; e != Mesh->numberofedges; e++)
	{
		found = true;
		mat = getMaterial(3, e);
		for (int b = 0; b != blk_id_length; b++)
		{
			if (mat == exo->ebids[b])
			{
				exo->num_elem_in_blk[b]++;
				exo->connections[b].append(Mesh->edgelist[e * 2 + 0] + 1);
				exo->connections[b].append(Mesh->edgelist[e * 2 + 1] + 1);
				found = false;
				break;
			}
		}
		if (found)
		{
			if (mat != -1)
			{
				exo->ebids[blk_id_length] = mat;
				exo->num_nodes_per_elem[pos++] = 2;
				exo->num_elem_in_blk[blk_id_length]++;
				exo->connections[blk_id_length].append(Mesh->edgelist[e * 2 + 0] + 1);
				exo->connections[blk_id_length].append(Mesh->edgelist[e * 2 + 1] + 1);
				blk_id_length++;
			}
		}
	}
}

void C_Model::EXODUS_nodes(C_Exodus *exo)
{
	double x_trans = (this->min.x() + this->max.x()) * 0.5;
	double y_trans = (this->min.y() + this->max.y()) * 0.5;
	double angle = double(this->ExportRotationAngle) / double(180) * MY_PI;
	double tmp_x;
	double tmp_y;
	for (int n = 0; n != exo->num_nodes; n++)
	{
		exo->x[n] = Mesh->pointlist[n * 3 + 0] - x_trans;
		exo->y[n] = Mesh->pointlist[n * 3 + 1] - y_trans;
		exo->z[n] = Mesh->pointlist[n * 3 + 2];
		tmp_x = exo->x[n] * cos(angle) - exo->y[n] * sin(angle);
		tmp_y = exo->x[n] * sin(angle) + exo->y[n] * cos(angle);
		exo->x[n] = tmp_x;
		exo->y[n] = tmp_y;
		exo->node_num_map[n] = n + 1;
	}
}

/*void C_Model::EXODUS_nodes(C_Exodus *exo)
{
	for (int n = 0; n != exo->num_nodes; n++)
	{
		exo->x[n] = Mesh->pointlist[n * 3 + 0];
		exo->y[n] = Mesh->pointlist[n * 3 + 1];
		exo->z[n] = Mesh->pointlist[n * 3 + 2];
		exo->node_num_map[n] = n + 1;
	}
}*/


#endif

void C_Model::ExportCOMSOL()
{
	this->tranformBackward();
    QFile file(FileNameTmp);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;
	QTextStream out(&file);
	out.setRealNumberPrecision(24); 
	out<<"# Created by MeshIT\n";
	out << "\n";
	out << "# Major & minor version\n";
	out<<"0 1\n";
	out << "\n";
	out<<"1 # number of tags\n";
	out<<"# Tags\n";
	out<<"5 MESH\n";
	out << "\n";
	out<<"1 # number of types\n";
	out<<"# Types\n";
	out<<"3 obj\n";
	out << "\n";
	out<<"# A 3D mesh object\n";
	out<<"0 0 1\n";
	out<<"4 Mesh # class\n";
	out<<"4 # version\n";
	out<<"3 # sdim\n";
	out << this->Mesh->numberofpoints << " # number of mesh points\n";
	out<<"0 # lowest mesh point index\n";
	out << "\n";

	out<<"# Mesh point coordinates\n";
	for (int n = 0; n != this->Mesh->numberofpoints; n++){
		out << this->Mesh->pointlist[n * 3 + 0] << " " << this->Mesh->pointlist[n * 3 + 1] << " " << this->Mesh->pointlist[n * 3 + 2] << "\n";
	}
	out<<"\n";

	out<<"4 # number of element types\n";
	out << "\n";


	out << "# Type #0\n";
	out << "3 vtx # type name\n";
	out << "1 # number of nodes per element\n";
	out << "0 # number of elements\n";
	out << "# Elements\n";
	out << "#...vertex IDs\n";
	out << "\n";

	out << "0 # number of geometric entity indices\n";
	out << "# Geometric entity indices\n";
	out << "#...domain markers\n";
	out << "\n";

	out << "# Type #1\n";
	out << "3 edg # type name\n";
	out << "2 # number of nodes per element\n";
	out << this->Mesh->numberofedges << " # number of elements\n";
	out << "# Elements\n";
	for (int f = 0; f != this->Mesh->numberofedges; f++){
		out << this->Mesh->edgelist[f * 2 + 0] << " " << this->Mesh->edgelist[f * 2 + 1] << " " << "\n";
	}
	out << "\n";

	out << this->Mesh->numberofedges << " # number of geometric entity indices\n";
	out << "# Geometric entity indices\n";
	for (int f = 0; f != this->Mesh->numberofedges; f++){
		out << this->Mesh->edgemarkerlist[f] << " " << "\n";
	}
	out << "\n";

	out << "# Type #3\n";
	out << "3 tri # type name\n";
	out << "3 # number of nodes per element\n";
	out << this->Mesh->numberoftriangles << " # number of elements\n";
	out<<"# Elements\n";
	for (int f = 0; f != this->Mesh->numberoftriangles; f++){
		out << this->Mesh->trianglelist[f * 3 + 0] << " " << this->Mesh->trianglelist[f * 3 + 1] << " " << this->Mesh->trianglelist[f * 3 + 2] << " " << "\n";
	}
	out<<"\n";

	out << this->Mesh->numberoftriangles << " # number of geometric entity indices\n";
	out << "# Geometric entity indices\n";
	for (int f = 0; f != this->Mesh->numberoftriangles; f++){
		out << this->Mesh->trianglemarkerlist[f] << " " << "\n";
	}
	out << "\n";

	out << "# Type #4\n";
	out << "3 tet # type name\n";
	out << "4 # number of nodes per element\n";
	out << this->Mesh->numberoftetrahedra << " # number of elements\n";
	out << "# Elements\n";
	for (int f = 0; f != this->Mesh->numberoftetrahedra; f++){
		out << this->Mesh->tetrahedronlist[f * 4 + 0] << " " << this->Mesh->tetrahedronlist[f * 4 + 1] << " " << this->Mesh->tetrahedronlist[f * 4 + 2] << " " << this->Mesh->tetrahedronlist[f * 4 + 3] << " " << "\n";
	}
	out << "\n";

	out << this->Mesh->numberoftetrahedra << " # number of geometric entity indices\n";
	out << "# Geometric entity indices\n";
	for (int f = 0; f != this->Mesh->numberoftetrahedra; f++){
		out << this->Mesh->tetrahedronmarkerlist[f]+1 << " " << "\n";
	}
	out << "\n";

	file.close();
	this->tranformForward(); 
}

void C_Model::ExportVTU2D(QString surfaceIDs){
	double point[3];
	int offset=0;
	int connection=0;
	int pos=0;
	int marker;
	this->tranformBackward();
	this->VTU.clear();  
	while (surfaceIDs.section(",",pos,pos).length()>0){ 
		marker=surfaceIDs.section(",",pos,pos).toInt();
		for (int f=0;f!=this->Mesh->trianglemarkerlist.length();f++){
			if (this->Mesh->trianglemarkerlist[f]==marker){
				this->VTU.NumberOfCells++;
				for (int p=0;p!=3;p++){
					this->Mesh->getCoordinates(this->Mesh->trianglelist[f*3+p],point);   
					this->VTU.Points.append(C_Vector3D(point[0],point[1],point[2]));
					this->VTU.connectivity.append(connection*3+p);
				}
				connection++;
				this->VTU.offsets.append(offset*3+3);
				offset++;
				this->VTU.types.append(5);  
			}			
		}
		pos++;
	}
	this->VTU.NumberOfPoints=this->VTU.Points.length();
	this->VTU.write(this->FileNameTmp);
	this->tranformForward();
}


void C_Model::ExportTIN(QString surfaceIDs){
	double point[3][3];
	int faceNumber=0;
	int pos=0;
	int marker;
	this->tranformBackward(); 
    QFile file(FileNameTmp);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;
	QTextStream out(&file);
	out.setRealNumberPrecision(24);
	while (surfaceIDs.section(",",pos,pos).length()>0){ 
		marker=surfaceIDs.section(",",pos,pos).toInt();
		for (int f=0;f!=this->Mesh->trianglemarkerlist.length();f++){
			if (this->Mesh->trianglemarkerlist[f]==marker){
				for (int p=0;p!=3;p++) this->Mesh->getCoordinates(this->Mesh->trianglelist[f*3+p],point[p]);   
				out<<faceNumber++<<" ";
				out<<point[0][0]<<" "<<point[0][1]<<" "<<point[0][2]<<" "; 
				out<<point[1][0]<<" "<<point[1][1]<<" "<<point[1][2]<<" "; 
				out<<point[2][0]<<" "<<point[2][1]<<" "<<point[2][2]<<"\n"; 
			}
		}
		pos++;
	}
	file.close();
	this->tranformForward(); 
}

QString extract(QString input, QString request){
	QString handle = input.section(request, -1, -1);
	return handle.section("\"", 1, 1);
}

void C_Model::Open(){
	QString type, name, part, filename;
	double size;
	int mat;
	QString line,tmp;
	QFile file(FileNameModel);
	QFileInfo fileNameModel_info(FileNameModel);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
	QTextStream in(&file);
	while (!in.atEnd()){
		line=in.readLine().simplified();
		if (line.contains("Props")){
			int pos;
			bool ok;
			tmp = line.right(line.length()-line.indexOf("method"));
			this->intAlgorythm=tmp.section("\"",1,1);

			/* Read premesh gradient if available. */
			pos = line.indexOf("gradient_premesh");
			if( pos >= 0 ) {
				QString premesh_gradient = line.right(line.length() - pos).section("\"", 1, 1);
				double radius = premesh_gradient.toDouble(&ok);
				if( ok )
					preMeshGradient = radius;
			}

			/* Read mesh gradient if available. */
			pos = line.indexOf("gradient_mesh");
			if( pos >= 0 ) {
				QString mesh_gradient = line.right(line.length() - pos).section("\"", 1, 1);
				double radius = mesh_gradient.toDouble(&ok);
				if( ok )
					meshGradient = radius;
			}
		}
		if (line.contains("DataSet")){
			type = extract(line, " group=");
			part = extract(line, " part=");
			filename = extract(line, " file=");
			name = extract(line, " name=");
			mat = extract(line, " material=").toInt();
			size = extract(line, " size=").toDouble();


			if (part=="SD"){
				this->VTU.read(fileNameModel_info.absoluteDir().absolutePath()+"/"+filename);
				if (type=="WELL"){
					this->VTU_SD_to_Polyline(name, type, mat, size); 
				}else{
					this->VTU_SD_to_Surface(name, type, mat, size); 
				}
			}
			if (part=="SEG"){
				this->VTU.read(fileNameModel_info.absoluteDir().absolutePath()+"/"+filename);
				this->VTU_SEG_to_Polyline(name, type); 
			}
			if (part=="CH"){
				this->VTU.read(fileNameModel_info.absoluteDir().absolutePath()+"/"+filename);
				this->VTU_CH_to_Surface(name, type); 
			}
			if (part=="TRI"){
				this->VTU.read(fileNameModel_info.absoluteDir().absolutePath()+"/"+filename);
				this->VTU_TRI_to_Surface(name, type); 
			}
			if (part=="CON"){
				this->VTU.read(fileNameModel_info.absoluteDir().absolutePath()+"/"+filename);
				if (type=="WELL"){
					this->VTU_CON_to_Polyline(name, type); 
				}else{
					this->VTU_CON_to_Surface(name, type); 
				}
			}
			if (type=="INTERSECTIONS"){
				this->VTU.read(fileNameModel_info.absoluteDir().absolutePath()+"/"+filename);
				this->VTU_INT_to_Model();
			}
			if (type=="MATS"){
				this->VTU.read(fileNameModel_info.absoluteDir().absolutePath()+"/"+filename);
				this->VTU_MAT_to_Model(); 
			}
			if (type=="TETS"){
				this->VTU.read(fileNameModel_info.absoluteDir().absolutePath()+"/"+filename);
				this->VTU_TET_to_Model(); 
			}
		}
	}
	this->calculate_min_max();
	this->tranformForward();
	this->calculate_size_of_intersections();
	this->calculate_size_of_constraints();
	file.close();
}

void C_Model::Save(){
	this->tranformBackward();
	QFileInfo fileNameModel_info(FileNameModel);
    QString ModelName = fileNameModel_info.completeBaseName();
    QString ModelPath = fileNameModel_info.absoluteDir().absolutePath() + "/" + fileNameModel_info.completeBaseName(); 
	QDir VTUDir;
	VTUDir.mkdir(ModelPath+"/"); 
    QFile file(this->FileNameModel);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;
	QTextStream out(&file);
	out.setRealNumberPrecision(24); 
	out<<"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"; 
	out<<"  <Model>\n";
	out << "    <Props method=\"" << this->intAlgorythm<< "\" "
	    << "gradient_premesh=\"" << preMeshGradient << "\" "
	    << "gradient_mesh=\"" << meshGradient << "\""
	    << "/>\n";
	out<<"  </Model>\n";
	out<<"  <Collection>\n";
	for (int s=0;s!=this->Surfaces.length();s++){
		if (this->Surfaces[s].SDs.length()!=0){ 
			out << "    <DataSet group=\"" + this->Surfaces[s].Type + "\" part=\"SD\" file=\"" + ModelName + "/" + Surfaces[s].Name + "_SD.vtu" + "\" name=\"" + this->Surfaces[s].Name + "\" material=\"" << this->Surfaces[s].MaterialID << "\" size=\"" << this->Surfaces[s].size << "\"/>\n";
			Surfaces[s].makeVTU_SD();
			Surfaces[s].VTU.write(ModelPath+"/"+Surfaces[s].Name+"_SD.vtu"); 
		}
		if (this->Surfaces[s].ConvexHull.Ns.length()!=0){ 
			out<<"    <DataSet group=\""+this->Surfaces[s].Type+"\" part=\"CH\" file=\""+ModelName+"/"+Surfaces[s].Name+"_CH.vtu"+"\" name=\""+this->Surfaces[s].Name+"\"/>\n";
			Surfaces[s].makeVTU_CH();
			Surfaces[s].VTU.write(ModelPath+"/"+Surfaces[s].Name+"_CH.vtu"); 
		}
		if (this->Surfaces[s].Ts.length()!=0){ 
			out<<"    <DataSet group=\""+this->Surfaces[s].Type+"\" part=\"TRI\" file=\""+ModelName+"/"+Surfaces[s].Name+"_TRI.vtu"+"\" name=\""+this->Surfaces[s].Name+"\"/>\n";
			Surfaces[s].makeVTU_TRI();
			Surfaces[s].VTU.write(ModelPath+"/"+Surfaces[s].Name+"_TRI.vtu"); 
		}
		if (this->Surfaces[s].Constraints.length()!=0){ 
			out<<"    <DataSet group=\""+this->Surfaces[s].Type+"\" part=\"CON\" file=\""+ModelName+"/"+Surfaces[s].Name+"_CON.vtu"+"\" name=\""+this->Surfaces[s].Name+"\"/>\n";
			Surfaces[s].makeVTU_CON();
			Surfaces[s].VTU.write(ModelPath+"/"+Surfaces[s].Name+"_CON.vtu"); 
		}
	}
	for (int p=0;p!=this->Polylines.length();p++){
		if (this->Polylines[p].SDs.length()!=0){ 
			out << "    <DataSet group=\"" + this->Polylines[p].Type + "\" part=\"SD\" file=\"" + ModelName + "/" + Polylines[p].Name + "_SD.vtu" + "\" name=\"" + this->Polylines[p].Name + "\" material=\"" << this->Polylines[p].MaterialID << "\" size=\"" << this->Polylines[p].size << "\"/>\n";
			Polylines[p].makeVTU_SD();
			Polylines[p].VTU.write(ModelPath+"/"+Polylines[p].Name+"_SD.vtu"); 
		}
		if (this->Polylines[p].Constraints.length()!=0){ 
			out<<"    <DataSet group=\""+this->Polylines[p].Type+"\" part=\"SEG\" file=\""+ModelName+"/"+Polylines[p].Name+"_SEG.vtu"+"\" name=\""+this->Polylines[p].Name+"\"/>\n";
			Polylines[p].makeVTU_SEG();
			Polylines[p].VTU.write(ModelPath+"/"+Polylines[p].Name+"_SEG.vtu"); 
		}
		if (this->Polylines[p].Constraints.length()!=0){ 
			out<<"    <DataSet group=\""+this->Polylines[p].Type+"\" part=\"CON\" file=\""+ModelName+"/"+Polylines[p].Name+"_CON.vtu"+"\" name=\""+this->Polylines[p].Name+"\"/>\n";
			Polylines[p].makeVTU_CON();
			Polylines[p].VTU.write(ModelPath+"/"+Polylines[p].Name+"_CON.vtu"); 
		}
	}
	if (this->Intersections.length()!=0){ 
		out<<"    <DataSet group=\"INTERSECTIONS\" part=\"\" file=\""+ModelName+"/"+ModelName+"_INT.vtu"+"\" name=\""+ModelName+"\"/>\n";
		this->makeVTU_INT();
		this->VTU.write(ModelPath+"/"+ModelName+"_INT.vtu"); 
	}
	if (this->Mats.length()!=0){ 
		out<<"    <DataSet group=\"MATS\" part=\"\" file=\""+ModelName+"/"+ModelName+"_MAT.vtu"+"\" name=\""+ModelName+"\"/>\n";
		this->makeVTU_MAT();
		this->VTU.write(ModelPath+"/"+ModelName+"_MAT.vtu"); 
	}
	if (this->Mesh){ 
		out<<"    <DataSet group=\"TETS\" part=\"\" file=\""+ModelName+"/"+ModelName+"_TET.vtu"+"\" name=\""+ModelName+"\"/>\n";
		this->makeVTU_TET();
		this->VTU.write(ModelPath+"/"+ModelName+"_TET.vtu"); 
	}
	out<<"  </Collection>\n"; 
	out<<"</VTKFile>\n"; 
	file.close();
	this->tranformForward(); 
}

void C_Model::tranformForward(){
	C_Vector3D diffMaxShift;
	double maxValue;
	this->shift=(this->min+this->max)/2;
	diffMaxShift=this->max-this->shift;
	maxValue=diffMaxShift.x();
	if (diffMaxShift.y()>maxValue) maxValue=diffMaxShift.y();
	if (diffMaxShift.z()>maxValue) maxValue=diffMaxShift.z();
	this->scale=1/maxValue;  
	for (int s=0;s!=this->Surfaces.length();s++){
		this->Surfaces[s].size *= scale;
		for (int sd=0;sd!=this->Surfaces[s].SDs.length();sd++){
			this->Surfaces[s].SDs[sd]-=shift;
			this->Surfaces[s].SDs[sd]*=this->scale;
		}
		for (int n=0;n!=this->Surfaces[s].Ns.length();n++){
			this->Surfaces[s].Ns[n]-=shift;
			this->Surfaces[s].Ns[n]*=this->scale;
		}
		for (int n=0;n!=this->Surfaces[s].ConvexHull.Ns.length();n++){
			this->Surfaces[s].ConvexHull.Ns[n]-=shift;
			this->Surfaces[s].ConvexHull.Ns[n]*=this->scale;
		}
		for (int c=0;c!=this->Surfaces[s].Constraints.length();c++){
			for (int n=0;n!=this->Surfaces[s].Constraints[c].Ns.length();n++){
				this->Surfaces[s].Constraints[c].Ns[n]-=shift;
				this->Surfaces[s].Constraints[c].Ns[n]*=this->scale;
			}
		}
	}
	for (int p=0;p!=Polylines.length();p++){
		this->Polylines[p].size *= scale;
		for (int sd = 0; sd != this->Polylines[p].SDs.length(); sd++){
			this->Polylines[p].SDs[sd]-=this->shift;
			this->Polylines[p].SDs[sd]*=this->scale;
		}
		for (int n=0;n!=this->Polylines[p].Path.Ns.length();n++){
			this->Polylines[p].Path.Ns[n]-=this->shift;
			this->Polylines[p].Path.Ns[n]*=this->scale;
		}
		for (int c=0;c!=this->Polylines[p].Constraints.length();c++){
			for (int n=0;n!=this->Polylines[p].Constraints[c].Ns.length();n++){
				this->Polylines[p].Constraints[c].Ns[n]-=this->shift;
				this->Polylines[p].Constraints[c].Ns[n]*=this->scale;
			}
		}
	}
	for (int i=0;i!=this->Intersections.length();i++){
		for (int n=0;n!=this->Intersections[i].Ns.length();n++){
			this->Intersections[i].Ns[n]-=shift;
			this->Intersections[i].Ns[n]*=this->scale;
		}
	}
	if (this->Mesh){
		for (int p=0;p!=this->Mesh->numberofpoints;p++){
			this->Mesh->pointlist[p*3+0]-=shift.x(); 
			this->Mesh->pointlist[p*3+0]*=scale; 
			this->Mesh->pointlist[p*3+1]-=shift.y(); 
			this->Mesh->pointlist[p*3+1]*=scale; 
			this->Mesh->pointlist[p*3+2]-=shift.z(); 
			this->Mesh->pointlist[p*3+2]*=scale; 
		}
	}
	for (int m=0;m!=this->Mats.length();m++){
		for (int l=0;l!=this->Mats[m].Locations.length();l++){
			this->Mats[m].Locations[l]-=shift;
			this->Mats[m].Locations[l]*=this->scale;
		}
	}
}

void C_Model::tranformBackward(){
	for (int s=0;s!=this->Surfaces.length();s++){
		this->Surfaces[s].size /= scale;
		for (int sd=0;sd!=this->Surfaces[s].SDs.length();sd++){
			this->Surfaces[s].SDs[sd]/=scale;
			this->Surfaces[s].SDs[sd]+=shift;
		}
		for (int n=0;n!=this->Surfaces[s].Ns.length();n++){
			this->Surfaces[s].Ns[n]/=scale;
			this->Surfaces[s].Ns[n]+=shift;
		}
		for (int n=0;n!=this->Surfaces[s].ConvexHull.Ns.length();n++){
			this->Surfaces[s].ConvexHull.Ns[n]/=scale;
			this->Surfaces[s].ConvexHull.Ns[n]+=shift;
		}
		for (int c=0;c!=this->Surfaces[s].Constraints.length();c++){
			for (int n=0;n!=this->Surfaces[s].Constraints[c].Ns.length();n++){
				this->Surfaces[s].Constraints[c].Ns[n]/=scale;
				this->Surfaces[s].Constraints[c].Ns[n]+=shift;
			}
		}
	}
	for (int p=0;p!=Polylines.length();p++){
		this->Polylines[p].size /= scale;
		for (int sd=0;sd!=this->Polylines[p].SDs.length();sd++){
			this->Polylines[p].SDs[sd]/=scale;
			this->Polylines[p].SDs[sd]+=shift;
		}
		for (int n=0;n!=this->Polylines[p].Path.Ns.length();n++){
			this->Polylines[p].Path.Ns[n]/=scale;
			this->Polylines[p].Path.Ns[n]+=shift;
		}
		for (int c=0;c!=this->Polylines[p].Constraints.length();c++){
			for (int n=0;n!=this->Polylines[p].Constraints[c].Ns.length();n++){
				this->Polylines[p].Constraints[c].Ns[n]/=scale;
				this->Polylines[p].Constraints[c].Ns[n]+=shift;
			}
		}
	}
	for (int i=0;i!=this->Intersections.length();i++){
		for (int n=0;n!=this->Intersections[i].Ns.length();n++){
			this->Intersections[i].Ns[n]/=scale;
			this->Intersections[i].Ns[n]+=shift;
		}
	}
	if (this->Mesh){
		for (int p=0;p!=this->Mesh->numberofpoints;p++){
			this->Mesh->pointlist[p*3+0]/=scale; 
			this->Mesh->pointlist[p*3+0]+=shift.x(); 
			this->Mesh->pointlist[p*3+1]/=scale; 
			this->Mesh->pointlist[p*3+1]+=shift.y(); 
			this->Mesh->pointlist[p*3+2]/=scale; 
			this->Mesh->pointlist[p*3+2]+=shift.z(); 
		}
	}
	for (int m=0;m!=this->Mats.length();m++){
		for (int l=0;l!=this->Mats[m].Locations.length();l++){
			this->Mats[m].Locations[l]/=this->scale;
			this->Mats[m].Locations[l]+=shift;
		}
	}
	this->shift=C_Vector3D(0.0,0.0,0.0);
	this->scale=1.0; 
}

void C_Model::calculate_min_max(){
	if (this->Surfaces.length()!=0){
		min.setX(Surfaces[0].SDs[0].x()); 
		min.setY(Surfaces[0].SDs[0].y()); 
		min.setZ(Surfaces[0].SDs[0].z()); 
	}else if (this->Polylines.length()!=0){
		min.setX(Polylines[0].SDs[0].x()); 
		min.setY(Polylines[0].SDs[0].y()); 
		min.setZ(Polylines[0].SDs[0].z()); 
	}
	max.setX(min.x()); 
	max.setY(min.y()); 
	max.setZ(min.z()); 
	for (int s=0;s!=this->Surfaces.length();s++){
		Surfaces[s].calculate_min_max();
		if (Surfaces[s].min.x()<min.x()) min.setX(Surfaces[s].min.x()); 
		if (Surfaces[s].min.y()<min.y()) min.setY(Surfaces[s].min.y());
		if (Surfaces[s].min.z()<min.z()) min.setZ(Surfaces[s].min.z());
		if (Surfaces[s].max.x()>max.x()) max.setX(Surfaces[s].max.x()); 
		if (Surfaces[s].max.y()>max.y()) max.setY(Surfaces[s].max.y());
		if (Surfaces[s].max.z()>max.z()) max.setZ(Surfaces[s].max.z());
	}
	for (int p=0;p!=this->Polylines.length();p++){
		Polylines[p].calculate_min_max();
		if (Polylines[p].min.x()<min.x()) min.setX(Polylines[p].min.x());
		if (Polylines[p].min.y()<min.y()) min.setY(Polylines[p].min.y());
		if (Polylines[p].min.z()<min.z()) min.setZ(Polylines[p].min.z());
		if (Polylines[p].max.x()>max.x()) max.setX(Polylines[p].max.x());
		if (Polylines[p].max.y()>max.y()) max.setY(Polylines[p].max.y());
		if (Polylines[p].max.z()>max.z()) max.setZ(Polylines[p].max.z());
	}
	emit ModelInfoChanged();
}

void C_Model::calculate_size_of_intersections(){
	// for intersection take smallest size of intersecting features like polyline-surface or surface-surface intersection
	for (int i = 0; i != this->Intersections.length(); i++){
		if (this->Intersections[i].Ns.length() == 1)
		{
			this->Intersections[i].size = this->Polylines[this->Intersections[i].Object[1]].size < this->Surfaces[this->Intersections[i].Object[0]].size ? this->Polylines[this->Intersections[i].Object[1]].size : this->Surfaces[this->Intersections[i].Object[0]].size;
		}
		else
		{
			this->Intersections[i].size = this->Surfaces[this->Intersections[i].Object[1]].size < this->Surfaces[this->Intersections[i].Object[0]].size ? this->Surfaces[this->Intersections[i].Object[1]].size : this->Surfaces[this->Intersections[i].Object[0]].size;
		}
	}
}

void C_Model::calculate_size_of_constraints(){
	// for intersection take smallest size of intersecting features like polyline-surface or surface-surface intersection
	for (int s1 = 0; s1 != this->Surfaces.length(); s1++){
		for (int c1 = 0; c1 != this->Surfaces[s1].Constraints.length(); c1++){
			for (int s2 = 0; s2 != this->Surfaces.length(); s2++){
				for (int c2 = 0; c2 != this->Surfaces[s2].Constraints.length(); c2++){
					if (this->Surfaces[s1].Constraints[c1].IsIdenticallyWith(this->Surfaces[s2].Constraints[c2])){
						this->Surfaces[s1].Constraints[c1].size = this->Surfaces[s1].size < this->Surfaces[s2].size ? this->Surfaces[s1].size : this->Surfaces[s2].size;
						this->Surfaces[s2].Constraints[c2].size = this->Surfaces[s1].size < this->Surfaces[s2].size ? this->Surfaces[s1].size : this->Surfaces[s2].size;
					}
				}
			}
			for (int p = 0; p != this->Polylines.length(); p++){
				for (int c2 = 0; c2 != this->Polylines[p].Constraints.length(); c2++){
					if (this->Surfaces[s1].Constraints[c1].IsIdenticallyWith(this->Polylines[p].Constraints[c2])){
						this->Surfaces[s1].Constraints[c1].size = this->Surfaces[s1].size < this->Polylines[p].size ? this->Surfaces[s1].size : this->Polylines[p].size;
						this->Polylines[p].Constraints[c2].size = this->Surfaces[s1].size < this->Polylines[p].size ? this->Surfaces[s1].size : this->Polylines[p].size;
					}
				}
			}
		}
	}
}

bool C_Model::PreTestIntersectionsPolylineSurface(const C_Polyline * P, const C_Surface * S) const {
	if (P->Path.max.x()<S->min.x()) return false;
	if (P->Path.max.y()<S->min.y()) return false;
	if (P->Path.max.z()<S->min.z()) return false;
	if (P->Path.min.x()>S->max.x()) return false;
	if (P->Path.min.y()>S->max.y()) return false;
	if (P->Path.min.z()>S->max.z()) return false;
	return true;
}

bool C_Model::PreTestIntersectionsSegmentSurface(const C_Vector3D * S1, const C_Vector3D * S2, const C_Surface * S) const {
	C_Vector3D minWithError=S->min-C_Vector3D(1e-12,1e-12,1e-12);
	C_Vector3D maxWithError=S->max+C_Vector3D(1e-12,1e-12,1e-12);
	if (S1->x()<minWithError.x() && S2->x()<minWithError.x()) return false;
	if (S1->x()>maxWithError.x() && S2->x()>maxWithError.x()) return false;
	if (S1->y()<minWithError.y() && S2->y()<minWithError.y()) return false;
	if (S1->y()>maxWithError.y() && S2->y()>maxWithError.y()) return false;
	if (S1->z()<minWithError.z() && S2->z()<minWithError.z()) return false;
	if (S1->z()>maxWithError.z() && S2->z()>maxWithError.z()) return false;
	return true;
}

bool C_Model::PreTestIntersectionsTrianglePolyline(const C_Triangle * T, const C_Polyline * P) const {
	if (T->max.x()<P->Path.min.x()) return false;
	if (T->max.y()<P->Path.min.y()) return false;
	if (T->max.z()<P->Path.min.z()) return false;
	if (T->min.x()>P->Path.max.x()) return false;
	if (T->min.y()>P->Path.max.y()) return false;
	if (T->min.z()>P->Path.max.z()) return false;
	return true;
}

bool C_Model::PreTestIntersectionsSegmentTriangle(const C_Vector3D * S1, const C_Vector3D * S2, const C_Triangle * T) const {
	C_Vector3D minWithError=T->min-C_Vector3D(1e-12,1e-12,1e-12);
	C_Vector3D maxWithError=T->max+C_Vector3D(1e-12,1e-12,1e-12);
	if (S1->x()<minWithError.x() && S2->x()<minWithError.x()) return false;
	if (S1->x()>maxWithError.x() && S2->x()>maxWithError.x()) return false;
	if (S1->y()<minWithError.y() && S2->y()<minWithError.y()) return false;
	if (S1->y()>maxWithError.y() && S2->y()>maxWithError.y()) return false;
	if (S1->z()<minWithError.z() && S2->z()<minWithError.z()) return false;
	if (S1->z()>maxWithError.z() && S2->z()>maxWithError.z()) return false;
	return true;
}

bool C_Box::tri_in_box(const C_Triangle * Tri) const {
	if ((Tri->Ns[0]->x() < this->min.x()) && (Tri->Ns[1]->x() < this->min.x()) && (Tri->Ns[2]->x() < this->min.x())) return false;
	if ((Tri->Ns[0]->x() > this->max.x()) && (Tri->Ns[1]->x() > this->max.x()) && (Tri->Ns[2]->x() > this->max.x())) return false;
	if ((Tri->Ns[0]->y() < this->min.y()) && (Tri->Ns[1]->y() < this->min.y()) && (Tri->Ns[2]->y() < this->min.y())) return false;
	if ((Tri->Ns[0]->y() > this->max.y()) && (Tri->Ns[1]->y() > this->max.y()) && (Tri->Ns[2]->y() > this->max.y())) return false;
	if ((Tri->Ns[0]->z() < this->min.z()) && (Tri->Ns[1]->z() < this->min.z()) && (Tri->Ns[2]->z() < this->min.z())) return false;
	if ((Tri->Ns[0]->z() > this->max.z()) && (Tri->Ns[1]->z() > this->max.z()) && (Tri->Ns[2]->z() > this->max.z())) return false;
	return true;
}

bool C_Box::seg_in_box(const C_Vector3D * V1, const C_Vector3D * V2) const {
	if ((V1->x() < this->min.x()) && (V2->x() < this->min.x())) return false;
	if ((V1->x() > this->max.x()) && (V2->x() > this->max.x())) return false;
	if ((V1->y() < this->min.y()) && (V2->y() < this->min.y())) return false;
	if ((V1->y() > this->max.y()) && (V2->y() > this->max.y())) return false;
	if ((V1->z() < this->min.z()) && (V2->z() < this->min.z())) return false;
	if ((V1->z() > this->max.z()) && (V2->z() > this->max.z())) return false;
	return true;
}


void C_Box::calculate_center(){
	this->center.setX(0.5*(this->min.x() + this->max.x()));
	this->center.setY(0.5*(this->min.y() + this->max.y()));
	this->center.setZ(0.5*(this->min.z() + this->max.z()));
}

bool C_Box::too_much_tri() const {
	if (this->T1s.length() > 48) return true;
	if (this->T2s.length() > 48) return true;
	return false;
}

bool C_Box::too_much_seg() const {
	if (this->N1s.length() > 48) return true;
	if (this->N2s.length() > 48) return true;
	return false;
}

void C_Box::generate_subboxes(){
	this->Box[0] = new C_Box;
	Box[0]->min.setX(this->min.x());
	Box[0]->min.setY(this->min.y());
	Box[0]->min.setZ(this->min.z());
	Box[0]->max.setX(this->center.x());
	Box[0]->max.setY(this->center.y());
	Box[0]->max.setZ(this->center.z());
	Box[0]->calculate_center();
	this->Box[1] = new C_Box;
	Box[1]->min.setX(this->center.x());
	Box[1]->min.setY(this->min.y());
	Box[1]->min.setZ(this->min.z());
	Box[1]->max.setX(this->max.x());
	Box[1]->max.setY(this->center.y());
	Box[1]->max.setZ(this->center.z());
	Box[1]->calculate_center();
	this->Box[2] = new C_Box;
	Box[2]->min.setX(this->min.x());
	Box[2]->min.setY(this->center.y());
	Box[2]->min.setZ(this->min.z());
	Box[2]->max.setX(this->center.x());
	Box[2]->max.setY(this->max.y());
	Box[2]->max.setZ(this->center.z());
	Box[2]->calculate_center();
	this->Box[3] = new C_Box;
	Box[3]->min.setX(this->min.x());
	Box[3]->min.setY(this->min.y());
	Box[3]->min.setZ(this->center.z());
	Box[3]->max.setX(this->center.x());
	Box[3]->max.setY(this->center.y());
	Box[3]->max.setZ(this->max.z());
	Box[3]->calculate_center();
	this->Box[4] = new C_Box;
	Box[4]->min.setX(this->center.x());
	Box[4]->min.setY(this->center.y());
	Box[4]->min.setZ(this->min.z());
	Box[4]->max.setX(this->max.x());
	Box[4]->max.setY(this->max.y());
	Box[4]->max.setZ(this->center.z());
	Box[4]->calculate_center();
	this->Box[5] = new C_Box;
	Box[5]->min.setX(this->center.x());
	Box[5]->min.setY(this->min.y());
	Box[5]->min.setZ(this->center.z());
	Box[5]->max.setX(this->max.x());
	Box[5]->max.setY(this->center.y());
	Box[5]->max.setZ(this->max.z());
	Box[5]->calculate_center();
	this->Box[6] = new C_Box;
	Box[6]->min.setX(this->min.x());
	Box[6]->min.setY(this->center.y());
	Box[6]->min.setZ(this->center.z());
	Box[6]->max.setX(this->center.x());
	Box[6]->max.setY(this->max.y());
	Box[6]->max.setZ(this->max.z());
	Box[6]->calculate_center();
	this->Box[7] = new C_Box;
	Box[7]->min.setX(this->center.x());
	Box[7]->min.setY(this->center.y());
	Box[7]->min.setZ(this->center.z());
	Box[7]->max.setX(this->max.x());
	Box[7]->max.setY(this->max.y());
	Box[7]->max.setZ(this->max.z());
	Box[7]->calculate_center();
	for (int b = 0; b != 8; b++){
		Box[b]->T1s.clear();
		Box[b]->T2s.clear();
		Box[b]->N1s.clear();
		Box[b]->N2s.clear();
	}
}

void C_Box::split_tri(C_Line * IntSegments){
	int coplanar;
	C_Vector3D isectpt1, isectpt2;
	this->generate_subboxes();

	for (int t1 = 0; t1 != this->T1s.length(); t1++){
		for (int b = 0; b != 8; b++){
			if (Box[b]->tri_in_box(this->T1s[t1])){
				Box[b]->T1s.append(this->T1s[t1]);
			}
		}
	}

	for (int t2 = 0; t2 != this->T2s.length(); t2++){
		for (int b = 0; b != 8; b++){
			if (Box[b]->tri_in_box(this->T2s[t2])){
				Box[b]->T2s.append(this->T2s[t2]);
			}
		}
	}


	for (int b = 0; b != 8; b++){
		if ((Box[b]->T1s.length() != 0) && (Box[b]->T2s.length() != 0)){
			if (Box[b]->too_much_tri()){
				Box[b]->split_tri(IntSegments);
			}
			else{
				for (int t1 = 0; t1 != Box[b]->T1s.length(); t1++){
					for (int t2 = 0; t2 != Box[b]->T2s.length(); t2++){
						if (tri_tri_intersect_with_isectline(*Box[b]->T1s[t1], *Box[b]->T2s[t2], &coplanar, &isectpt1, &isectpt2) != 0){
							IntSegments->appendNonExistingSegment(isectpt1, isectpt2);
						}
					}
				}
				delete Box[b];
			}
		}else{
			delete Box[b];
		}
	}
}

void C_Box::split_seg(QList<C_Vector3D*> * TPs, int I1, int I2){
	C_Vector3D * TP;

	this->generate_subboxes();

	for (int n1 = 0; n1 < this->N1s.length()-1; n1+=2){
		for (int b = 0; b != 8; b++){
			if (Box[b]->seg_in_box(this->N1s[n1], this->N1s[n1+1])){
				Box[b]->N1s.append(this->N1s[n1]);
				Box[b]->N1s.append(this->N1s[n1+1]);
			}
		}
	}

	for (int n2 = 0; n2 < this->N2s.length()-1; n2+=2){
		for (int b = 0; b != 8; b++){
			if (Box[b]->seg_in_box(this->N2s[n2], this->N2s[n2+1])){
				Box[b]->N2s.append(this->N2s[n2]);
				Box[b]->N2s.append(this->N2s[n2+1]);
			}
		}
	}

	for (int b = 0; b != 8; b++){
		if ((Box[b]->N1s.length() != 0) && (Box[b]->N2s.length() != 0)){
			if (Box[b]->too_much_seg()){
				Box[b]->split_seg(TPs, I1, I2);
			}
			else{
				for (int n1 = 0; n1 < Box[b]->N1s.length()-1; n1+=2){
					for (int n2 = 0; n2 < Box[b]->N2s.length()-1; n2+=2){
						C_Line connector = connector.calculateSkewLineTransversal(*Box[b]->N1s[n1], *Box[b]->N1s[n1 + 1], *Box[b]->N2s[n2], *Box[b]->N2s[n2 + 1]);
						if (connector.Ns.length() == 2 && lengthSquared(connector.Ns[0] - connector.Ns[1])<1e-24){
							TP = new C_Vector3D((connector.Ns[0] + connector.Ns[1]) / 2);
							TP->intID = I1;

							/* The following statement is not thread-safe but
							 * will be executed concurrently by multiple threads.
							 * Therefore, it is important to guard this statement
							 * with a mutex. */
							mutex.lock();
							TPs->append(TP);
							mutex.unlock();

							TP = new C_Vector3D((connector.Ns[0] + connector.Ns[1]) / 2);
							TP->intID = I2;

							/* Protect list access against race-conditions. */
							mutex.lock();
							TPs->append(TP);
							mutex.unlock();
						}
					}
				}
				delete Box[b];
			}
		}else{
			delete Box[b];
		}
	}
}


void C_Model::calculate_int_polyline(int s1, int s2){
	C_Line IntSegments;
	int coplanar;
	C_Vector3D isectpt1, isectpt2;
	C_Box Box;

	const C_Surface& cs1 = this->Surfaces[s1];
	const C_Surface& cs2 = this->Surfaces[s2];

	if (cs1.min.x() < cs2.min.x()) Box.min.setX(cs1.min.x()); else Box.min.setX(cs2.min.x());
	if (cs1.max.x() > cs2.max.x()) Box.max.setX(cs1.max.x()); else Box.max.setX(cs2.max.x());
	if (Box.min.x() > Box.max.x()) return;
	if (cs1.min.y() < cs2.min.y()) Box.min.setY(cs1.min.y()); else Box.min.setY(cs2.min.y());
	if (cs1.max.y() > cs2.max.y()) Box.max.setY(cs1.max.y()); else Box.max.setY(cs2.max.y());
	if (Box.min.y() > Box.max.y()) return;
	if (cs1.min.z() < cs2.min.z()) Box.min.setZ(cs1.min.z()); else Box.min.setZ(cs2.min.z());
	if (cs1.max.z() > cs2.max.z()) Box.max.setZ(cs1.max.z()); else Box.max.setZ(cs2.max.z());
	if (Box.min.z() > Box.max.z()) return;

	for (int t1 = 0; t1 != cs1.Ts.length(); t1++){
		if (Box.tri_in_box(&cs1.Ts[t1])) Box.T1s.append(&cs1.Ts[t1]);
	}

	for (int t2 = 0; t2 != cs2.Ts.length(); t2++){
		if (Box.tri_in_box(&cs2.Ts[t2])) Box.T2s.append(&cs2.Ts[t2]);
	}

	if ((Box.T1s.length() != 0) && (Box.T2s.length() != 0)) Box.split_tri(&IntSegments);

	while (IntSegments.Ns.length()!=0){
		C_Line newInt;
		newInt.GenerateFirstSplineOfSegments(&IntSegments);
		newInt.Object[0]=s1; 
		newInt.Object[1]=s2;
		newInt.calculate_min_max();

		/* Protect list accesses against race-conditions. */
		mutex.lock();
		Intersections.append(newInt);
		Surfaces[s1].Intersections.append(&Intersections.last());   
		Surfaces[s2].Intersections.append(&Intersections.last());   
		mutex.unlock();
	}
}

bool C_Model::PointOnTriangle(const C_Vector3D& point, const C_Triangle& tri) const
{
	C_Vector3D v0 = C_Vector3D(tri.Ns[1]->x() - tri.Ns[0]->x(), tri.Ns[1]->y() - tri.Ns[0]->y(), tri.Ns[1]->z() - tri.Ns[0]->z());
	C_Vector3D v1 = C_Vector3D(tri.Ns[2]->x() - tri.Ns[0]->x(), tri.Ns[2]->y() - tri.Ns[0]->y(), tri.Ns[2]->z() - tri.Ns[0]->z());
	C_Vector3D v2 = C_Vector3D(point.x() - tri.Ns[0]->x(), point.y() - tri.Ns[0]->y(), point.z() - tri.Ns[0]->z());
	double dot00 = dot(v0, v0);
	double dot01 = dot(v0, v1);
	double dot02 = dot(v0, v2);
	double dot11 = dot(v1, v1);
	double dot12 = dot(v1, v2);
	double d_1 = 1.0 / (dot00*dot11 - dot01*dot01);
	double u = (dot11 * dot02 - dot01*dot12) * d_1;
	double v = (dot00 * dot12 - dot01*dot02) * d_1;
	return (u >= 0.0) && (v >= 0.0) && (u + v < 1.0);
}

void C_Model::calculate_int_point(int p, int s){
	QList<C_Vector3D*> Segments;
	QList<C_Triangle*> TriMesh;
	C_Line newInt;
	C_Vector3D isectpt;
	C_Vector3D point_minus;
	C_Vector3D point_plus;
	if (Polylines[p].Path.Ns.length() < 2)
	{
		if (PreTestIntersectionsPolylineSurface(&this->Polylines[p], &this->Surfaces[s]))
		{
			point_minus = Polylines[p].Path.Ns[0] - C_Vector3D(1e-12, 0.0, 1e-12);
			point_plus = Polylines[p].Path.Ns[0] + C_Vector3D(1e-12, 0.0, +1e-12);
			if (this->PreTestIntersectionsSegmentSurface(&point_minus, &point_plus, &Surfaces[s]))
			{
				for (int t = 0; t != Surfaces[s].Ts.length(); t++)
				{
					if (this->PreTestIntersectionsTrianglePolyline(&Surfaces[s].Ts[t], &Polylines[p]))
					{
						if (PointOnTriangle(Polylines[p].Path.Ns[0], Surfaces[s].Ts[t]))
						{
							isectpt.setType("INTERSECTION_POINT");
							newInt.Ns.append(Polylines[p].Path.Ns[0]);
							newInt.Object[0] = s;
							newInt.Object[1] = p;

							/* Protect list accesses against race-conditions. */
							mutex.lock();
							this->Intersections.append(newInt);
							Polylines[p].Intersections.append(&this->Intersections.last());
							Surfaces[s].Intersections.append(&this->Intersections.last());
							mutex.unlock();

							return;
						}
					}
				}
			}
		}

	}
	else
	{
		if (PreTestIntersectionsPolylineSurface(&this->Polylines[p], &this->Surfaces[s])){
			for (int n = 0; n != Polylines[p].Path.Ns.length() - 1; n++){
				if (this->PreTestIntersectionsSegmentSurface(&this->Polylines[p].Path.Ns[n], &this->Polylines[p].Path.Ns[n + 1], &Surfaces[s])){
					Segments.append(&this->Polylines[p].Path.Ns[n]);
					Segments.append(&this->Polylines[p].Path.Ns[n + 1]);
				}
			}
			for (int t = 0; t != Surfaces[s].Ts.length(); t++){
				if (this->PreTestIntersectionsTrianglePolyline(&Surfaces[s].Ts[t], &Polylines[p])) TriMesh.append(&Surfaces[s].Ts[t]);
			}
			for (int se = 0; se != Segments.length(); se = se + 2){
				for (int t = 0; t != TriMesh.length(); t++){
					if (PreTestIntersectionsSegmentTriangle(Segments[se], Segments[se + 1], TriMesh[t])){
						if (triangle_ray_intersection(*TriMesh[t], *Segments[se], *Segments[se + 1] - *Segments[se], &isectpt)){
							isectpt.setType("INTERSECTION_POINT");
							newInt.Ns.append(isectpt);
							newInt.Object[0] = s;
							newInt.Object[1] = p;

							/* Protect list accesses against race-conditions. */
							mutex.lock();
							this->Intersections.append(newInt);
							Polylines[p].Intersections.append(&this->Intersections.last());
							Surfaces[s].Intersections.append(&this->Intersections.last());
							mutex.unlock();

							return;
						}
					}
				}
			}
		}
		Segments.clear();
		TriMesh.clear();
	}
}

void C_Model::insert_int_triplepoints(){
	//cleaning TPs
	QList<C_Vector3D*> memTPs;
	bool insert;
	for (int t = 0; t != this->TPs.length(); t++){
		insert = true;
		for (int m = 0; m != memTPs.length(); m++){
			if (memTPs[m]->intID == TPs[t]->intID){
				if (lengthSquared(*memTPs[m] - *TPs[t])<1e-24){
					insert = false;
					break;
				}
			}
		}
		if (insert) memTPs.append(this->TPs[t]);
	}
	TPs = memTPs;

	//inserting to intersection
	for (int t=0;t!=this->TPs.length();t++){
		this->TPs[t]->setType("TRIPLE_POINT");
		Intersections[this->TPs[t]->intID].AddPoint(*this->TPs[t]); 
	}
	for (int i=0;i!=this->Intersections.length();i++){
		this->Intersections[i].CleanIdenticalPoints(); 
		this->Intersections[i].AddPosition();
		this->Intersections[i].RefineByLength(this->Intersections[i].size);
	}
}

/*! \brief Calculates the intersection point of two polylines and added this point to both.
	\details 
		Go through all intersection polylines and calculates the shortest transversal between the polyline parts, by calling 
		C_Line::calculateSkewLineTransversal(const C_Vector3D& P1, const C_Vector3D& P2, const C_Vector3D& Q1, const C_Vector3D& Q2).
		If the length of the shortest transversal is smaller than 1e-12, then the average of the transversal start and end point are added to the belonging polyline segments.
		These points are marked as "TRIPLE_POINT"..

*/
void C_Model::calculate_int_triplepoints(int I1, int I2){

	C_Box Box;
	const QList<C_Line>& Intersecs = Intersections;
	const C_Line& line1 = Intersecs[I1];
	const C_Line& line2 = Intersecs[I2];

	if (line1.min.x() < line2.min.x()) Box.min.setX(line2.min.x()); else Box.min.setX(line1.min.x());
	if (line1.max.x() > line2.max.x()) Box.max.setX(line2.max.x()); else Box.max.setX(line1.max.x());
	if (Box.min.x() > Box.max.x()) return;
	if (line1.min.y() < line2.min.y()) Box.min.setY(line2.min.y()); else Box.min.setY(line1.min.y());
	if (line1.max.y() > line2.max.y()) Box.max.setY(line2.max.y()); else Box.max.setY(line1.max.y());
	if (Box.min.y() > Box.max.y()) return;
	if (line1.min.z() < line2.min.z()) Box.min.setZ(line2.min.z()); else Box.min.setZ(line1.min.z());
	if (line1.max.z() > line2.max.z()) Box.max.setZ(line2.max.z()); else Box.max.setZ(line1.max.z());
	if (Box.min.z() > Box.max.z()) return;
	Box.min -= C_Vector3D(1e-12, 1e-12, 1e-12);
	Box.max += C_Vector3D(1e-12, 1e-12, 1e-12);

	for (int n1 = 0; n1 != line1.Ns.length() - 1; n1++){
		if (Box.seg_in_box(&line1.Ns[n1], &line1.Ns[n1 + 1])){
			Box.N1s.append(&(line1.Ns[n1]));
			Box.N1s.append(&line1.Ns[n1 + 1]);
		}
	}

	for (int n2 = 0; n2 != line2.Ns.length() - 1; n2++){
		if (Box.seg_in_box(&line2.Ns[n2], &line2.Ns[n2 + 1])){
			Box.N2s.append(&line2.Ns[n2]);
			Box.N2s.append(&line2.Ns[n2 + 1]);
		}
	}

	if ((Box.N1s.length() != 0) && (Box.N2s.length() != 0)) Box.split_seg(&this->TPs, I1, I2);
}

void C_Model::analyze_self_intersections(const tetgenio& out)
{
	selfIntersections.clear();

	for(int i = 0; i < out.numberoftrifaces*3; i+=3)
	{
		C_Triangle tri;
		QSet<const C_Surface*> surfaces;

		for(int j = 0; j != Surfaces.length(); j++)
		{
			for(int k = 0; k < Surfaces[j].Ns.length(); k++)
			{
				C_Vector3D& p = Surfaces[j].Ns[k];
				if( p.tetID == out.trifacelist[i] )
					tri.Ns[0] = &p;
				else if( p.tetID == out.trifacelist[i+1] )
					tri.Ns[1] = &p;
				else if( p.tetID == out.trifacelist[i+2] )
					tri.Ns[2] = &p;
				else
					continue;

				surfaces.insert(&Surfaces[j]);
			}
		}

		if( tri.Ns[0] && tri.Ns[1] && tri.Ns[2] )
			selfIntersections.append(std::make_pair(tri, surfaces));
	}
}

void C_Model::calculate_tets(QString switches){
	QList <C_Vector3D> points;
	QList<double> pointlist;
	QList<int> pointtetIDlist;
	C_Vector3D diff;
	int tetID=0;
	int mem_tetID;

	for (int s = 0; s != this->Surfaces.length(); s++){
		for (int n = 0; n != this->Surfaces[s].duplicates; n++){
			mem_tetID = -1;
			for (int p = 0; p != pointtetIDlist.length(); p++){
				if (FABS(Surfaces[s].Ns[n].x() - pointlist[p * 3 + 0])<1e-12 && FABS(Surfaces[s].Ns[n].y() - pointlist[p * 3 + 1])<1e-12 && FABS(Surfaces[s].Ns[n].z() - pointlist[p * 3 + 2])<1e-12){
					mem_tetID = pointtetIDlist[p];
					break;
				}
			}
			if (mem_tetID == -1){
				pointlist.append(Surfaces[s].Ns[n].x());
				pointlist.append(Surfaces[s].Ns[n].y());
				pointlist.append(Surfaces[s].Ns[n].z());
				pointtetIDlist.append(tetID);
				Surfaces[s].Ns[n].tetID = tetID++;
			}
			else{
				Surfaces[s].Ns[n].tetID = mem_tetID;
			}
		}
	}
	for (int s = 0; s != this->Surfaces.length(); s++){
		for (int n = this->Surfaces[s].duplicates; n != this->Surfaces[s].Ns.length(); n++){
			pointlist.append(Surfaces[s].Ns[n].x());
			pointlist.append(Surfaces[s].Ns[n].y());
			pointlist.append(Surfaces[s].Ns[n].z());
			pointtetIDlist.append(tetID);
			Surfaces[s].Ns[n].tetID = tetID++;
		}
	}

	for (int po=0;po!=Polylines.length();po++){
		for (int n=0;n!=Polylines[po].Path.Ns.length();n++){
			mem_tetID=-1;
			for (int p=0;p!=pointtetIDlist.length();p++){
				if (FABS(Polylines[po].Path.Ns[n].x()-pointlist[p*3+0])<1e-12 && FABS(Polylines[po].Path.Ns[n].y()-pointlist[p*3+1])<1e-12 && FABS(Polylines[po].Path.Ns[n].z()-pointlist[p*3+2])<1e-12){
					mem_tetID=pointtetIDlist[p];
					break;
				}
			}
			if (mem_tetID==-1){
				pointlist.append(Polylines[po].Path.Ns[n].x()); 
				pointlist.append(Polylines[po].Path.Ns[n].y()); 
				pointlist.append(Polylines[po].Path.Ns[n].z()); 
				pointtetIDlist.append(tetID); 
				Polylines[po].Path.Ns[n].tetID=tetID++;
			}else{
				Polylines[po].Path.Ns[n].tetID=mem_tetID; 
			}
		}
	}

	int number_faces=0;
	for (int s=0;s!=this->Surfaces.length();s++){
		number_faces+=Surfaces[s].Ts.length();
	}
	int number_edges=0;
	for (int p=0;p!=Polylines.length();p++){
		if (this->Polylines[p].Path.Ns.length()>=2){
			number_edges+=Polylines[p].Path.Ns.length()-1;
		}
	}
	tetgenio in, out, tmp;
	tetgenio::facet *f;
	tetgenio::polygon *p;
// All indices start from 0
	in.firstnumber = 0;
	in.numberofpoints = pointtetIDlist.length();
	in.pointlist = new REAL[in.numberofpoints * 3];

	for (int p=0;p!=pointtetIDlist.length();p++){
		in.pointlist[p*3+0]=pointlist[p*3+0];
		in.pointlist[p*3+1]=pointlist[p*3+1];
		in.pointlist[p*3+2]=pointlist[p*3+2];
	}

	in.numberoffacets = number_faces;
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];

	int face_number=0;
	for (int s=0;s!=this->Surfaces.length();s++){
		for (int t=0;t!=Surfaces[s].Ts.length();t++){
			f = &in.facetlist[face_number++];
			f->numberofpolygons = 1;
			f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
			f->numberofholes = 0;
			f->holelist = NULL;
			p = &f->polygonlist[0];
			p->numberofvertices = 3;
			p->vertexlist = new int[p->numberofvertices];
			p->vertexlist[0] = Surfaces[s].Ts[t].Ns[0]->tetID;
			p->vertexlist[1] = Surfaces[s].Ts[t].Ns[1]->tetID;
			p->vertexlist[2] = Surfaces[s].Ts[t].Ns[2]->tetID;
		}
	}

	if (face_number==0) return;
	face_number=0;
	for (int s=0;s!=this->Surfaces.length();s++){
		for (int t=0;t!=Surfaces[s].Ts.length();t++){
			in.facetmarkerlist[face_number++] = s;
		}
	}

	in.numberofedges =number_edges;
	in.edgelist  = new int[in.numberofedges*2];
	in.edgemarkerlist = new int[in.numberofedges];

	int edgelistcount=0;
	int edgemarkerlistcount=0;
	for (int p=0;p!=Polylines.length();p++){
		if (this->Polylines[p].Path.Ns.length()>=2){
			for (int n=0;n!=Polylines[p].Path.Ns.length()-1;n++){
				in.edgelist[edgelistcount++]=this->Polylines[p].Path.Ns[n].tetID;   
				in.edgelist[edgelistcount++]=this->Polylines[p].Path.Ns[n+1].tetID;   
				in.edgemarkerlist[edgemarkerlistcount++]=p+2; //currently "0" ad "1" are not working!!! bug of tetgen
			}
		}
	}

	in.numberofregions = 0;
	for (int m=0;m!=this->Mats.length();m++) in.numberofregions+=this->Mats[m].Locations.length();
	in.regionlist = new REAL[in.numberofregions * 5];
	int currentRegion=0;
	for (int m=0;m!=this->Mats.length();m++){
		for (int l=0;l!=this->Mats[m].Locations.length();l++){ 
			in.regionlist[currentRegion*5+0]=this->Mats[m].Locations[l].x();
			in.regionlist[currentRegion*5+1]=this->Mats[m].Locations[l].y(); 
			in.regionlist[currentRegion*5+2]=this->Mats[m].Locations[l].z();
			in.regionlist[currentRegion*5+3]=m; 
			in.regionlist[currentRegion*5+4]=0;
			currentRegion++;
		}
	}

	// Output the PLC to files 'barin.node' and 'barin.poly'.
	if( QFileInfo(QDir::currentPath()).isWritable() )
	{
		in.save_nodes(const_cast<char*>("in"));
		in.save_edges(const_cast<char*>("in"));
		in.save_poly(const_cast<char*>("in"));
	}

	// Tetrahedralize the PLC
	QString Attr=switches;
	try {
		tetrahedralize((char*)Attr.toLatin1().data(), &in, &out,NULL,NULL);
	} catch (int x) {
		printf("tetgen failed\n");
		switch (x) {
		case 1: // Out of memory.
			printf("Error:  Out of memory.\n");
			emit PrintError("Out of memory.");
			break;
		case 2: // Encounter an internal error.
			printf("Please report this bug to Hang.Si@wias-berlin.de. Include\n");
			printf("  the message above, your input data set, and the exact\n");
			printf("  command line you used to run this program, thank you.\n");
			emit PrintError("An internal error has occured.");
			break;
		case 3:
			printf("A self-intersection was detected. Program stopped.\n");
			printf("Hint: use -d option to detect all self-intersections.\n");
			emit PrintError("A self-intersection was detected.");
			emit PrintError("Analyzing all self-intersections...");
			/* Run tetgen once again (using -d commandline option) to detect
			 * all self-intersections. */
			Attr = Attr + "d";
			tetrahedralize((char*)Attr.toLatin1().data(), &in, &out, NULL, NULL);
			analyze_self_intersections(out);
			emit ErrorInfoChanged("A self-intersection was detected.");
			break;
		case 4:
			printf("A very small input feature size was detected. Program stopped.\n");
			emit PrintError("A very small input feature size was detected.");
			break;
		case 5:
			printf("Two very close input facets were detected. Program stopped.\n");
			printf("Hint: use -Y option to avoid adding Steiner points in boundary.\n");
			emit PrintError("Two very close input facets were detected.");
			break;
		case 10:
			printf("An input error was detected. Program stopped.\n");
			emit PrintError("Input error, please verify your input data.");
			break;
		  } // switch (x)
		return;
	}

	// Output mesh to files 'barout.node', 'barout.ele' and 'barout.face'.
	if( QFileInfo(QDir::currentPath()).isWritable() )
	{
		out.save_nodes(const_cast<char*>("out"));
		// out.save_elements("out");
		// out.save_faces("out");
		out.save_edges(const_cast<char*>("out"));
		// out.save_poly("out");
		out.save_faces2smesh(const_cast<char*>("out"));
	}

	if (this->Mesh){
		delete this->Mesh;
		this->Mesh = 0;
	}
	this->Mesh = new C_Mesh3D;

	Mesh->numberofpoints=out.numberofpoints;
	for (int p = 0; p < out.numberofpoints; p++){
		this->Mesh->pointlist.append(out.pointlist[p*3+0]); 
		this->Mesh->pointlist.append(out.pointlist[p*3+1]); 
		this->Mesh->pointlist.append(out.pointlist[p*3+2]); 
	}

	for (int e = 0; e < out.numberofedges; e++){
		if (out.edgemarkerlist[e]>=2 && out.edgemarkerlist[e]<(this->Polylines.length()+2)){  //currently "0" ad "1" are not working!!! bug of tetgen 
			this->Mesh->edgelist.append(out.edgelist[2*e+0]); 
			this->Mesh->edgelist.append(out.edgelist[2*e+1]); 
			this->Mesh->edgemarkerlist.append(out.edgemarkerlist[e]-2);   //currently "0" ad "1" are not working!!! bug of tetgen 
		}
	}
	Mesh->numberofedges=this->Mesh->edgemarkerlist.length();

	for (int f = 0; f < out.numberoftrifaces; f++){
		if (out.trifacemarkerlist[f]>=0 && out.trifacemarkerlist[f]<this->Surfaces.length()){
			this->Mesh->trianglelist.append(out.trifacelist[3 * f + 0]);
			this->Mesh->trianglelist.append(out.trifacelist[3*f+1]);  
			this->Mesh->trianglelist.append(out.trifacelist[3*f+2]);  
			this->Mesh->trianglemarkerlist.append(out.trifacemarkerlist[f]);
		}
	}
	Mesh->numberoftriangles=this->Mesh->trianglemarkerlist.length();

	for (int t = 0; t < out.numberoftetrahedra; t++){
		if (out.tetrahedronattributelist[t]>=0 && out.tetrahedronattributelist[t]<this->Mats.length()){
			this->Mesh->tetrahedronlist.append(out.tetrahedronlist[4*t+0]);  
			this->Mesh->tetrahedronlist.append(out.tetrahedronlist[4*t+1]);  
			this->Mesh->tetrahedronlist.append(out.tetrahedronlist[4*t+2]);  
			this->Mesh->tetrahedronlist.append(out.tetrahedronlist[4*t+3]);
			this->Mesh->tetrahedronmarkerlist.append(out.tetrahedronattributelist[t]);
		}
	}
	Mesh->numberoftetrahedra=this->Mesh->tetrahedronmarkerlist.length();

	emit ModelInfoChanged();
}

/* Combine polyline and surface constraints into a single list. */
void C_Model::get_constraints(std::list<C_Line*>& res)
{
	for (int s = 0; s != Surfaces.length(); s++)
		for (int c = 0; c != Surfaces[s].Constraints.length(); c++)
			res.push_back(&Surfaces[s].Constraints[c]);

	for (int p = 0; p != Polylines.length(); p++)
		for (int c = 0; c != Polylines[p].Constraints.length(); c++)
			res.push_back(&Polylines[p].Constraints[c]);
}

/* Check whether at least one constraint is marked. */
bool C_Model::has_selected_constraints()
{
	std::list<C_Line*> constraints;
	get_constraints(constraints);

	std::list<C_Line*>::iterator it;
	for(it = constraints.begin(); it != constraints.end(); ++it)
	{
		C_Line& line = **it;
		if( line.Type != "UNDEFINED" )
			return true;
	}

	return false;
}

/* Mark all constraints. */
void C_Model::set_all_constraints(const QString& type)
{
	std::list<C_Line*> constraints;
	get_constraints(constraints);

	std::list<C_Line*>::iterator it;
	for(it = constraints.begin(); it != constraints.end(); ++it)
		(*it)->Type = type;
}

/* Mark all constraints. */
void C_Model::select_all_constraints()
{
	return set_all_constraints("SEGMENTS");
}

void C_Model::deselect_all_constraints()
{
	return set_all_constraints("UNDEFINED");
}

/* Apply selections based on the assigned materials. This requires a preceding
 * tetgen run (without refinement switch -q). */
void C_Model::material_selections()
{
	std::list<C_Line*> constraints;
	get_constraints(constraints);

	/* Walk through all the constraints and check if the generated materials
	 * reach all the line segments. */
	std::list<C_Line*>::iterator it;
	for(it = constraints.begin(); it != constraints.end(); ++it)
	{
		C_Line& line = **it;

		bool has_line = false;

		QList<C_Vector3D>::const_iterator pit;
		for(pit = line.Ns.cbegin(); pit != line.Ns.cend(); ++pit)
		{
			bool has_match = false;

			/* Check triangles. */
			for(int f = 0; f < this->Mesh->numberoftriangles && !has_match; f++)
			{
				if( getMaterial(5, f) == -1 )
					continue;

				for(int p = 0; p != 3; p++)
				{
					double point[3];
					Mesh->getCoordinates(Mesh->trianglelist[f*3+p], point);

					if( fabs(point[0] - pit->x()) < 1e-10 &&
					    fabs(point[1] - pit->y()) < 1e-10 &&
					    fabs(point[2] - pit->z()) < 1e-10 )
					{
						has_match = true;
						break;
					}
				}
			}

			/* Check tetrahedrons. */
			for(int t = 0; t < this->Mesh->numberoftetrahedra && !has_match; t++)
			{
				if( getMaterial(10, t) == -1 )
					continue;

				for(int p = 0; p != 4; p++)
				{
					double point[3];
					Mesh->getCoordinates(Mesh->tetrahedronlist[t*4+p], point);

					if( fabs(point[0] - pit->x()) < 1e-10 &&
					    fabs(point[1] - pit->y()) < 1e-10 &&
					    fabs(point[2] - pit->z()) < 1e-10 )
					{
						has_match = true;
						break;
					}
				}
			}

			if( ! has_match )
			{
				has_line = false;
				break;
			}

			has_line = true;
		}

		/* Mark line accordingly. */
		line.Type = has_line ? "SEGMENTS" : "UNDEFINED";
	}
}

bool C_Model::verify_materials()
{
	bool state = true;
	for(int i = 0; i < Mats.size(); ++i)
	{
		for(int j = 0; j < Mats[i].Locations.size(); ++j)
		{
			const C_Vector3D& pos = Mats[i].Locations[j];
			if( ! Mesh || ! Mesh->isPointInsideAnyMaterial(pos) )
			{
				std::stringstream ss;
				ss << "Material point lies outside the meshed area: " <<
				      "(Material: " << i << ", Location: " << j << ").";
				emit PrintError(QString::fromStdString(ss.str()));
				state = false;
			}
		}
	}

	return state;
}

void C_Model::makeTets(bool xCutEnable, double xCutValue, bool xDirection, bool yCutEnable, double yCutValue, bool yDirection, bool zCutEnable, double zCutValue, bool zDirection){

	int xDir,yDir,zDir;
	int mat;
	bool drawThis;
	double point[4][3];
	double center[3];
	double normal[3];

	if (xDirection){xDir=1;}else{xDir=-1;}
	if (yDirection){yDir=1;}else{yDir=-1;}
	if (zDirection){zDir=1;}else{zDir=-1;}
	xCutValue/=256.0;
	yCutValue/=256.0;
	zCutValue/=256.0;

	glDeleteLists(this->listTets,1); 

	/* Clear last error. */
	glGetError();

	GLuint list = glGenLists(1);
	glNewList(list, GL_COMPILE);
	glDisable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.White);

	glBegin(GL_LINES);
	for (int f = 0; f!=this->Mesh->numberoftriangles;f++){
		mat = getMaterial(5, f);
		drawThis=true;
		if (!this->Surfaces[this->Mesh->trianglemarkerlist[f]].drawMatEdges) drawThis=false;
		if (mat==-1) drawThis=false;
		this->Mesh->getCenterOfTriangle(f,center); 
		if (xCutEnable && (center[0]-xCutValue)*xDir<0) drawThis=false;
		if (yCutEnable && (center[1]-yCutValue)*yDir<0) drawThis=false;
		if (zCutEnable && (center[2]-zCutValue)*zDir<0) drawThis=false;
		if (drawThis){
			int mod=mat%6; 
			if (mod==0) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Red);
			if (mod==1) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Green);
			if (mod==2) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Blue);
			if (mod==3) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Yellow);
			if (mod==4) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Cyan);
			if (mod==5) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Magenta);
			for (int p=0;p!=3;p++) this->Mesh->getCoordinates(this->Mesh->trianglelist[f*3+p],point[p]); 
			glVertex3d(point[0][0],point[0][1],point[0][2]);
			glVertex3d(point[1][0],point[1][1],point[1][2]);
			glVertex3d(point[1][0],point[1][1],point[1][2]);
			glVertex3d(point[2][0],point[2][1],point[2][2]);
			glVertex3d(point[2][0],point[2][1],point[2][2]);
			glVertex3d(point[0][0],point[0][1],point[0][2]);
		}
	}
	for (int t = 0; t!=this->Mesh->numberoftetrahedra;t++){
		mat = getMaterial(10, t);
		drawThis = true;
		if (!this->Mats[this->Mesh->tetrahedronmarkerlist[t]].drawMatEdges) drawThis=false;
		if (mat==-1) drawThis=false;
		this->Mesh->getCenterOfTetrahedron(t,center); 
		if (xCutEnable && (center[0]-xCutValue)*xDir<0) drawThis=false;
		if (yCutEnable && (center[1]-yCutValue)*yDir<0) drawThis=false;
		if (zCutEnable && (center[2]-zCutValue)*zDir<0) drawThis=false;
		if (drawThis){
			int mod=mat%6; 
			if (mod==0) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Red);
			if (mod==1) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Green);
			if (mod==2) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Blue);
			if (mod==3) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Yellow);
			if (mod==4) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Cyan);
			if (mod==5) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Magenta);
			for (int p=0;p!=4;p++) this->Mesh->getCoordinates(this->Mesh->tetrahedronlist[t*4+p],point[p]); 
			glVertex3d(point[0][0],point[0][1],point[0][2]);
			glVertex3d(point[1][0],point[1][1],point[1][2]);
			glVertex3d(point[1][0],point[1][1],point[1][2]);
			glVertex3d(point[2][0],point[2][1],point[2][2]);
			glVertex3d(point[2][0],point[2][1],point[2][2]);
			glVertex3d(point[0][0],point[0][1],point[0][2]);
			glVertex3d(point[0][0],point[0][1],point[0][2]);
			glVertex3d(point[3][0],point[3][1],point[3][2]);
			glVertex3d(point[1][0],point[1][1],point[1][2]);
			glVertex3d(point[3][0],point[3][1],point[3][2]);
			glVertex3d(point[2][0],point[2][1],point[2][2]);
			glVertex3d(point[3][0],point[3][1],point[3][2]);
		}
	}
    glEnd();

	glEnable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.Grey);

	glBegin(GL_TRIANGLES);
	for (int f = 0; f!=this->Mesh->numberoftriangles;f++){
		mat = getMaterial(5, f);
		drawThis=true;
		if (!this->Surfaces[this->Mesh->trianglemarkerlist[f]].drawMatFaces) drawThis=false;
		if (mat==-1) drawThis=false;
		this->Mesh->getCenterOfTriangle(f,center); 
		if (xCutEnable && (center[0]-xCutValue)*xDir<0) drawThis=false;
		if (yCutEnable && (center[1]-yCutValue)*yDir<0) drawThis=false;
		if (zCutEnable && (center[2]-zCutValue)*zDir<0) drawThis=false;
		if (drawThis){
			int mod=mat%6; 
			if (mod==0) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.RedTrans);
			if (mod==1) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.GreenTrans);
			if (mod==2) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.BlueTrans);
			if (mod==3) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.YellowTrans);
			if (mod==4) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.CyanTrans);
			if (mod==5) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.MagentaTrans);
			this->Mesh->getNormalOfTriangle(f,normal);  
			glNormal3d(normal[0],normal[1],normal[2]);
			for (int p=0;p!=3;p++) this->Mesh->getCoordinates(this->Mesh->trianglelist[f*3+p],point[p]); 
			glVertex3d(point[0][0],point[0][1],point[0][2]);
			glVertex3d(point[1][0],point[1][1],point[1][2]);
			glVertex3d(point[2][0],point[2][1],point[2][2]);
		}
	}
	for (int t = 0; t!=this->Mesh->numberoftetrahedra;t++){
		mat = getMaterial(10, t);
		drawThis=true;
		if (!this->Mats[this->Mesh->tetrahedronmarkerlist[t]].drawMatFaces) drawThis=false;
		if (mat==-1) drawThis=false;
		this->Mesh->getCenterOfTetrahedron(t,center); 
		if (xCutEnable && (center[0]-xCutValue)*xDir<0) drawThis=false;
		if (yCutEnable && (center[1]-yCutValue)*yDir<0) drawThis=false;
		if (zCutEnable && (center[2]-zCutValue)*zDir<0) drawThis=false;
		if (drawThis){
			int mod=mat%6; 
			if (mod==0) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.RedTrans);
			if (mod==1) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.GreenTrans);
			if (mod==2) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.BlueTrans);
			if (mod==3) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.YellowTrans);
			if (mod==4) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.CyanTrans);
			if (mod==5) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.MagentaTrans);
			for (int p=0;p!=4;p++) this->Mesh->getCoordinates(this->Mesh->tetrahedronlist[t*4+p],point[p]); 
			this->Mesh->getNormalOfTetrahedron(t,0,normal);  
			glNormal3d(normal[0],normal[1],normal[2]);
			glVertex3d(point[0][0],point[0][1],point[0][2]);
			glVertex3d(point[1][0],point[1][1],point[1][2]);
			glVertex3d(point[2][0],point[2][1],point[2][2]);
			this->Mesh->getNormalOfTetrahedron(t,1,normal);  
			glNormal3d(normal[0],normal[1],normal[2]);
			glVertex3d(point[0][0],point[0][1],point[0][2]);
			glVertex3d(point[1][0],point[1][1],point[1][2]);
			glVertex3d(point[3][0],point[3][1],point[3][2]);
			this->Mesh->getNormalOfTetrahedron(t,2,normal);  
			glNormal3d(normal[0],normal[1],normal[2]);
			glVertex3d(point[1][0],point[1][1],point[1][2]);
			glVertex3d(point[2][0],point[2][1],point[2][2]);
			glVertex3d(point[3][0],point[3][1],point[3][2]);
			this->Mesh->getNormalOfTetrahedron(t,3,normal);  
			glNormal3d(normal[0],normal[1],normal[2]);
			glVertex3d(point[2][0],point[2][1],point[2][2]);
			glVertex3d(point[0][0],point[0][1],point[0][2]);
			glVertex3d(point[3][0],point[3][1],point[3][2]);
		}
	}
	glEnd();

	glEndList();
	listTets = list;

	std::string err_msg;
	if( ! check_opengl_error(err_msg) )
		emit PrintError(QString::fromStdString(err_msg));
}

int C_Model::getMaterial(int vtkType, long listID)
{
	switch (vtkType)
	{
	case 10:
		return this->Mesh->tetrahedronmarkerlist[listID];
		break;
	case 5:
		return this->Surfaces[this->Mesh->trianglemarkerlist[listID]].MaterialID;
		break;
	case 3:
		return this->Polylines[this->Mesh->edgemarkerlist[listID]].MaterialID;
		break;
	case 1:
		break;
	}
	return -1;
}

void C_Model::makeMats(int Material, int Location){
	double x,y,z;
	GLuint list = glGenLists(1);
	glNewList(list, GL_COMPILE);
	glDisable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.White);
	if (Material!=-1 && Location!=-1){
		x = this->Mats[Material].Locations[Location].x();
		y = this->Mats[Material].Locations[Location].y();
		z = this->Mats[Material].Locations[Location].z();

		glEnable(GL_LINE_STIPPLE);

		glBegin(GL_LINES);
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Yellow);
			glVertex3d(-2.0,y,z);
			glVertex3d(x-0.02,y,z);
			glVertex3d(x+0.02,y,z);
			glVertex3d(2.0,y,z);
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Red);
			glVertex3d(x,-2.0,z);
			glVertex3d(x,y-0.02,z);
			glVertex3d(x,y+0.02,z);
			glVertex3d(x,+2.0,z);
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Green);
			glVertex3d(x,y,-2.0);
			glVertex3d(x,y,z-0.02);
			glVertex3d(x,y,z+0.02);
			glVertex3d(x,y,+2.0);
		glEnd();

		glDisable(GL_LINE_STIPPLE);

		glBegin(GL_POINTS);
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Blue);
			glVertex3d(x,y,z);
		glEnd();
	}
	if (Material!=-1){
		for (int l=0;l!=this->Mats[Material].Locations.length();l++){
			if (l!=Location){
				glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.Red);
				glBegin(GL_POINTS);
					glVertex3d(this->Mats[Material].Locations[l].x(),this->Mats[Material].Locations[l].y(),this->Mats[Material].Locations[l].z());
				glEnd();
			}
			glWrite(QString::number(Material) + "-" + QString::number(l),this->Mats[Material].Locations[l].x()+0.02,this->Mats[Material].Locations[l].y(),this->Mats[Material].Locations[l].z(),0.02); 
		}
	}
	glEnable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.Grey);

	//here 3D
	glEndList();
	listMats = list;
}

void C_Model::makeMarkers(const C_Line* selectedLine, const QString& comp)
{
	GLuint list = glGenLists(1);
	glNewList(list, GL_COMPILE);

	glDisable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.White);

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, Cols.Orange);

	if( selectedLine )
	{
		const C_Line& line = *selectedLine;
		glBegin(GL_LINE_STRIP);
			for(int n = 0; n < line.Ns.length(); n++)
				glVertex3f(line.Ns[n].x(), line.Ns[n].y(), line.Ns[n].z());
		glEnd();

		if( line.Ns.length() == 1 )
		{
			glBegin(GL_POINTS);
			glVertex3f(line.Ns[0].x(), line.Ns[0].y(), line.Ns[0].z());
			glEnd();
		}
	}
	else
	{
		for(int i = 0; i < Surfaces.length(); ++i)
		{
			const C_Surface& surface = Surfaces[i];
			if( !comp.isEmpty() && surface.Name != comp )
				continue;

			for(int j = 0; j < surface.Constraints.length(); j++)
			{
				const C_Line& line = surface.Constraints[j];
				if( line.Type != "UNDEFINED" )
					continue;

				glBegin(GL_LINE_STRIP);
					for(int n = 0; n < line.Ns.length(); n++)
						glVertex3f(line.Ns[n].x(), line.Ns[n].y(),
						           line.Ns[n].z());
				glEnd();

				if( line.Ns.length() == 1 )
				{
					glBegin(GL_POINTS);
					glVertex3f(line.Ns[0].x(), line.Ns[0].y(), line.Ns[0].z());
					glEnd();
				}
			}
		}
	}

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.Grey);
	glEnable(GL_LIGHT0);

	glEndList();
	listMarkers = list;
}

void C_Model::makeErrorMarkers(const C_Line& line)
{
	GLuint list = glGenLists(1);
	glNewList(list, GL_COMPILE);

	glDisable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.White);

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, Cols.Orange);

	glBegin(GL_LINE_STRIP);
		for(int n = 0; n < line.Ns.length(); n++)
			glVertex3f(line.Ns[n].x(), line.Ns[n].y(), line.Ns[n].z());
	glEnd();

	if( line.Ns.length() == 1 )
	{
		glBegin(GL_POINTS);
		glVertex3f(line.Ns[0].x(), line.Ns[0].y(), line.Ns[0].z());
		glEnd();
	}

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Cols.Grey);
	glEnable(GL_LIGHT0);

	glEndList();
	listErrorMarkers = list;
}

void C_Model::glWrite(QString string, double x,double y,double z,double scale){
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Cols.White);
	glBegin(GL_LINES);
	for (int s=0; s!=string.length();s++){
		if (string[s]=='x'){
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
		}
		if (string[s]=='y'){
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.0*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.0*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x+0.0*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x+0.0*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
		}
		if (string[s]=='z'){
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
		}
		if (string[s]=='-'){
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
		}
		if (string[s]=='0'){
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
		}
		if (string[s]=='1'){
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
		}
		if (string[s]=='2'){
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
		}
		if (string[s]=='3'){
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z-0.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-0.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
		}
		if (string[s]=='4'){
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
		}
		if (string[s]=='5'){
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
		}
		if (string[s]=='6'){
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
		}
		if (string[s]=='7'){
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
		}
		if (string[s]=='8'){
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
		}
		if (string[s]=='9'){
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+0.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z+1.0*scale);
			glVertex3d(x+0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
			glVertex3d(x-0.5*scale+s*1.5*scale, y+0.0*scale, z-1.0*scale);
		}
	}
	glEnd();
}

C_Vector3D C_Model::getGeomCenter()
{
	C_Vector3D geom(0, 0, 0);
	std::list<const QList<C_Vector3D>*> objects;

	for(int s = 0; s != Surfaces.length(); s++)
	{
		const C_Surface& surface = Surfaces[s];

		if( surface.drawScatteredData ||
		    surface.drawConvexHull ||
		    surface.drawFaces ||
		    surface.drawEdges )
		{
			if( surface.SDs.length() > 0 )
				objects.push_back(&surface.SDs);
		}

		if( surface.drawIntEdges ||
		    surface.drawIntVertices ||
		    surface.drawConstraints )
		{
			for(int n = 0; n < surface.Intersections.length(); n++)
			{
				const C_Line *line = surface.Intersections[n];
				if( line->Ns.length() > 0 )
					objects.push_back(&line->Ns);
			}
		}

		if( surface.drawConstraints )
		{
			for(int n = 0; n < surface.Constraints.length(); n++)
			{
				const C_Line *line = &surface.Constraints[n];
				if( line->Ns.length() > 0 )
					objects.push_back(&line->Ns);
			}
		}
	}

	for(int p = 0; p != Polylines.length(); p++)
	{
		const C_Polyline& polyline = Polylines[p];

		if( polyline.drawScatteredData ||
		    polyline.drawEdges ||
		    polyline.drawVertices)
		{
			if( polyline.SDs.length() > 0 )
				objects.push_back(&polyline.SDs);
		}

		if( polyline.drawIntVertices )
		{
			for(int n = 0; n < polyline.Intersections.length(); n++)
			{
				const C_Line *line = polyline.Intersections[n];
				if( line->Ns.length() > 0 )
					objects.push_back(&line->Ns);
			}
		}

		if( polyline.drawConstraints )
		{
			for(int n = 0; n < polyline.Constraints.length(); n++)
			{
				const C_Line *line = &polyline.Constraints[n];
				if( line->Ns.length() > 0 )
					objects.push_back(&line->Ns);
			}
		}
	}

	std::list<const QList<C_Vector3D>*>::const_iterator it;
	for(it = objects.cbegin(); it != objects.cend(); ++it)
	{
		const QList<C_Vector3D>& points = **it;

		C_Vector3D center(0, 0, 0);
		for(int n = 0; n < points.length(); n++)
			center += points[n];

		center /= points.length();
		geom += center;
	}

	if( objects.size() == 0 )
		return geom;

	return geom / objects.size();
}

bool C_Tetrahedron::is_on_same_side(const C_Vector3D& v1, const C_Vector3D& v2,
                                    const C_Vector3D& v3, const C_Vector3D& v4,
                                    const C_Vector3D &p) const
{
	C_Vector3D normal;
	cross(v2 - v1, v3 - v1, &normal);

	return sign(dot(normal, v4 - v1)) == sign(dot(normal, p - v1));
}

bool C_Tetrahedron::encloses_point(const C_Vector3D &p) const
{
	const C_Vector3D& v1 = Ns[0];
	const C_Vector3D& v2 = Ns[1];
	const C_Vector3D& v3 = Ns[2];
	const C_Vector3D& v4 = Ns[3];

	return is_on_same_side(v1, v2, v3, v4, p) &&
	       is_on_same_side(v2, v3, v4, v1, p) &&
	       is_on_same_side(v3, v4, v1, v2, p) &&
	       is_on_same_side(v4, v1, v2, v3, p);
}

void C_Mesh3D::getCoordinates(long pointnumber, C_Vector3D& point) const
{
	double p[3];
	getCoordinates(pointnumber, p);
	point.setX(p[0]);
	point.setY(p[1]);
	point.setZ(p[2]);
}

void C_Mesh3D::getCoordinates(long pointnumber, double * point) const {
	point[0]=this->pointlist[pointnumber*3+0];
	if (fabs(point[0]) < 1e-12) point[0] = 0;
	point[1]=this->pointlist[pointnumber*3+1]; 
	if (fabs(point[1]) < 1e-12) point[1] = 0;
	point[2]=this->pointlist[pointnumber*3+2];
	if (fabs(point[2]) < 1e-12) point[2] = 0;
}

void C_Mesh3D::getNormalOfTriangle(long trianglenumber, double * normal){
	double point[3][3];
	double diff[2][3];
	double cross[3];
	double norm;
	for (int p=0;p!=3;p++) this->getCoordinates(this->trianglelist[trianglenumber*3+p], point[p]); 
	for (int c=0;c!=3;c++) diff[0][c]=point[1][c]-point[0][c];
	for (int c=0;c!=3;c++) diff[1][c]=point[2][c]-point[0][c];
	cross[0]=diff[0][1]*diff[1][2]-diff[0][2]*diff[1][1];
	cross[1]=diff[0][2]*diff[1][0]-diff[0][0]*diff[1][2];
	cross[2]=diff[0][0]*diff[1][1]-diff[0][1]*diff[1][0];
	norm=sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
	normal[0]=cross[0]/norm;
	normal[1]=cross[1]/norm;
	normal[2]=cross[2]/norm;
}

void C_Mesh3D::getNormalOfTetrahedron(long tetrahedronnumber, int face, double * normal){
	double point[4][3];
	double diff[2][3];
	double cross[3];
	double norm;
	for (int p=0;p!=4;p++) this->getCoordinates(this->tetrahedronlist[tetrahedronnumber*4+p], point[p]); 
	if (face==0){
		diff[0][0]=point[1][0]-point[0][0];
		diff[0][1]=point[1][1]-point[0][1];
		diff[0][2]=point[1][2]-point[0][2];
		diff[1][0]=point[2][0]-point[0][0];
		diff[1][1]=point[2][1]-point[0][1];
		diff[1][2]=point[2][2]-point[0][2];
	}
	if (face==1){
		diff[0][0]=point[1][0]-point[0][0];
		diff[0][1]=point[1][1]-point[0][1];
		diff[0][2]=point[1][2]-point[0][2];
		diff[1][0]=point[3][0]-point[0][0];
		diff[1][1]=point[3][1]-point[0][1];
		diff[1][2]=point[3][2]-point[0][2];
	}
	if (face==2){
		diff[0][0]=point[2][0]-point[1][0];
		diff[0][1]=point[2][1]-point[1][1];
		diff[0][2]=point[2][2]-point[1][2];
		diff[1][0]=point[3][0]-point[1][0];
		diff[1][1]=point[3][1]-point[1][1];
		diff[1][2]=point[3][2]-point[1][2];
	}
	if (face==3){
		diff[0][0]=point[0][0]-point[2][0];
		diff[0][1]=point[0][1]-point[2][1];
		diff[0][2]=point[0][2]-point[2][2];
		diff[1][0]=point[3][0]-point[2][0];
		diff[1][1]=point[3][1]-point[2][1];
		diff[1][2]=point[3][2]-point[2][2];
	}
	cross[0]=diff[0][1]*diff[1][2]-diff[0][2]*diff[1][1];
	cross[1]=diff[0][2]*diff[1][0]-diff[0][0]*diff[1][2];
	cross[2]=diff[0][0]*diff[1][1]-diff[0][1]*diff[1][0];
	norm=sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
	normal[0]=cross[0]/norm;
	normal[1]=cross[1]/norm;
	normal[2]=cross[2]/norm;
}

bool C_Mesh3D::isPointInsideAnyMaterial(const C_Vector3D &p) const
{
	for(int t = 0; t < numberoftetrahedra; t++)
	{
		C_Tetrahedron tetrahedron;
		for(int p = 0; p != 4; p++)
			getCoordinates(tetrahedronlist[t*4+p], tetrahedron.Ns[p]);

		if( tetrahedron.encloses_point(p) )
			return true;
	}

	return false;
}

void C_Mesh3D::getCenterOfTriangle(long trianglenumber, double * center){
	double point[3][3];
	for (int p=0;p!=3;p++) this->getCoordinates(this->trianglelist[trianglenumber*3+p], point[p]);
	center[0]=(point[0][0]+point[1][0]+point[2][0])/3;
	center[1]=(point[0][1]+point[1][1]+point[2][1])/3;
	center[2]=(point[0][2]+point[1][2]+point[2][2])/3;
}

void C_Mesh3D::getCenterOfTetrahedron(long tetrahedronnumber, double * center){
	double point[4][3];
	for (int p=0;p!=4;p++) this->getCoordinates(this->tetrahedronlist[tetrahedronnumber*4+p], point[p]); 
	center[0]=(point[0][0]+point[1][0]+point[2][0]+point[3][0])/4;
	center[1]=(point[0][1]+point[1][1]+point[2][1]+point[3][1])/4;
	center[2]=(point[0][2]+point[1][2]+point[2][2]+point[3][2])/4;
}

C_Mesh3D::C_Mesh3D(){
	this->numberofpoints=0;
	this->numberofedges=0;
	this->numberoftriangles=0;
	this->numberoftetrahedra=0;
}


C_Mesh3D::~C_Mesh3D(){
	this->pointlist.clear(); 
	this->edgelist.clear();
	this->edgemarkerlist.clear();
	this->trianglelist.clear();
	this->trianglemarkerlist.clear();
	this->tetrahedronlist.clear();
	this->tetrahedronmarkerlist.clear();
	this->numberofpoints = 0;
	this->numberofedges = 0;
	this->numberoftriangles = 0;
	this->numberoftetrahedra = 0;
}
