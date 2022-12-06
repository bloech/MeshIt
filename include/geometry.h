/* MeshIt - a 3D mesh generator for fractured reservoirs
 *
 * Copyright (C) 2020
 *
 * Mauro Cacace (GFZ, cacace@gfz-potsdam.de),
 * Guido Bl�cher (GFZ, bloech@gfz-potsdam.de),
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

#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include <QtCore/QtCore>
#include <QtGui/QtGui>
#include <QtOpenGL/QGLWidget>

#include <iostream>
#include <list>

#include "c_vector.h"
#include "feflow.h"
#include "tetgen.h"
//	To compile MeshIt (Visual Studio) without having Exodus libraries included uncomment the following definition
// #define NOEXODUS

#ifndef NOEXODUS
#include "exodus.h"
#include "exodusII.h"
#endif

extern "C"
{
#define REAL double
#define VOID void
#define ANSI_DECLARATORS
#include "triangle.h"
}

class C_Eigenvalue
{
public:
	void ComputeEigenvalue();
	void Tridiagonal();
	bool QLAlgorithm();
	void DecreasingSort();
	void GuaranteeRotation();

	double Element[3][3];
	double Diag[3];
	double Subd[3];
	bool IsRotation;
};

class C_Colors
{
public:
	C_Colors();

	GLfloat LightRed[4];
	GLfloat Red[4];
	GLfloat RedTrans[4];
	GLfloat LightGreen[4];
	GLfloat Green[4];
	GLfloat GreenTrans[4];
	GLfloat LightBlue[4];
	GLfloat Blue[4];
	GLfloat BlueTrans[4];
	GLfloat LightMagenta[4];
	GLfloat Magenta[4];
	GLfloat MagentaTrans[4];
	GLfloat Yellow[4];
	GLfloat YellowTrans[4];
	GLfloat Cyan[4];
	GLfloat CyanTrans[4];
	GLfloat White[4];
	GLfloat WhiteTrans[4];
	GLfloat Grey[4];
	GLfloat MintTulip[4];
	GLfloat Orange[4];
};

class C_Material
{
public:
	C_Material();

	QList<C_Vector3D> Locations;
	bool drawMatFaces;
	bool drawMatEdges;
};

class C_VTU
{
public:
	void write(QString);
	void read(QString);
	void clear();

	int NumberOfPoints, NumberOfCells;
	QList<C_Vector3D> Points;
	QList<int> connectivity;
	QList<int> offsets;
	QList<int> types;
	QList<int> matType;
	QList<int> matR;
	QList<int> matG;
	QList<int> matB;
	QList<int> Object1;
	QList<int> Object2;
	QList<int> pointData;
};

class C_Line
{
public:
	void AddPosition();
	void AddPoint(C_Vector3D);
	bool isInside(C_Vector3D);
	void calculate_min_max();
	void Invert();
	void EraseSpecialPoints();
	void MakeCornersSpecial();
	bool SortByType(QString);
	void CleanIdenticalPoints();
	void RefineByLength(double);
	bool IsIdenticallyWith(C_Line);
	void appendNonExistingSegment(C_Vector3D, C_Vector3D);
	void GenerateFirstSplineOfSegments(C_Line *);
	C_Line calculateSkewLineTransversal(const C_Vector3D &, const C_Vector3D &, const C_Vector3D &, const C_Vector3D &);
	C_Vector3D getPointAtPos(double);

	QList<C_Vector3D> Ns;
	QList<double> NsPos;
	QString Type;
	double size;
	int Object[2];
	float RGB[3];
	C_Vector3D min;
	C_Vector3D max;
};

class C_Segment
{
public:
	C_Vector3D *Ns[2];
	int marker;
	int polyline;
};

class C_Polyline
{
public:
	C_Polyline();
	void calculate_segments(bool);
	void calculate_Constraints();
	void calculate_min_max();
	void makeScatteredData();
	void makeEdges();
	void makeVertices();
	void makeIntVertices();
	void makeConstraints();
	void makeVTU_SD();
	void makeVTU_SEG();
	void makeVTU_CON();

	QList<C_Vector3D> SDs;
	C_Line Path;
	QString Name;
	QString Type;
	double size;
	QList<C_Line *> Intersections;
	QList<C_Line> Constraints;
	C_Vector3D min;
	C_Vector3D max;
	C_Colors Cols;
	bool drawScatteredData;
	bool drawEdges;
	bool drawVertices;
	bool drawIntVertices;
	bool drawConstraints;
	GLuint listScatteredData;
	GLuint listEdges;
	GLuint listVertices;
	GLuint listIntVertices;
	GLuint listConstraints;
	C_VTU VTU;
	int MaterialID;
};

class C_Triangle
{
public:
	void calculate_min_max();
	void setNormalVector();

	C_Vector3D *Ns[3];
	C_Vector3D normal_vector;
	C_Vector3D min;
	C_Vector3D max;
	int marker;
	int mesh;
};

class C_Tetrahedron
{
public:
	C_Vector3D Ns[4];
	C_Triangle Tr[4];
	C_Vector3D Center;
	int region;

	bool encloses_point(const C_Vector3D &p) const;

private:
	bool is_on_same_side(const C_Vector3D &v1, const C_Vector3D &v2,
											 const C_Vector3D &v3, const C_Vector3D &v4,
											 const C_Vector3D &p) const;
};

class C_Box
{
public:
	C_Vector3D min, max, center;
	QList<const C_Triangle *> T1s;
	QList<const C_Triangle *> T2s;
	QList<const C_Vector3D *> N1s;
	QList<const C_Vector3D *> N2s;
	C_Line IntSegments;
	bool tri_in_box(const C_Triangle *Tri) const;
	bool seg_in_box(const C_Vector3D *V1, const C_Vector3D *V2) const;
	bool too_much_tri() const;
	bool too_much_seg() const;
	void split_tri(C_Line *IntSegments);
	void split_seg(QList<C_Vector3D *> *TPs, int I1, int I2);
	void calculate_center();
	void generate_subboxes();

private:
	C_Box *Box[8];
};

class C_Mesh3D
{
public:
	long numberofpoints, numberofedges, numberoftriangles, numberoftetrahedra;
	long numberofpointswithmaterial, numberofedgeswithmaterial, numberoftriangleswithmaterial, numberoftetrahedrawithmaterial;
	void getCoordinates(long pointnumber, C_Vector3D &point) const;
	void getCoordinates(long pointNumber, double *point) const;
	void getCenterOfTriangle(long trianglenumber, double *center);
	void getCenterOfTetrahedron(long tetrahedronnumber, double *center);
	void getNormalOfTriangle(long trianglenumber, double *normal);
	void getNormalOfTetrahedron(long tetrahedronnumber, int face, double *normal);
	bool isPointInsideAnyMaterial(const C_Vector3D &) const;
	QList<double> pointlist;
	QList<long> edgelist;
	QList<int> edgemarkerlist;
	QList<long> trianglelist;
	QList<int> trianglemarkerlist;
	QList<long> tetrahedronlist;
	QList<int> tetrahedronmarkerlist;

	// temporary generated for export
	long numberofpoints_export, numberofedges_export, numberoftriangles_export, numberoftetrahedra_export;
	QList<double> pointlist_export;
	QList<long> edgelist_export;
	QList<int> edgemarkerlist_export;
	QList<long> trianglelist_export;
	QList<int> trianglemarkerlist_export;
	QList<long> tetrahedronlist_export;
	QList<int> tetrahedronmarkerlist_export;

	C_Mesh3D();
	~C_Mesh3D();
};

/*! \class C_Surface
 *	\ingroup geometry
 *	\brief Definition file for Class C_Surface
 *	\details Class C_Surface is the template class for storing and handling all types of geometric objects defining a surface (FAULT, UNIT or BORDER) of the model.\n
 *	It inherits properties from classes C_Vector3D, C_Line and C_Triangle.\n
 *	In addition it provides all basic operations to be performed on a 2D surface.\n
 *	Those comprises:
 *	\arg C_Surface::calculate_Constraints();
 *	\arg C_Surface::separate_Constraints();
 *	\arg C_Surface::calculate_normal_vector();
 *	\arg C_Surface::calculate_min_max();
 *	\arg C_Surface::calculate_convex_hull();
 *	\arg C_Surface::interpolation();
 *	\arg C_Surface::alignIntersectionsToConvexHull();
 *	\arg C_Surface::calculate_triangles_coarse();
 *	\arg C_Surface::calculate_triangles_fine();
 *	\arg C_Surface::rotate_horizontal();
 *	\arg C_Surface::rotate_back();
 *	\arg C_Surface::PointInPolygon();
 */
class C_Surface
{
private:
	/*! \brief Double array to perform planar rotation around cartesian axes.
	 *	\details [0] => rotation is around the Z-axis;\n
	 *	[1] => rotation is around the X-axis;
	 */
	double sinp[2], cosp[2]; //[0] rotation around z and [1] rotation around x
public:
	QList<QList<double>> MatrixSpline;
	QList<double> VectorSpline;
	QList<C_Vector3D *> selectedSDs;
	void fillSpline();
	void SolveSpline();
	double *MatrixKriging;
	double *VectorKriging;
	double *KrigingWeights;
	int *KrigingPivot;
	double KrigingBeta;
	void fillKriging();
	void Crout_LU_Decomposition_with_Pivoting(double *A, int *pivot);
	void Crout_LU_with_Pivoting_Solve(double *LU, double *B, double *x, int *pivot);
	/// \brief Instance of Class C_VTU
	C_VTU VTU;
	void makeVTU_SD();	// VTU for Scattered Data Points
	void makeVTU_CH();	// VTU for Convex Hull
	void makeVTU_TRI(); // VTU for Triangles
	void makeVTU_CON(); // VTU for Segments and Holes
	/// \brief String identifing the Name of the 2D surface (instance of Class C_Surface) of the model.
	QString Name;
	/*! \brief String defining the Type of an instance of Class C_Surface.
	 *	\details Admissible Types are:
	 *	\arg UNIT => the 2D surface defines a geoloigcal unit of the model.
	 *	\arg FAULT => the 2D surface defines a fault of the model.
	 *	\arg BORDER => the 2D surface constrains a lateral boundary of the model.
	 */
	QString Type;
	/// \brief Double defines the maximum possible edge length for trinagles.
	double size;
	/// \brief Instance of Class C_Vector3D which defines the minimum extension of the model window.
	C_Vector3D min;
	/// \brief Instance of Class C_Vector3D which defines the maximum extension of the model window.
	C_Vector3D max;
	/// \brief Instance of Class C_vector3D which defines the normal vector of a two dimensional geometric object.
	C_Vector3D normal_vector;
	/// \brief List of instances of Class C_Vector3D defining the Scattered Input Data.
	QList<C_Vector3D> SDs;
	/// \brief List of instances of Class C_Vector3D defining the vertices of the triangles of the 2D surfaces.
	QList<C_Vector3D> Ns;
	int duplicates;
	/// \brief List of instances of Class C_Triangle defining the triangles for all  2D surfaces.
	QList<C_Triangle> Ts;
	/// \brief instance of Class C_Line to represent a Convex Hull.
	C_Line ConvexHull;
	/// \brief List of instances of Class C_Line to store the intersection polylines between two or more surfaces (reference to Class Model::Intersections).
	QList<C_Line *> Intersections;
	/// \brief List of instances of Class C_Line to store the segments used to constrain the final 3D tetrahedralization.
	QList<C_Line> Constraints;
	/// \brief List of instances of Class C_Vector3D to store the coordinates of points lying inside closed C_Line definining a hole for a surface in the 2D triangulation.
	QList<C_Vector3D> HoleCoords;
	/*!	\ingroup ImportRead
	 *	\brief Calculate the minimum and maximum extensions to set the model window.
	 *	\details MAX is given by a C_Vector3D having as (x,y,z)-coodinates the absolute maximum between all coordinates among all points defining all surfaces of the model.\n
	 *	MIN is given by a C_Vector3D having as (x,y,z)-coodinates the absolute minimum between all coordinates among all points defining all surfaces of the model.
	 */
	void calculate_min_max();
	/*! \ingroup PreMesh
	 *	\brief Calculate the best fitting plane based on a set of C_Vector3D Points.
	 *	\details The plane is described by three point and a normal vector. This is the "ordinary least squares" fit.
	 *	\note This is the "ordinary least squares" fit which is appropriate only when the z is expected to be a linear function of x and y.\n
	 *	A more general form for a "best fit plane" in 3-spaceis given by a "geometric" least squares.
	 *	However, even the latter will fail if your points are in a plane parallel to z-direction. For this case the z-component of the normal vector
	 *	is 0 and the x,y- components are calculated by use of the x and y values of the scattered data points.
	 */
	void calculate_normal_vector();
	/*! \ingroup PreMesh
		\brief Perform a rotation of a set of 3D scattered data points onto a two dimensional plane by a direction vector of the projection.
		\details Let M be the vector normal to the best fitting plane \sa void calculate_normal_vector(); and N be the vector normal to the plane you want to rotate into

			Calculate the rotation angle as

		costheta = dot(M,N)/(normalize(M)*normalize(N))

		Calculate the rotation axis as

		cross(M, N, axis)
		normalize (axis)

			where cross is a function that performs the cross product, and normalizes is the fuction building a unit vector.

			Compute the rotation matrix from the axis and angle as

			c = costheta
			s = sqrt(1-c*c)
			C = 1-c

		r1 = [ x*x*C+c    x*y*C-z*s  x*z*C+y*s ]
		r2 = [ y*x*C+z*s  y*y*C+c    y*z*C-x*s ]
			r3 = [ z*x*C-y*s  z*y*C+x*s  z*z*C+c   ]

			For each point (Ns, SDs, Convexhull, Intersections, Constaints), compute its corresponding point on the new plane as

		newpoint.x() = dot(r1, oldpoint)
		newpoint.y() = dot(r2, oldpoint)
		newpoint.z() = dot(r3, oldpoint)
	*/
	void rotate(bool onto_z);
	/*!	\ingroup PreMesh
	 *	\brief Determine the Convex Hull (two dimensional) of the scattered data points for each surface of the model.
	 *	\details It determines the Convex Hull of the clouds of points projected on their Best Fitting Plane (BFP),
	 *	\sa C_Surface::calculate_normal_vector() and C_Surface::rotate_horizontal().\n
	 *	In a second stage, it marks the corner points of the Convex Hull by type (\sa C_Vector3D::setType).
	 *	These points defines the constrained parts of the Convex Hull.\n
	 *	As a last phase, it performs a refinement of the parts of the Convex Hull by inserting additional points separated by a user
	 *	defined length parameter (\a Size).
	 *	\sa C_Line::EraseSpecialPoints() , C_Line::MakeCornersSpecial() , C_Line::SortByType(),
	 *	C_Line::AddPosition() , C_Line::RefineByLength().
	 */
	void calculate_convex_hull();
	/*!	\ingroup PreMesh
	 *	\brief Perform an IDW interpolation of the points describing the Convex Hulls to the scattered data input points.
	 */
	void interpolation(QString object /**< [in] String defined the type of object to interpolate*/, QString method);
	/*!\ingroup PreMesh
	 *	\brief Merges the first and the last point of intersections to belonging convexhull.
	 *	\details If the first or the last point of the intersection is closer than 1e-12 to one special point of the convexhull,
	 *	then the point gets the coordinates of this convexhull point.\n
	 *	If the first or last point of the intersection is close to a segment of the convexhull (smaller than 1e-12) and not close to one special point of the convex hull,
	 *	then this point is added to the convex hull as a special point.
	 */
	void alignIntersectionsToConvexHull();
	/*!	\ingroup Simplify
	 *	\brief Fill a list of C_Line (C_Surface::SegmentAndHoles) used to define the segments to constain the final 3d tetrahedralization.
	 *	\details It loops through the Convex Hull and Intersections points for each surface identifying the finite parts of the C_Line (defined as the parts between two non Default points)
	 *	and filling the pixel color array (C_Line::RGB) with increasing values for each part.
	 *	At the end all parts are characterized by having a unique value of their C_Line::RGB array.
	 */
	void calculate_Constraints();
	/*!	\ingroup Simplify
	 *	\brief Determine the coordinates of an interior point in a hole used as a constrain for the 2D triangulation.
	 *	\details The C_line defining a hole must be close.
	 *	In case the Hole is defined by intersection of different segments, it first connect all the segments in a closed line.\n
	 *	It then calculate the coordinates of a generic point lying inside the closed line which will be used as a reference for constraining the hole during the 2D triangulation.
	 */
	void separate_Constraints();
	/*!	\ingroup Mesh
	 *	\brief It triangulates
	 */
	void calculate_triangles(bool withConstraints, double gradient = 2.0);
	double sin(int) const;
	double cos(int) const;
	double IDW(double x, double y);
	/*!	\ingroup PreMesh
	 *	\brief return an interpolate value based on a SPLINE algorithm.
	 *	\sa C_Surface::interpolation()
	 */
	double SPLINE(double x, double y);
	/*!	\ingroup PreMesh
	 *	\brief return an interpolate value based on a KRIGING algorithm.
	 *	\sa C_Surface::interpolation()
	 */
	double KRIGING(double x, double y);
	void setSin(int, double);
	void setCos(int, double);
	C_Colors Cols;
	void clearScatteredData();
	bool drawScatteredData;
	void makeScatteredData();
	GLuint listScatteredData;
	bool drawConvexHull;
	void makeConvexHull();
	GLuint listConvexHull;
	bool drawFaces;
	void makeFaces();
	GLuint listFaces;
	bool drawEdges;
	void makeEdges();
	GLuint listEdges;
	bool drawIntEdges;
	void makeIntEdges();
	GLuint listIntEdges;
	bool drawIntVertices;
	void makeIntVertices();
	GLuint listIntVertices;
	bool drawConstraints;
	void makeConstraints(bool selectionMode = false);
	GLuint listConstraints;
	bool drawMatFaces;
	bool drawMatEdges;
	// bool isMaterial;
	int MaterialID;
	C_Surface();
	~C_Surface();
};

inline double C_Surface::sin(int axis) const { return double(sinp[axis]); }
inline double C_Surface::cos(int axis) const { return double(cosp[axis]); }

inline void C_Surface::setSin(int axis, double aSin) { sinp[axis] = aSin; }
inline void C_Surface::setCos(int axis, double aCos) { cosp[axis] = aCos; }

typedef std::pair<C_Triangle, QSet<const C_Surface *>> SelfIntersection;

/*! \class C_Model
 *	\ingroup geometry
 *
 */
class C_Model : public QObject
{
	Q_OBJECT
signals:
	void ModelInfoChanged();
	void PrintError(QString);
	void ErrorInfoChanged(QString);

private:
	void analyze_self_intersections(const tetgenio &out);

public:
	void makeVTU_INT();
	void makeVTU_MAT();
	void makeVTU_TET();
	C_VTU VTU;
	void VTU_SD_to_Polyline(QString name, QString type, int mat, double size);
	void VTU_SD_to_Surface(QString name, QString type, int mat, double size);
	void VTU_CH_to_Surface(QString name, QString type);
	void VTU_SEG_to_Polyline(QString name, QString type);
	void VTU_TRI_to_Surface(QString name, QString type);
	void VTU_CON_to_Polyline(QString name, QString type);
	void VTU_CON_to_Surface(QString name, QString type);
	void VTU_INT_to_Model();
	void VTU_MAT_to_Model();
	void VTU_TET_to_Model();
	C_Vector3D min, max;
	QList<C_Vector3D> Ns;						// Vertices of Tets
	QList<C_Triangle> Fs;						// All faces of final mesh
	QList<C_Triangle> FsWithMat;		// faces with material of final mesh
	QList<C_Segment> Ss;						// All segments of final mesh
	QList<C_Segment> SsWithMat;			// faces with material of final mesh
	QList<C_Tetrahedron> Ts;				// Tets of final mesh
	QList<C_Tetrahedron> TsWithMat; // Tets of final mesh
	C_Mesh3D *Mesh;
	void ExportLists_make();
	void ExportLists_clean();
	// C_Mesh3D * Mesh;
	int getMaterial(int vtkType, long listID);
	QList<C_Vector3D *> TPs;		 // Vertices of TRIPLEPOINTS
	QList<C_Surface> Surfaces;	 // faults and units
	QList<C_Polyline> Polylines; // wells spline
	QList<C_Line> Intersections; // all intersection splines
	bool PreTestIntersectionsSegmentTriangle(const C_Vector3D *S1, const C_Vector3D *S2, const C_Triangle *T) const;
	bool PreTestIntersectionsPolylineSurface(const C_Polyline *P, const C_Surface *S) const;
	bool PreTestIntersectionsSegmentSurface(const C_Vector3D *S1, const C_Vector3D *S2, const C_Surface *S) const;
	bool PreTestIntersectionsTrianglePolyline(const C_Triangle *T, const C_Polyline *P) const;
	bool PointOnTriangle(const C_Vector3D &, const C_Triangle &) const;
	/*! \ingroup PreMesh
	 *	\brief Merging of all convexhull points that share the same MergeID.
	 *	\details Based on a common MergeID, it searches for all the points within all the convexhulls sharing the same MergeID and merge all of them
	 *	in one point in the space by taking their summed averages.
	 *
	 *	It marks the new points as special of type 4.
	 *
	 *	It then clear the convexhulls by their duplicated points and then perform a refinement by means of a reference given length (\a minDistance) of the convexhull taking all special points as fixed points.
	 */
	void calculate_int_polyline(int s1, int s2);
	void calculate_int_point(int p, int s);
	void calculate_int_triplepoints(int I1, int I2);
	void insert_int_triplepoints();
	void calculate_tets(QString switches);
	void get_constraints(std::list<C_Line *> &res);
	bool has_selected_constraints();
	void set_all_constraints(const QString &type);
	void select_all_constraints();
	void deselect_all_constraints();
	void material_selections();
	bool verify_materials();
	QString FileNameModel;
	QString FileNameTmp;
	QStringList FileNamesTmp;
	QString FilePath;
	QString intAlgorythm;
	double preMeshGradient;
	double meshGradient;
	double ExportRotationAngle;
	void Open();
	void Save();
	void ReadGocadFile();
	void AddSurface(QString type);
	void AddPolyline(QString type);
	void DeleteSurface(QString surfaceName);
	void ExportFeFlow();
	void ExportOGS();
	void ExportCOMSOL();
	void ExportABAQUS(QString borderIDs);
	void ExportVTU3D();
	void ExportTIN(QString surfaceID);
	void ExportVTU2D(QString surfaceID);
#ifndef NOEXODUS
	void ExportEXODUS(QString borderIDs);
	void EXODUS_sides(C_Exodus *, QString);
	void EXODUS_element(C_Exodus *);
	void EXODUS_nodes(C_Exodus *);
#endif
	void calculate_min_max();
	void calculate_size_of_intersections();
	void calculate_size_of_constraints();
	void tranformForward();
	void tranformBackward();
	C_Vector3D shift;
	double scale;
	void calculateNumberWithMaterials();

	QList<C_Material> Mats;
	C_Colors Cols;
	bool drawTets;
	void makeTets(bool xCutEnable, double xCutValue, bool xDirection, bool yCutEnable, double yCutValue, bool yDirection, bool zCutEnable, double zCutValue, bool zDirection);
	GLuint listTets;
	void glWrite(QString string, double x, double y, double z, double scale);
	bool drawMats;
	void makeMats(int Material, int Location);
	GLuint listMats;
	void makeMarkers(const C_Line *selectedLine = NULL,
									 const QString &comp = "");
	GLuint listMarkers;
	void makeErrorMarkers(const C_Line &line);
	GLuint listErrorMarkers;

	C_Vector3D getGeomCenter();
	C_Polyline *findPolyline(const QString &name);
	C_Surface *findSurface(const QString &name);

	QList<SelfIntersection> selfIntersections;

	C_Model();
	~C_Model();
};

<<<<<<< HEAD
class C_PLC : public QObject
{
	Q_OBJECT

signals:
	void PrintError(QString);
=======
class C_PLC: public QObject
{
    Q_OBJECT

signals:
    void PrintError(QString);
>>>>>>> 0f66d14e58f8b34f68dd8559c728b92ccddf2bb4

public:
	/// @brief Name of the file to import
	QString fileName;
	/// @brief Name of the path where the file is stores
	QString filePath;
	/// @brief List of instances of Class C_Vector3D defining the vertices of the triangles of making the PLC
	QList<C_Vector3D> Ns;
	/// @brief List of instances of Class C_Triangle defining the triangles making the PLC
	QList<C_Triangle> Ts;
	/// @brief Pointer to an instance of the class C_Mesh3D to store and output the final 3D mesh
	C_Mesh3D *Mesh;
	/// @brief basic constructor and destructor
	C_PLC() {}
	~C_PLC() {}
	///	@brie read the PLC file
	void readPLCFile();
	/// @brief modify the options to give to tetgen for the final tetrahedralization();();
	void modifyTetgenOptions();
	/// @brief Perform the final 3D mesh on the PLC
	/// @param switches tetgen input functions
	void meshPLC(QString switches);
};

#endif // _GEOMETRY_H_
