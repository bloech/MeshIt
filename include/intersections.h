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

#ifndef _INTERSECTIONS_H_
#define _INTERSECTIONS_H_

#include "c_vector.h"

/********************************************************************************/
/*	Triangle/triangle intersection test routine, by Tomas Moller, 1997.			*/
/*	See article "A Fast Triangle-Triangle Intersection Test",					*/
/*	Journal of Graphics Tools, 2(2), 1997										*/
/*	updated: 2001-06-20 (added line of intersection)							*/
/*																				*/
/*	int tri_tri_intersect(double V0[3],double V1[3],double V2[3],				*/
/*						  double U0[3],double U1[3],double U2[3])				*/
/*	parameters: vertices of triangle 1: V0,V1,V2								*/
/*				vertices of triangle 2: U0,U1,U2								*/
/* result    :	returns 1 if the triangles intersect, otherwise 0				*/
/*																				*/
/*	Here is a version withouts divisions (a little faster)						*/
/*	int NoDivTriTriIsect(double V0[3],double V1[3],double V2[3],				*/
/*						 double U0[3],double U1[3],double U2[3]);				*/
/*																				*/
/*	This version computes the line of intersection as well						*/
/*	(if they are not coplanar):													*/
/*	int tri_tri_intersect_with_isectline(double V0[3],double V1[3],				*/
/*										 double V2[3], double U0[3],			*/
/*										 double U1[3], double U2[3],			*/
/*										 int *coplanar, double isectpt1[3],		*/
/*										 double isectpt2[3]);					*/
/* coplanar returns whether the tris are coplanar								*/
/* isectpt1, isectpt2 are the endpoints of the line of intersection				*/
/********************************************************************************/

// if USE_EPSILON_TEST is true then we do a check: 
// if |dv|<EPSILON then dv=0.0;
// else no check is done (which is less robust)

#define USE_EPSILON_TEST 1
#define EPSILON 1e-12

// some macros
// this edge to edge test is based on Franlin Antonio's gem:
// "Faster Line Segment Intersection", in Graphics Gems III, pp. 199-202

#define EDGE_EDGE_TEST(V0,U0,U1)                      \
  Bx=U0[i0]-U1[i0];                                   \
  By=U0[i1]-U1[i1];                                   \
  Cx=V0[i0]-U0[i0];                                   \
  Cy=V0[i1]-U0[i1];                                   \
  f=Ay*Bx-Ax*By;                                      \
  d=By*Cx-Bx*Cy;                                      \
  if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))  \
  {                                                   \
    e=Ax*Cy-Ay*Cx;                                    \
    if(f>0)                                           \
    {                                                 \
      if(e>=0 && e<=f) return 1;                      \
    }                                                 \
    else                                              \
    {                                                 \
      if(e<=0 && e>=f) return 1;                      \
    }                                                 \
  }                                

#define EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2) \
{                                              \
  double Ax,Ay,Bx,By,Cx,Cy,e,d,f;               \
  Ax=V1[i0]-V0[i0];                            \
  Ay=V1[i1]-V0[i1];                            \
  /* test edge U0,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U0,U1);                    \
  /* test edge U1,U2 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U1,U2);                    \
  /* test edge U2,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U2,U0);                    \
}

#define POINT_IN_TRI(V0,U0,U1,U2)           \
{                                           \
  double a,b,c,d0,d1,d2;                     \
  /* is T1 completly inside T2? */          \
  /* check if V0 is inside tri(U0,U1,U2) */ \
  a=U1[i1]-U0[i1];                          \
  b=-(U1[i0]-U0[i0]);                       \
  c=-a*U0[i0]-b*U0[i1];                     \
  d0=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U2[i1]-U1[i1];                          \
  b=-(U2[i0]-U1[i0]);                       \
  c=-a*U1[i0]-b*U1[i1];                     \
  d1=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U0[i1]-U2[i1];                          \
  b=-(U0[i0]-U2[i0]);                       \
  c=-a*U2[i0]-b*U2[i1];                     \
  d2=a*V0[i0]+b*V0[i1]+c;                   \
  if(d0*d1>0.0)                             \
  {                                         \
    if(d0*d2>0.0) return 1;                 \
  }                                         \
}

// sort so that a<=b
#define SORT2(a,b,smallest)       \
             if(a>b)       \
             {             \
               double c;    \
               c=a;        \
               a=b;        \
               b=c;        \
               smallest=1; \
             }             \
             else smallest=0;

inline void
isect2(C_Vector3D vtx0, C_Vector3D vtx1, C_Vector3D vtx2,
	double vv0, double vv1, double vv2,
	double d0, double d1, double d2,
	double *isect0, double *isect1,
	C_Vector3D *isectpoint0, C_Vector3D *isectpoint1)
{
	double tmp = d0 / (d0 - d1);
	C_Vector3D diff;
	*isect0 = vv0 + (vv1 - vv0)*tmp;
	diff = (vtx1 - vtx0);
	diff *= tmp;
	*isectpoint0 = diff + vtx0;
	tmp = d0 / (d0 - d2);
	*isect1 = vv0 + (vv2 - vv0)*tmp;
	diff = (vtx2 - vtx0);
	diff *= tmp;
	*isectpoint1 = vtx0 + diff;
}

inline int
compute_intervals_isectline(C_Triangle T, 
	double vv0, double vv1, double vv2,
	double d0, double d1, double d2,
	double d0d1, double d0d2, 
	double* isect0, double* isect1,
	C_Vector3D *isectpoint0, C_Vector3D *isectpoint1)
{
	if(d0d1>0.0f)
	{
		// here we know that d0d2<=0.0 that is d0, d1 are on the same side of p2, 
		// while d2 is on the other side of p2 or on p2 itself 
		isect2(*T.Ns[2], *T.Ns[0], *T.Ns[1], vv2, vv0, vv1, d2, d0, d1, isect0, isect1, isectpoint0, isectpoint1);
	}
	else if(d0d2>0.0f)
	{
		// here we know that d0d1<=0.0 that is do and d2 are on the same side of p2,
		// while d1 is on the other side of p2 or on p2 itself
		isect2(*T.Ns[1], *T.Ns[0], *T.Ns[2], vv1, vv0, vv2, d1, d0, d2, isect0, isect1, isectpoint0, isectpoint1);
	}
	else if (d1*d2>0.0f || d0!=0.0f)
	{
		// here we know that dod1<=0.0 or that do!=0.0 that is that d1 and d2 are on the same side of p2,
		// while d0 is on the other side or on p2 itself
		isect2(*T.Ns[0], *T.Ns[1], *T.Ns[2], vv0, vv1, vv2, d0, d1, d2, isect0, isect1, isectpoint0, isectpoint1);
	}
	else if(d1!=0.0f)
	{
		// here we know that d0d1<=0.0 that is do and d2 are on the same side of p2,
		// while d1 is on the other side of p2 or on p2 itself
		isect2(*T.Ns[1], *T.Ns[0], *T.Ns[2], vv1, vv0, vv2, d1, d0, d2, isect0, isect1, isectpoint0, isectpoint1);
	}
	else if(d2!=0.0f)
	{
		// here we know that d0d2<=0.0 that is d0, d1 are on the same side of p2, 
		// while d2 is on the other side of p2 or on p2 itself 
		isect2(*T.Ns[2], *T.Ns[0], *T.Ns[1], vv2, vv0, vv1, d2, d0, d1, isect0, isect1, isectpoint0, isectpoint1);
	}
	else
	{
		// triangles are coplanar
		return 1;
	}
	return 0;
}

int
coplanar_tri_tri(C_Vector3D N, C_Triangle T1, C_Triangle T2)
{
	double a[3];
	double V0[3],V1[3],V2[3];
	double U0[3],U1[3],U2[3];
	short i0,i1;
	// first project onto an axis-aligned plane, that maximizes the area of the triangles,
	//compute indices: i0,i1.
	a[0] = FABS(N.x());
	a[1] = FABS(N.y());
	a[2] = FABS(N.z());
	if(a[0]>a[1])
	{
		if(a[0]>a[2])
		{
			// a[0] is the greatest
			i0=1;
			i1=2;
		}
		else
		{
			// a[2] is the greatest
			i0=0;
			i1=1;
		}
	}
	else
	{
		// a[0]<=a[1]
		if(a[2]>a[1])
		{
			// a[2] is thegreatest
			i0=0;
			i1=1;
		}
		else
		{
			// a[1] is the greatest
			i0=0;
			i1=2;
		}
	}
	// test all edges of triangle T1 against the edges of triangle T2
	V0[0]=T1.Ns[0]->x();
	V0[1]=T1.Ns[0]->y();
	V0[2]=T1.Ns[0]->z();
	V1[0]=T1.Ns[1]->x();
	V1[1]=T1.Ns[1]->y();
	V1[2]=T1.Ns[1]->z();
	V2[0]=T1.Ns[2]->x();
	V2[1]=T1.Ns[2]->y();
	V2[2]=T1.Ns[2]->z();
	U0[0]=T2.Ns[0]->x();
	U0[1]=T2.Ns[0]->y();
	U0[2]=T2.Ns[0]->z();
	U1[0]=T2.Ns[1]->x();
	U1[1]=T2.Ns[1]->y();
	U1[2]=T2.Ns[1]->z();
	U2[0]=T2.Ns[2]->x();
	U2[1]=T2.Ns[2]->y();
	U2[2]=T2.Ns[2]->z();
	EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2);
	EDGE_AGAINST_TRI_EDGES(V1,V2,U0,U1,U2);
	EDGE_AGAINST_TRI_EDGES(V2,V0,U0,U1,U2);
	// finally, test if T1 is totally contained in T2 or vice versa 
	POINT_IN_TRI(V0,U0,U1,U2);
	POINT_IN_TRI(U0,V0,V1,V2);
	return 0;
}

int 
tri_tri_intersect_with_isectline(C_Triangle T1, C_Triangle T2, int* coplanar, C_Vector3D *isectpt1, C_Vector3D *isectpt2)
{
//	1.	compute plane equation (p1) of triangle T1=(V0,V1,V2) 
//		p1: N1.X+d1=0 
	C_Vector3D E1=(*T1.Ns[1]-*T1.Ns[0]);
	C_Vector3D E2=(*T1.Ns[2]-*T1.Ns[0]);
	C_Vector3D N1; 
	cross(E1, E2, &N1);
	double d1 = -dot(N1,*T1.Ns[0]);
//	2.a compute signed distance of triangle T2=(U0,U1,U2) to plane p1
	double du0=dot(N1,*T2.Ns[0])+d1; 
	double du1=dot(N1,*T2.Ns[1])+d1;
	double du2=dot(N1,*T2.Ns[2])+d1;
//	2.b coplanarity robustness check 
//	If all points of T2 are coplanar 
//	(they lie on the same side of p1 = all signed distances are non equal zero and of the same sign)
//	then reject the case as trivial since no intersection occurs
#if USE_EPSILON_TEST==1
	if (FABS(du0) < EPSILON)
		du0 = 0.0;
	if (FABS(du1) < EPSILON)
		du1 = 0.0;
	if (FABS(du2) < EPSILON)
		du2 = 0.0;
#endif
		double du0du1 = du0*du1;
		double du0du2 = du0*du2;
	if(du0du1>0.0f && du0du2>0.0f)
		return 0;
// 3 compute plane equation (p2)of triangle T2=(U0,U1,U2)
//	 p2: N2.X+d2=0
	E1 = (*T2.Ns[1] - *T2.Ns[0]);
	E2 = (*T2.Ns[2] - *T2.Ns[0]);
	C_Vector3D N2;
	cross(E1, E2, &N2);
	double d2 = -dot(N2, *T2.Ns[0]);
// 4.a Compute signed distance of triangle T1=(V0,V1,V2) to plane p2
	double dv0 = dot(N2, *T1.Ns[0]) + d2;
	double dv1 = dot(N2, *T1.Ns[1]) + d2;
	double dv2 = dot(N2, *T1.Ns[2]) + d2;
// 4.b	coplanarity robustness check
//		If all points of T1 are coplanar 
//		(they lie on the same side of p2 = all signed distances are non equal zero and of the same sign)
//		then reject the case as trivial since no intersection occurs
#if USE_EPSILON_TEST==1
	if (FABS(dv0) < EPSILON)
		dv0 = 0.0;
	if (FABS(dv1) < EPSILON)
		dv1 = 0.0;
	if (FABS(dv2) < EPSILON)
		dv2 = 0.0;
#endif
	double dv0dv1 = dv0*dv1;
	double dv0dv2 = dv0*dv2;
	if(dv0dv1>0.0f && dv0dv2>0.0f)
		return 0;
// 5.a compute direction vector (d) of intersection line (L=o+td) between T1 and T2 
	C_Vector3D *d = new C_Vector3D; 
	cross(N1, N2, d);
// 5.b compute simplified projection onto the largest component of d 
	double max = FABS(d->x());
	short index = 0;
	double b = FABS(d->y());
	double c = FABS(d->z());
	double vp0, vp1, vp2;
	double up0, up1, up2;
	if(b>max) max=b,index=1;
	if(c>max) max=c,index=2;
	if(index==0)
	{
		vp0=T1.Ns[0]->x();
		vp1=T1.Ns[1]->x();
		vp2=T1.Ns[2]->x();
		up0=T2.Ns[0]->x();
		up1=T2.Ns[1]->x();
		up2=T2.Ns[2]->x();
	}
	else if(index==1)
	{
		vp0=T1.Ns[0]->y();
		vp1=T1.Ns[1]->y();
		vp2=T1.Ns[2]->y();
		up0=T2.Ns[0]->y();
		up1=T2.Ns[1]->y();
		up2=T2.Ns[2]->y();
	}
	else if(index==2)
	{
		vp0=T1.Ns[0]->z();
		vp1=T1.Ns[1]->z();
		vp2=T1.Ns[2]->z();
		up0=T2.Ns[0]->z();
		up1=T2.Ns[1]->z();
		up2=T2.Ns[2]->z();
	}
// 6.a compute interval for triangle T1
	double isect1[2];
	C_Vector3D isectpointA1,isectpointA2;
	*coplanar = compute_intervals_isectline(T1, vp0, vp1, vp2, dv0, dv1, dv2, dv0dv1, dv0dv2, &isect1[0], &isect1[1], &isectpointA1, &isectpointA2);
	if(*coplanar) 
		return coplanar_tri_tri(N1,T1,T2);
// 6.a compute interval for triangle T2
	double isect2[2];
	C_Vector3D isectpointB1,isectpointB2;
	compute_intervals_isectline(T2, up0, up1, up2, du0, du1, du2, du0du1, du0du2, &isect2[0], &isect2[1], &isectpointB1, &isectpointB2);
	int smallest1, smallest2;
	SORT2(isect1[0], isect1[1], smallest1);
	SORT2(isect2[0], isect2[1], smallest2);
	if(isect1[1]<isect2[0] || isect2[1]<isect1[0])
		return 0;
// now we know that the two triangles T1 and T2 intersect
	if(isect2[0]<isect1[0])
	{
		if (smallest1 == 0)
			*isectpt1 = isectpointA1;
		else
			*isectpt1 = isectpointA2;
		if(isect2[1]<isect1[1])
		{
			if (smallest2 == 0)
				*isectpt2 = isectpointB2;
			else
				*isectpt2 = isectpointB1;
		}
		else
		{
			if (smallest1 == 0)
				*isectpt2 = isectpointA2;
			else
				*isectpt2 = isectpointA1;
		}
	}
	else
	{
		if (smallest2 == 0)
			*isectpt1 = isectpointB1;
		else
			*isectpt1 = isectpointB2;
		if(isect2[1]>isect1[1])
		{
			if (smallest1 == 0)
				*isectpt2 = isectpointA2;
			else
				*isectpt2 = isectpointA1;
		}
		else
		{
			if (smallest2 == 0)
				*isectpt2 = isectpointB2;
			else
				*isectpt2 = isectpointB1;
		}
	}
	return 1;
}

int
triangle_ray_intersection(C_Triangle TRI, C_Vector3D O, C_Vector3D D, C_Vector3D * isectpt)
{
	double det, inv_det, u, v;
	double t;
	//Find vectors for two edges sharing V1
	C_Vector3D e1 = (*TRI.Ns[1] - *TRI.Ns[0]);
	C_Vector3D e2 = (*TRI.Ns[2] - *TRI.Ns[0]);
	//Begin calculating determinant - also used to calculate u parameter
	C_Vector3D P;
	cross(D, e2, &P);
	//if determinant is near zero, ray lies in plane of triangle
	det = dot(e1, P);
	//NOT CULLING
	if (det > -EPSILON && det < EPSILON)
		return 0;
	inv_det = 1.f / det;
	//calculate distance from V1 to ray origin
	C_Vector3D T = (O - *TRI.Ns[0]);
	//Calculate u parameter and test bound
	u = dot(T, P) * inv_det;
	//The intersection lies outside of the triangle
	if (u < 0.f || u > 1.f)
		return 0;
	//Prepare to test v parameter
	C_Vector3D Q;
	cross(T, e1, &Q);
	//Calculate V parameter and test bound
	v = dot(D, Q) * inv_det;
	//The intersection lies outside of the triangle
	if (v < 0.f || u + v  > 1.f)
		return 0;
	t = dot(e2, Q) * inv_det;
	if (t < 0.f || t > 1.f)
		return 0;
	else
	{
		*isectpt = O + t*D;
		return 1;
	}
	// No hit, no win
	return 0;
}

#endif	// _INTERSECTIONS_H_
