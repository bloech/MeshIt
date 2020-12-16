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
 
#include <math.h>
#include "c_vector.h"

void C_Vector3D::rotX(double sin, double cos){
	C_Vector3D tmp;
	tmp.setY(y()*cos-z()*sin);
	tmp.setZ(y()*sin+z()*cos);
	setY(tmp.y());
	setZ(tmp.z());
}

void C_Vector3D::rotZ(double sin, double cos){
	C_Vector3D tmp;
	tmp.setX(x()*cos-y()*sin);
	tmp.setY(x()*sin+y()*cos);
	setX(tmp.x());
	setY(tmp.y());
}

C_Vector3D C_Vector3D::normalized() const
{
    // Need some extra precision if the length is very small.
    double len = xp * xp +
                 yp * yp +
                 zp * zp;
    return *this / sqrt(len);
}

/*! \brief
    Returns the dot product of \a v1 and \a v2.
*/
double C_Vector3D::dotProduct(const C_Vector3D& v1, const C_Vector3D& v2)
{
    return v1.xp * v2.xp + v1.yp * v2.yp + v1.zp * v2.zp;
}

C_Vector3D C_Vector3D::crossProduct(const C_Vector3D& v1, const C_Vector3D& v2)
{
    return C_Vector3D(v1.yp * v2.zp - v1.zp * v2.yp,
					  v1.zp * v2.xp - v1.xp * v2.zp,
                      v1.xp * v2.yp - v1.yp * v2.xp);
}

/*! \brief
	Returns the projection \a x on a line segment defined by \a v1 and \a v2.
	If the projection is not on the line segment \a v1, \a v2, the null vector is returned.

	Let a line in three dimensions be specified by two points v1=(x_1,y_1,z_1) and v2=(x_2,y_2,z_2) lying on it, 
	so a vector along the line is given by
	
	(1) v=v1+t(v2-v1) 	
	
	The squared distance between a point on the line with parameter t and a point x=(x_0,y_0,z_0) is therefore
	
	(2) d^2=[(v1-x)+(v2-v1)t]^2
	
	To minimize the distance, set d(d^2)/dt=0 and solve for t to obtain
	
	(3) t=-((v1-x)dotProduct(v2-v1))/(|v2-v1|^2). 
	
	If \a t is in the range of 0 to 1 the projected point is calculated and returned.
	
	If \a t is not in the range of 0 to 1 the null vector is returned.
*/
C_Vector3D C_Vector3D::projectionTo(const C_Vector3D& v1, const C_Vector3D& v2)
{
	double t=-this->dotProduct((v1-*this),(v2-v1))/(v2-v1).lengthSquared();
	if (0<=t && t<1){
		return v1+t*(v2-v1);
	}
	return C_Vector3D(0,0,0);
}

C_Vector3D C_Vector3D::normal
        (const C_Vector3D& v1, const C_Vector3D& v2, const C_Vector3D& v3)
{
    return crossProduct((v2 - v1), (v3 - v1)).normalized();
}

double C_Vector3D::length() const
{
    return sqrt(xp * xp + yp * yp + zp * zp);
}

/*! \brief
    Returns the squared length of the vector from the origin.
    This is equivalent to the dot product of the vector with itself.

    \sa length(), dotProduct()
*/
double C_Vector3D::lengthSquared() const
{
    return xp * xp + yp * yp + zp * zp;
}

/*! \brief
    Returns the distance from this vertex to a plane defined by
    the vertex \a plane and a \a normal unit vector.  The \a normal
    parameter is assumed to have been normalized to a unit vector.

    The return value will be negative if the vertex is below the plane,
    or zero if it is on the plane.

    \sa normal(), distanceToLine()
*/
double C_Vector3D::distanceToPlane
        (const C_Vector3D& plane, const C_Vector3D& normal) const
{
    return dotProduct(*this - plane, normal);
}


/*! \brief 
	Returns the distance that this vertex is from a line defined
    by \a point and the unit vector \a direction.

    If \a direction is a null vector, then it does not define a line.
    In that case, the distance from \a point to this vertex is returned.

    \sa distanceToPlane()
*/
double C_Vector3D::distanceToLine
        (const C_Vector3D& point, const C_Vector3D& direction) const
{
    if (direction.isNull())
        return (*this - point).length();
    C_Vector3D p = point + dotProduct(*this - point, direction) * direction;
    return (*this - p).length();
}
