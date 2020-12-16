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

#ifndef _C_VECTOR_H_
#define _C_VECTOR_H_

#include <QtGui/QtGui>
#include <math.h>

#define FABS(x) ((double)fabs(x))	// implement as is fastest on your machine

#define radians(x) ((x) * M_PI / 180.0)
#define degrees(x) ((x) * 180.0 / M_PI)
#define sign(x) ((x > 0) - (x < 0))

class C_Vector3D
{
public:
	C_Vector3D();
	C_Vector3D(double xpos, double ypos, double zpos);
	double x() const;
	double y() const;
	double z() const;
	QString type() const;
	void setX(double x);
	void setY(double y);
	void setZ(double z);
	void setType(QString type);

	C_Vector3D &operator+=(const C_Vector3D &vector);
	C_Vector3D &operator-=(const C_Vector3D &vector);
	C_Vector3D &operator*=(double factor);
	C_Vector3D &operator*=(const C_Vector3D &vector);
	C_Vector3D &operator/=(double divisor);

	friend inline const C_Vector3D operator+(const C_Vector3D &v1, const C_Vector3D &v2);
	friend inline const C_Vector3D operator-(const C_Vector3D &v1, const C_Vector3D &v2);
	friend inline const C_Vector3D operator/(const C_Vector3D &vector, double divisor);
	friend inline const C_Vector3D operator*(double factor, const C_Vector3D &vector);
	friend inline bool operator==(const C_Vector3D &v1, const C_Vector3D &v2);

	inline bool equals(const C_Vector3D &vector, double tolerance = 0.0) const;

	int intID;
	int triID;
	int tetID;
	int matID;

private:
	QString typep;
	double xp;
	double yp;
	double zp;
	double normp;
};

inline C_Vector3D::C_Vector3D() : xp(0), yp(0), zp(0), typep("DEFAULT")
{}
inline C_Vector3D::C_Vector3D(double xpos, double ypos, double zpos) : xp(xpos), yp(ypos), zp(zpos), typep("DEFAULT")
{}

inline double C_Vector3D::x() const { return double(xp); }
inline double C_Vector3D::y() const { return double(yp); }
inline double C_Vector3D::z() const { return double(zp); }
inline QString C_Vector3D::type() const { return QString(typep); }
inline void C_Vector3D::setX(double aX) { xp = aX; }
inline void C_Vector3D::setY(double aY) { yp = aY; }
inline void C_Vector3D::setZ(double aZ) { zp = aZ; }
inline void C_Vector3D::setType(QString aType) { typep = aType; }

inline C_Vector3D &C_Vector3D::operator+=(const C_Vector3D &vector)
{
	xp += vector.xp;
	yp += vector.yp;
	zp += vector.zp;
	return *this;
}

inline C_Vector3D &C_Vector3D::operator-=(const C_Vector3D &vector)
{
	xp -= vector.xp;
	yp -= vector.yp;
	zp -= vector.zp;
	return *this;
}

inline C_Vector3D &C_Vector3D::operator*=(double factor)
{
	xp *= factor;
	yp *= factor;
	zp *= factor;
	return *this;
}

inline C_Vector3D &C_Vector3D::operator*=(const C_Vector3D &vector)
{
	xp *= vector.xp;
	yp *= vector.yp;
	zp *= vector.zp;
	return *this;
}

inline C_Vector3D &C_Vector3D::operator/=(double divisor)
{
	xp /= divisor;
	yp /= divisor;
	zp /= divisor;
	return *this;
}

inline const C_Vector3D operator+(const C_Vector3D &v1, const C_Vector3D &v2)
{
	return C_Vector3D(v1.xp + v2.xp, v1.yp + v2.yp, v1.zp + v2.zp);
}

inline const C_Vector3D operator-(const C_Vector3D &v1, const C_Vector3D &v2)
{
	return C_Vector3D(v1.xp - v2.xp, v1.yp - v2.yp, v1.zp - v2.zp);
}

inline const C_Vector3D operator/(const C_Vector3D &vector, double divisor)
{
	return C_Vector3D(vector.xp / divisor, vector.yp / divisor, vector.zp / divisor);
}

inline const C_Vector3D operator*(double factor, const C_Vector3D &vector)
{
	return C_Vector3D(vector.xp * factor, vector.yp * factor, vector.zp * factor);
}

inline bool operator==(const C_Vector3D &v1, const C_Vector3D &v2)
{
	return v1.x() == v2.x() && v1.y() == v2.y() && v1.z() == v2.z();
}

bool C_Vector3D::equals(const C_Vector3D &vector, double tolerance) const
{
    return fabs(vector.x() - x()) < tolerance &&
           fabs(vector.y() - y()) < tolerance &&
           fabs(vector.z() - z()) < tolerance;
}

inline double length(C_Vector3D v1)
{
	return sqrt(v1.x() * v1.x() + v1.y() * v1.y() + v1.z() * v1.z());
}

inline double lengthSquared(C_Vector3D v1)
{
	return v1.x() * v1.x() + v1.y() * v1.y() + v1.z() * v1.z();
}

inline double dot(C_Vector3D v1, C_Vector3D v2)
{
	return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}

inline void cross(C_Vector3D v1, C_Vector3D v2, C_Vector3D *n)
{
	n->setX(v1.y()*v2.z() - v2.y()*v1.z());
	n->setY(v1.z()*v2.x() - v2.z()*v1.x());
	n->setZ(v1.x()*v2.y() - v2.x()*v1.y());
}

inline void normalize(C_Vector3D *n)
{
	double len = length(*n);
	n->setX(n->x()/len);
	n->setY(n->y() / len);
	n->setZ(n->z() / len);
}

inline void normal(C_Vector3D v1, C_Vector3D v2, C_Vector3D v3, C_Vector3D *n)
{
	cross((v2 - v1), (v3 - v1), n);
	normalize(n);
}

/*
Returns the projection v1 on a line segment defined by s1 and s2.
If the projection is not on the line segment [s1, s2] --> the null vector is returned.
Let a line in three dimensions be specified by two points s1=(x_1,y_1,z_1) and s2=(x_2,y_2,z_2) lying on it,
so a vector along the line is given by:
(1) n = s1 + t*(s2-s1)
The squared distance between a point on the line with parameter t and a point v1=(x_0,y_0,z_0) is therefore:
(2) d^2 = [(s1-v1) + (s2-s1)*t]^2
To minimize the distance, set d_(d^2)/dt=0 and solve for t to obtain:
(3) t = -((s1-v1)dotProduct(s2-s1))/(|s2-s1|^2).
If 
	t is in the range of 0 to 1 the projected point is calculated and returned.
else If 
	t is not in the range of 0 to 1 the null vector is returned.
*/
inline void projectTo(C_Vector3D v1, C_Vector3D s1, C_Vector3D s2, C_Vector3D *n)
{
	double t = -dot((s1 - v1), (s2 - s1)) / lengthSquared(s2 - s1);
	if (0 <= t && t<1)
	{
		n->setX(s1.x() + t*(s2.x() - s1.x()));
		n->setY(s1.y() + t*(s2.y() - s1.y()));
		n->setZ(s1.z() + t*(s2.z() - s1.z()));
		return;
	}
	n->setX(0);
	n->setY(0);
	n->setZ(0);
}

inline double getAngleOnX(const C_Vector3D& vec)
{
	double a = acos(vec.y());
	return sin(a) == vec.z() ? degrees(a) : degrees(-a);
}

inline double getAngleOnZ(const C_Vector3D& vec)
{
	double a = asin(vec.y());
	return cos(a) == vec.x() ? degrees(a) : degrees(-a);
}

inline C_Vector3D stripX(const C_Vector3D& v)
{
	C_Vector3D res(0, v.y(), v.z());
	normalize(&res);
	return res;
}

inline C_Vector3D stripZ(const C_Vector3D& v)
{
	C_Vector3D res(v.x(), v.y(), 0);
	normalize(&res);
	return res;
}

inline C_Vector3D rotateAroundZ(const C_Vector3D& vec, double angle)
{
	C_Vector3D res;
	res.setX(cos(radians(angle)) * vec.x() -
	         sin(radians(angle)) * vec.y());
	res.setY(sin(radians(angle)) * vec.x() +
	         cos(radians(angle)) * vec.y());
	res.setZ(vec.z());
	return res;
}

inline C_Vector3D rotateAroundX(const C_Vector3D& vec, double angle)
{
	C_Vector3D res;
	res.setX(vec.x());
	res.setY(cos(radians(angle)) * vec.y() -
	         sin(radians(angle)) * vec.z());
	res.setZ(sin(radians(angle)) * vec.y() +
	         cos(radians(angle)) * vec.z());
	return res;
}

#endif	// _C_VECTOR_H_
