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

#include <QPolygonF>
#include <QThreadStorage>
#include <QtOpenGL/QtOpenGL>

#include <sstream>

#include "core.h"
#include "glwidget.h"

C_Colors Colors;

QLineF PolygonSelection::lastLine() const
{
	return QLineF(vertices[vertices.length()-2], vertices.last());
}

void PolygonSelection::reset()
{
	vertices.clear();
}

bool PolygonSelection::add(const QPointF &p)
{
	if( ! isValid() )
		return false;

	vertices.append(p);
	return true;
}

bool PolygonSelection::update(const QPointF &p)
{
	if( isEmpty() )
		return false;

	vertices.back() = p;

	if( isClosed() )
		vertices.back() = vertices.first();

	return true;
}

bool PolygonSelection::isEmpty() const
{
	return vertices.isEmpty();
}

bool PolygonSelection::isValid() const
{
	if( isClosed() )
		return true;

	for(int i = 1; i < vertices.length(); ++i)
	{
		QLineF line(vertices[i-1], vertices[i]);
		QPointF p;
		if( lastLine().intersect(line, &p) == QLineF::BoundedIntersection
		        && p != lastLine().p1() ) {
			return false;
		}
	}

	return true;
}

bool PolygonSelection::isClosed() const
{
	if( vertices.length() < 4 )
		return false;

	const QPointF& p1 = vertices.last();
	const QPointF& p2 = vertices.first();

	return sqrt(pow((p1.x()-p2.x()), 2) + pow((p1.y()-p2.y()), 2)) < 5;
}

bool PolygonSelection::getPoint(QPointF &p) const
{
	QPolygonF polygon;
	for(int i = 1; i < vertices.length(); ++i)
		polygon << vertices[i];

	for(int i = 2; i < vertices.length(); ++i)
	{
		qreal x = (vertices[i].x() + vertices[i-1].x() + vertices[i-2].x()) / 3;
		qreal y = (vertices[i].y() + vertices[i-1].y() + vertices[i-2].y()) / 3;

		QPointF centroid(x, y);
		if( polygon.containsPoint(centroid, Qt::OddEvenFill) )
		{
			p = centroid;
			return true;
		}
	}
	return false;
}

void PolygonSelection::draw(const GLWidget *glwidget) const
{
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glClear(GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	glOrtho(0.f, glwidget->width(), glwidget->height(), 0.f, 0.f, 1.f);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glDisable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Colors.White);

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, Colors.MintTulip);
	glColor3f(0, 0, 0);

	glBegin(GL_LINE_STRIP);
	for(int i = 0; i < vertices.length() - 1; ++i)
		glVertex3d(vertices[i].x(), vertices[i].y(), 0.0);
	glEnd();

	const GLfloat *color = isValid() ? Colors.MintTulip : Colors.Red;
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);

	glBegin(GL_LINE_STRIP);
	if( vertices.length() > 1 ) {
		glVertex3d(lastLine().p1().x(), lastLine().p1().y(), 0.0);
		glVertex3d(lastLine().p2().x(), lastLine().p2().y(), 0.0);
	}
	glEnd();

	glEnable(GL_LIGHT0);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
}


GradientControl &GradientControl::getInstance()
{
	/* Provide a global instance per thread. */
	static QThreadStorage<GradientControl> storage;
	if( ! storage.hasLocalData() )
		storage.setLocalData(GradientControl());
	return storage.localData();
}

double GradientControl::gradient() const
{
	return _gradient;
}

double GradientControl::meshsize() const
{
	return _meshsize;
}

int GradientControl::npoints() const
{
	return _npoints;
}

const double* GradientControl::pointlist() const
{
	return _pointlist;
}

const double* GradientControl::refinesize() const
{
	return _refinesize;
}

void GradientControl::update(double gradient, double meshsize, int npoints,
                             const double *pointlist, const double *refinesize)
{
	_gradient = gradient;
	_meshsize = meshsize;
	_npoints = npoints;
	_pointlist = pointlist;
	_refinesize = refinesize;
}


typedef REAL *vertex;
#define ONETHIRD 0.333333333333333333333333333333333333333333333333333333333333

extern "C" {
int triunsuitable(vertex triorg, vertex tridest, vertex triapex, REAL area)
{
	const GradientControl& g = GradientControl::getInstance();

	REAL grad = g.gradient();
	REAL sq_grad = grad*grad;
	REAL sq_meshsize = g.meshsize() * g.meshsize();

	REAL dxoa = triorg[0] - triapex[0];
	REAL dyoa = triorg[1] - triapex[1];
	REAL dxda = tridest[0] - triapex[0];
	REAL dyda = tridest[1] - triapex[1];
	REAL dxod = triorg[0] - tridest[0];
	REAL dyod = triorg[1] - tridest[1];

	/* Find the squares of the lengths of the triangle's three edges. */
	REAL oalen = dxoa * dxoa + dyoa * dyoa;
	REAL dalen = dxda * dxda + dyda * dyda;
	REAL odlen = dxod * dxod + dyod * dyod;

	/* Find the square of the length of the mean edge. */
	REAL sq_meanlen = ONETHIRD*(oalen+dalen+odlen);
	if( sq_meanlen > sq_meshsize )
		return 1;

	for(int v = 0; v < g.npoints(); v++)
	{
		REAL sq_refinesize = g.refinesize()[v] * g.refinesize()[v];
		REAL cx = (triorg[0] + tridest[0] + triapex[0])*ONETHIRD;
		REAL cy = (triorg[1] + tridest[1] + triapex[1])*ONETHIRD;

		const double& px = g.pointlist()[v * 2];
		const double& py = g.pointlist()[v * 2 + 1];
		REAL sq_dist = (cx - px)*(cx - px) + (cy - py)*(cy - py);

		if( sq_dist < sq_grad * (sq_meshsize - sq_refinesize) )
		{
			if( sq_meanlen > (sq_dist / sq_grad + sq_refinesize) )
				return 1;
		}
	}
	return 0;
}
}

bool check_opengl_error(std::string& err_msg)
{
	GLenum err = glGetError();
	if( err == GL_NO_ERROR )
		return true;

	if( err == GL_OUT_OF_MEMORY )
	{
		err_msg = "OpenGL error: Out of memory!";
	}
	else
	{
		std::ostringstream ss;
		ss << err;
		err_msg = "OpenGL error: Error code " + ss.str() + "!";
	}

	return false;
}
