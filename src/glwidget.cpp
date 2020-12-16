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

#include <QtGui/QtGui>
#include <QtOpenGL/QtOpenGL>
#include <math.h>
#include "core.h"
#include "glwidget.h"
#include <array>
#include <map>
#include <queue>

#ifdef _WIN32
#include <GL/glu.h>
#elif __LINUX__
#include <GL/glu.h>
#elif __APPLE__
#include <glu.h>
#endif

GLWidget::GLWidget(QWidget *parent) : QGLWidget(parent)
{
	this->xRot=90;
	this->yRot = this->zRot = 0;
	this->xTrans = this->yTrans = this->zTrans = 0;
	this->Perspec = true;
	this->drawAxis = true;
	this->rotationMode = RotationMode::MODEL;
	this->prevCenter = C_Vector3D(0, 0, 0);
	this->zoomSensitivity = 0.4;
}

GLWidget::~GLWidget()
{
	makeCurrent();
}

void
GLWidget::drawObject()
{
	glPushMatrix();

	/* Draw markers (e.g. missing selections). */
	glCallList(Model->listMarkers);

	/* Draw error markers. */
	glCallList(Model->listErrorMarkers);

	// draw surface objects
	for (int s = 0; s != this->Model->Surfaces.length(); s++)
	{
		if (Model->Surfaces[s].drawScatteredData)
			glCallList(Model->Surfaces[s].listScatteredData);
		if (Model->Surfaces[s].drawConvexHull)
			glCallList(Model->Surfaces[s].listConvexHull);
		if (Model->Surfaces[s].drawFaces)
			glCallList(Model->Surfaces[s].listFaces);
		if (Model->Surfaces[s].drawEdges)
			glCallList(Model->Surfaces[s].listEdges);
		if (Model->Surfaces[s].drawIntEdges)
			glCallList(Model->Surfaces[s].listIntEdges);
		if (Model->Surfaces[s].drawIntVertices)
			glCallList(Model->Surfaces[s].listIntVertices);
		if (Model->Surfaces[s].drawConstraints)
			glCallList(Model->Surfaces[s].listConstraints);
	}
	// draw polyline objects
	for (int p = 0; p != Model->Polylines.length(); p++)
	{
		if (Model->Polylines[p].drawScatteredData)
			glCallList(Model->Polylines[p].listScatteredData);
		if (Model->Polylines[p].drawEdges)
			glCallList(Model->Polylines[p].listEdges);
		if (Model->Polylines[p].drawVertices)
			glCallList(Model->Polylines[p].listVertices);
		if (Model->Polylines[p].drawIntVertices)
			glCallList(Model->Polylines[p].listIntVertices);
		if (Model->Polylines[p].drawConstraints)
			glCallList(Model->Polylines[p].listConstraints);
	}
	// draw materials
	if (Model->drawMats)
		glCallList(Model->listMats);
	// draw tets
	if (Model->drawTets)
		glCallList(Model->listTets);
	glPopMatrix();
}

void
GLWidget::setPerspec(bool pers)
{
	Perspec = pers;
	glDraw();
}

void
GLWidget::resetView()
{
	dist = 4;
	xRot = 90.0;
	yRot = zRot = xTrans = yTrans = zTrans = 0.0;
	glDraw();
}

double GLWidget::currentDist() const
{
	return dist;
}

void GLWidget::computePosition(double &xp, double &yp, double delta) const
{
	double alpha[2];
	double radius = sqrt(xp*xp + yp*yp);

	alpha[0] = acos(xp / radius) * 180.0 / M_PI;
	alpha[1] = -alpha[0];

	double a = 0;
	for(int i = 0; i < 2; ++i)
	{
		if( fabs(radius * sin(alpha[i] * M_PI / 180.0) - yp) < 1e-10 ) {
			a = alpha[i];
			break;
		}
	}

	a = a + delta;
	xp = radius * cos(a * M_PI / 180.0);
	yp = radius * sin(a * M_PI / 180.0);
}

void GLWidget::changeViewport(GLWidget::ViewportChange change, double delta)
{
	switch(change) {
	case Z_DIST:
		dist += delta;
		break;
	case X_TRANS:
		xTrans += delta;
		break;
	case Y_TRANS:
		yTrans += delta;
		break;
	case X_ROT:
		if( rotationMode == RotationMode::VIEWPORT )
			computePosition(yTrans, zTrans, delta);
		xRot -= delta;
		if( xRot > 360 )
			xRot -= 360;
		if( xRot < 0 )
			xRot += 360;
		break;
	case Z_ROT:
		if( rotationMode == RotationMode::VIEWPORT ) {
			computePosition(xTrans, zTrans, delta);
			if( xRot > 180 ) delta *= -1;
		}
		zRot -= delta;
		if( zRot > 360 )
			zRot -= 360;
		if( zRot < 0 )
			zRot += 360;
		break;
	}

	glDraw();
}

void GLWidget::moveViewport(const C_Vector3D& position, double dist,
                            double xRot, double zRot)
{
	this->xRot = xRot;
	this->zRot = zRot;

	RotationMode mode = rotationMode;
	rotationMode = RotationMode::MODEL;
	glDraw();

	xTrans = 0;
	yTrans = 0;
	zTrans = 0;
	setRotationCenter(position);
	this->dist = dist;
	glDraw();

	rotationMode = mode;
	glDraw();
}

void GLWidget::moveViewport(const C_Vector3D &position, double dist)
{
	return moveViewport(position, dist, xRot, zRot);
}

void GLWidget::moveViewport(const C_Vector3D &position, double dist,
                            const C_Vector3D &normalVector)
{
	const C_Vector3D& nvZ = stripZ(normalVector);
	double z_angle = getAngleOnZ(nvZ) + 270;

	const C_Vector3D& rotatedNormal = rotateAroundZ(normalVector, z_angle);

	const C_Vector3D& nvX = stripX(rotatedNormal);
	double x_angle = getAngleOnX(nvX) + 90;

	return moveViewport(position, dist, x_angle, z_angle);
}

void GLWidget::setRotationCenter(const C_Vector3D &center)
{
	C_Vector3D res = rotateAroundX(rotateAroundZ(center, zRot), -xRot);
	xTrans += res.x();
	yTrans += res.y();
	zTrans += res.z();
}

void
GLWidget::saveScreen(QString FileName)
{
	QFileInfo fileNameModel_info(FileName);
	QString Name = fileNameModel_info.absoluteFilePath();
	QImage image = this->grabFrameBuffer(false);
	image.save(Name);
}

void
GLWidget::makeAxis()
{
	glPushMatrix();
	glClear(GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	glTranslated(-aspect, -1, -3);
	glRotated(-xRot, 1.0, 0.0, 0.0);
	glRotated(zRot, 0.0, 0.0, 1.0);

	glDisable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Model->Cols.White);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Model->Cols.Grey);

	glBegin(GL_LINES);
		glVertex3d(0.0, 0.0, 0.0);
		glVertex3d(0.1, 0.0, 0.0);
		glVertex3d(0.0, 0.0, 0.0);
		glVertex3d(0.0, 0.1, 0.0);
		glVertex3d(0.0, 0.0, 0.0);
		glVertex3d(0.0, 0.0, 0.1);
	glEnd();

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Model->Cols.Yellow);
	Model->glWrite("x",0.2,0.0,0.0,0.02);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Model->Cols.Red);
	Model->glWrite("y",0.0,0.2,0.0,0.02);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Model->Cols.Green);
	Model->glWrite("z",0.0,0.0,0.2,0.02);

	glEnable(GL_LIGHT0);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Model->Cols.Grey);

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Model->Cols.BlueTrans);
	glBegin(GL_QUADS);
		glNormal3d(1.0,0.0,0.0);
		glVertex3d(0.00,0.00,0.00);
		glVertex3d(0.00,0.05,0.00);
		glVertex3d(0.00,0.05,0.05);
		glVertex3d(0.00,0.00,0.05);
		glNormal3d(1.0,0.0,0.0);
		glVertex3d(0.05,0.00,0.00);
		glVertex3d(0.05,0.05,0.00);
		glVertex3d(0.05,0.05,0.05);
		glVertex3d(0.05,0.00,0.05);
		glNormal3d(0.0,-1.0,0.0);
		glVertex3d(0.00,0.00,0.00);
		glVertex3d(0.05,0.00,0.00);
		glVertex3d(0.05,0.00,0.05);
		glVertex3d(0.00,0.00,0.05);
		glNormal3d(0.0,-1.0,0.0);
		glVertex3d(0.00,0.05,0.00);
		glVertex3d(0.05,0.05,0.00);
		glVertex3d(0.05,0.05,0.05);
		glVertex3d(0.00,0.05,0.05);
		glNormal3d(0.0,0.0,1.0);
		glVertex3d(0.00,0.00,0.00);
		glVertex3d(0.05,0.00,0.00);
		glVertex3d(0.05,0.05,0.00);
		glVertex3d(0.00,0.05,0.00);
		glNormal3d(0.0,0.0,1.0);
		glVertex3d(0.00,0.00,0.05);
		glVertex3d(0.05,0.00,0.05);
		glVertex3d(0.05,0.05,0.05);
		glVertex3d(0.00,0.05,0.05);
	glEnd();
	
	double incAngle = 30;
	double incPI = incAngle*0.017453293;
	double radius = 0.01;
	double nx,ny;
	double x1,y1;
	double x2,y2;

	glBegin(GL_TRIANGLES);
	for (double rad=0; rad<360/incAngle; rad++)
	{
		nx = sin(rad*incPI + 0.5*incPI);
		ny = cos(rad*incPI + 0.5*incPI);
		x1 = sin(rad*incPI)*radius;
		y1 = cos(rad*incPI)*radius;
		x2 = sin(rad*incPI + incPI)*radius;
		y2 = cos(rad*incPI + incPI)*radius;
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Model->Cols.GreenTrans);
		glNormal3d(-nx,-ny,0.1);
		glVertex3d(0.0,0.0,0.15);
		glVertex3d(x1,y1,0.1);
		glVertex3d(x2,y2,0.1);
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Model->Cols.RedTrans);
		glNormal3d(nx,0.1,ny);
		glVertex3d(0.0,0.15,0.0);
		glVertex3d(x1,0.1,y1);
		glVertex3d(x2,0.1,y2);
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,Model->Cols.YellowTrans);
		glNormal3d(-0.1,-nx,-ny);
		glVertex3d(0.15,0.0,0.0);
		glVertex3d(0.1,x1,y1);
		glVertex3d(0.1,x2,y2);
	}
	glEnd();
	glPopMatrix();
}

void
GLWidget::initializeGL()
{
	GLfloat lightPos[4] = { 0.0, 0.0, 1.0, 0.0 };
	glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_FLAT);

	glPointSize(4.0f);
	glLineWidth(2.0f);
	glLineStipple(2, 0xAAAA);

	glEnable(GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_NORMALIZE);
	glClearColor(0.318f, 0.341f, 0.431f, 1.0f);   

	Model = new C_Model;
}

void
GLWidget::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	if (Perspec)
		qgluPerspective(45.0, aspect, 0.01, 200.0);
	else
		glOrtho(-2 * aspect, 2 * aspect, -2, 2, -2000, 2000);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

	C_Vector3D center(0, 0, 0);
	if( rotationMode == RotationMode::OBJECT )
		center = Model->getGeomCenter();

	if( center.x() != prevCenter.x() ||
	    center.y() != prevCenter.y() ||
	    center.z() != prevCenter.z() )
	{
		C_Vector3D reset(-prevCenter.x(), -prevCenter.y(), -prevCenter.z());
		setRotationCenter(reset);

		if( rotationMode == RotationMode::OBJECT )
			setRotationCenter(center);

		prevCenter = center;
	}

	glTranslated(0.0, 0.0, -dist);
	glTranslated(xTrans, yTrans, zTrans);
	glRotated(-xRot, 1.0, 0.0, 0.0);
	glRotated(zRot, 0.0, 0.0, 1.0);
	glTranslated(-center.x(), -center.y(), -center.z());

	drawObject();
	if (this->drawAxis)
		this->makeAxis();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	polygonSelection.draw(this);
}

void
GLWidget::resizeGL(int width, int height)
{
	aspect = double(width) / double(height);
	glViewport(0, 0, width, height);
}

void
GLWidget::mousePressEvent(QMouseEvent *event)
{
	lastGlobalPos = event->globalPos();
	this->globalRect.setTopLeft(this->mapToGlobal(this->pos()));
	this->globalRect.setWidth(this->width());
	this->globalRect.setHeight(this->height());

	int x = event->globalX() - globalRect.x();
	int y = globalRect.y() + globalRect.height() - event->globalY();

	if (event->buttons() & Qt::LeftButton && selectionState == SINGLE)
	{
		int size = 5;
		applySelection(x - (size - 1) / 2, y - (size - 1) / 2, size, size);
	}
	if (event->buttons() & Qt::LeftButton && selectionState == MULTI)
	{
		if (!rubberBand)
			rubberBand = new QRubberBand(QRubberBand::Rectangle);
		rubberBand->setGeometry(event->globalX(), event->globalY(), 0, 0);
		rubberBand->show();
	}
	if (event->buttons() & Qt::LeftButton && selectionState == ALL)
		applySelection(0, 0, this->width(), this->height());

	if (event->buttons() & Qt::LeftButton && selectionState == BUCKET)
		applySelection(0, 0, this->width(), this->height(), x, y);

	if (event->buttons() & Qt::LeftButton && selectionState == POLYGON)
	{
		if( polygonSelection.isEmpty() )
		{
			setMouseTracking(true);
			polygonSelection.add(event->pos());
		}

		if( polygonSelection.add(event->pos()) )
		{
			if( polygonSelection.isClosed() )
			{
				setMouseTracking(false);

				QPointF p;
				if( polygonSelection.getPoint(p) ) {
					applySelection(0, 0, this->width(), this->height(),
					               p.x(), globalRect.height() - p.y());
				}
				polygonSelection.reset();

				glDraw();
			}
		}
	}

	if (event->buttons() & Qt::LeftButton && selectionState == INVERT)
		applySelection(0, 0, this->width(), this->height());
}

void
GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
	if (this->rubberBand && this->rubberBand->isVisible())
	{
		rubberBand->hide();
		int x = this->rubberBand->x() - this->globalRect.x();
		int y = this->globalRect.y() + this->globalRect.height() - (this->rubberBand->y() + this->rubberBand->height());
		this->applySelection(x, y, this->rubberBand->width(), this->rubberBand->height());
	}
}

void
GLWidget::mouseMoveEvent(QMouseEvent *event)
{
	bool enableRotation = selectionState != MULTI && selectionState != BUCKET
	        && selectionState != POLYGON;
	bool enableZoom = polygonSelection.isEmpty();
	bool enableShift = polygonSelection.isEmpty();

	/* Adjust sensitivity based on current zoom level. */
	double scale = std::exp((dist + 1) / 5.0) - 0.5;
	scale = std::min(std::max(scale * zoomSensitivity, 0.0), 10.0);

	currentGlobalPos = event->globalPos();
	if (event->buttons() & Qt::LeftButton && enableRotation)
	{
		double zd = -(double)(event->globalX() - this->lastGlobalPos.x());
		double xd = (double)(event->globalY() - this->lastGlobalPos.y());
		this->lastGlobalPos = event->globalPos();
		changeViewport(X_ROT, xd * scale);
		changeViewport(Z_ROT, zd * scale);
	}

	if (event->buttons() & Qt::RightButton && enableZoom)
	{
		/* Adjust sensitivity based on current zoom level. */
		double zoom_scale = std::exp((dist + 1) / 2.0) - 0.75;
		zoom_scale = std::min(std::max(zoom_scale * zoomSensitivity, 0.0), 10.0);

		double delta = -(double)(event->globalY() - this->lastGlobalPos.y()) / 100.;
		this->lastGlobalPos = event->globalPos();
		changeViewport(Z_DIST, delta * zoom_scale);
	}

	if (event->buttons() & Qt::MidButton && enableShift)
	{
		double xd = (double)(event->globalX() - this->lastGlobalPos.x()) / 400.;
		double yd = -(double)(event->globalY() - this->lastGlobalPos.y()) / 400.;
		this->lastGlobalPos = event->globalPos();
		changeViewport(X_TRANS, xd * scale);
		changeViewport(Y_TRANS, yd * scale);
	}

	if (event->buttons() & Qt::LeftButton && selectionState == MULTI)
	{
		if (event->globalX() < this->globalRect.left())
			this->currentGlobalPos.setX(this->globalRect.left());
		if (event->globalX() > this->globalRect.right())
			this->currentGlobalPos.setX(this->globalRect.right());
		if (event->globalY() < this->globalRect.top())
			this->currentGlobalPos.setY(this->globalRect.top());
		if (event->globalY() > this->globalRect.bottom())
			this->currentGlobalPos.setY(this->globalRect.bottom());
		rubberBand->setGeometry(QRect(this->lastGlobalPos, this->currentGlobalPos).normalized());
	}

	if (selectionState == POLYGON)
	{
		polygonSelection.update(event->pos());
		glDraw();
	}
}

void
GLWidget::qgluPerspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar)
{
	const GLdouble ymax = zNear * tan(fovy * M_PI / 360.0);
	const GLdouble ymin = -ymax;
	const GLdouble xmin = ymin * aspect;
	const GLdouble xmax = ymax * aspect;
	glFrustum(xmin, xmax, ymin, ymax, zNear, zFar);
}

void
GLWidget::addRGB(unsigned char R, unsigned char G, unsigned char B)
{
	if (R == 81 && G == 87 && B==110)
		return;
	for (int s = 0; s != this->selectedRGBs.length(); s = s + 3)
		if (this->selectedRGBs[s + 0] == R && this->selectedRGBs[s + 1] == G && this->selectedRGBs[s + 2] == B)
			return;
	selectedRGBs.append(R);
	selectedRGBs.append(G);
	selectedRGBs.append(B);
}

unsigned char* GLWidget::selectionReadPixel(int x, int y, int width, int height)
{
	int pixelRatio = devicePixelRatio();
	int bufferSize = pixelRatio * width * pixelRatio * height * 4;
	unsigned char *buffer = new unsigned char[bufferSize];

	/* Draw constraints in a special selection mode. */
	for (int s = 0; s != Model->Surfaces.length(); s++)
		Model->Surfaces[s].makeConstraints(true);

	glDisable(GL_LIGHTING);
	this->paintGL();

	glReadPixels(pixelRatio*x, pixelRatio*y, pixelRatio*width, pixelRatio*height,
	             GL_RGBA, GL_UNSIGNED_BYTE, buffer);

	for (int s = 0; s != Model->Surfaces.length(); s++)
		Model->Surfaces[s].makeConstraints(false);

	glEnable(GL_LIGHTING);
	this->updateGL();

	return buffer;
}

void GLWidget::selectionSearchAll(unsigned char *buffer, int width, int height)
{
	int pixelRatio = devicePixelRatio();
	for (int p = 0; p != pixelRatio*width*pixelRatio*height; p++)
		this->addRGB(buffer[4 * p + 0], buffer[4 * p + 1], buffer[4 * p + 2]);
}

void GLWidget::selectionFloodFill(unsigned char *buffer, int x, int y,
                                  int width, int height, FillMode mode = NORMAL)
{
#define IDX(p) 4 * ( (p).y() * width + (p).x() )

	typedef std::array<unsigned char, 3> RGB;
	std::map<RGB, int> screenColorsCount;
	std::map<RGB, int> floodCount;

	int pixelRatio = devicePixelRatio();

	x = x * pixelRatio;
	y = y * pixelRatio;
	width = width * pixelRatio;
	height = height * pixelRatio;

	for(int i = 0; i < height*width; i++)
	{
		RGB key = {buffer[4*i], buffer[4*i+1], buffer[4*i+2]};
		screenColorsCount[key]++;
	}

	std::queue<QPoint> Q;
	QPoint n(x, y);

	int xdir[4] = { -1, 1,  0, 0 };
	int ydir[4] = {  0, 0, -1, 1 };

	if( n.x() >= 0 && n.x() < width && n.y() >= 0 && n.y() < height )
	{
		buffer[IDX(n)] = 0;
		buffer[IDX(n) + 1] = 0;
		buffer[IDX(n) + 2] = 0;

		Q.push(n);
	}

	while( ! Q.empty() )
	{
		QPoint n = Q.front();
		Q.pop();

		for(int i = 0; i < 4; i++) {

			QPoint p(n.x() + xdir[i], n.y() + ydir[i]);

			if( p.x() < 0 || p.x() >= width || p.y() < 0 || p.y() >= height )
				continue;

			unsigned char *rgb = buffer + IDX(p);
			if( rgb[0] == 0 && rgb[1] == 0 && rgb[2] == 0 )
				continue;

			if( rgb[0] != 81 || rgb[1] != 87 || rgb[2] != 110 )
			{
				RGB key = {rgb[0], rgb[1], rgb[2]};
				floodCount[key]++;
				if( mode != GREEDY )
					continue;
			}

			rgb[0] = 0; rgb[1] = 0; rgb[2] = 0;
			Q.push(p);
		}
	}

	std::map<RGB, int>::const_iterator it;
	for(it = floodCount.begin(); it != floodCount.end(); ++it)
	{
		if( it->second == screenColorsCount[it->first] || it->second > 2 )
			this->addRGB(it->first[0], it->first[1], it->first[2]);
	}

#undef IDX
}

void
GLWidget::applySelection(int x, int y, int width, int height, int mx, int my)
{
	unsigned char *buffer = selectionReadPixel(x, y, width, height);

	selectedRGBs.clear();

	if( selectionState == BUCKET || selectionState == POLYGON )
	{
		FillMode mode = selectionState == BUCKET ? NORMAL : GREEDY;
		selectionFloodFill(buffer, mx, my, width, height, mode);
	}
	else {
		selectionSearchAll(buffer, width, height);
	}

	for (int s = 0; s != this->selectedRGBs.length(); s = s + 3)
		emit SelectionWasMade(this->selectedRGBs[s + 0],
		                      this->selectedRGBs[s + 1],
		                      this->selectedRGBs[s + 2]);

	delete[] buffer;
}

/*void GLWidget::mousePressEvent(QMouseEvent *event)
{
this->origin = event->globalPos();
this->GLWidgetGlobal.setX(event->globalX()-event->x());
this->GLWidgetGlobal.setWidth(this->geometry().width());
this->GLWidgetGlobal.setY(event->globalY() - event->y());
this->GLWidgetGlobal.setHeight(this->geometry().height());
if (!rubberBand)
rubberBand = new QRubberBand(QRubberBand::Rectangle);
rubberBand->setGeometry(QRect(origin, QSize()));
rubberBand->show();
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
this->current = event->globalPos();
if (event->x() < 0) current.setX(GLWidgetGlobal.x());
if (event->x() > GLWidgetGlobal.width()) current.setX(GLWidgetGlobal.x()+GLWidgetGlobal.width());
if (event->y() < 0) current.setY(GLWidgetGlobal.y());
if (event->y() > GLWidgetGlobal.height()) current.setY(GLWidgetGlobal.y() + GLWidgetGlobal.height());
rubberBand->setGeometry(QRect(origin, current).normalized());
}*/
