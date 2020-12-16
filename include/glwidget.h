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

#ifndef _GLWIDGET_H_
#define _GLWIDGET_H_

#include <QtOpenGL/QGLWidget>
#include <QtWidgets/QRubberBand>
#include "core.h"
#include "c_vector.h"
#include "geometry.h"

class GLWidget : public QGLWidget
{
	Q_OBJECT
public:
	GLWidget(QWidget *parent);
	~GLWidget();
	void drawObject();
	void setPerspec(bool);
	bool isPerspec(){return Perspec;}
	void resetView();
	void saveScreen(QString);
	void makeAxis();
	double currentDist() const;

	C_Model *Model;
	bool drawAxis;
	GLuint listAxis;

	RotationMode rotationMode;
	SelectionState selectionState;
	PolygonSelection polygonSelection;

	C_Vector3D prevCenter;

	double zoomSensitivity;

private:
	enum FillMode { NORMAL, GREEDY };

signals:
	void SelectionWasMade(unsigned char,unsigned char,unsigned char);
public:
	enum ViewportChange { Z_DIST, X_ROT, Z_ROT, X_TRANS, Y_TRANS };
	void changeViewport(ViewportChange change, double delta);
	void setRotationCenter(const C_Vector3D& center);
	void moveViewport(const C_Vector3D& position, double dist, double xRot,
	                  double zRot);
	void moveViewport(const C_Vector3D &position, double dist);
	void moveViewport(const C_Vector3D& position, double dist,
	                  const C_Vector3D& normalVector);
protected:
	void initializeGL();
	void paintGL();
	void resizeGL(int, int);
	void mousePressEvent(QMouseEvent * );
	void mouseMoveEvent(QMouseEvent * );
	void mouseReleaseEvent(QMouseEvent * );
	void qgluPerspective(GLdouble, GLdouble, GLdouble, GLdouble);
	void computePosition(double& xp, double& yp, double delta) const;
private slots:
private:
	unsigned char *selectionReadPixel(int, int, int, int);
	void selectionSearchAll(unsigned char *buffer, int width, int height);
	void selectionFloodFill(unsigned char *buffer, int x, int y,
	                        int width, int height, FillMode mode);
	void applySelection(int, int, int, int, int mx = 0, int my = 0);
	void addRGB(unsigned char, unsigned char, unsigned char);

	double xRot,yRot,zRot;
	bool Perspec;
	double xTrans,yTrans,zTrans;
	double left,right,top,bottom;
	double aspect,dist;
	QPoint lastGlobalPos,currentGlobalPos;
	QRect globalRect;
	QList <float> selectedRGBs;
	QRubberBand *rubberBand{ rubberBand = NULL };
};

#endif	// _GLWIDGET_H_
