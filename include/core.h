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

#ifndef _CORE_H_
#define _CORE_H_

#include <QLineF>
#include <QList>
#include <QPointF>

enum RotationMode { MODEL = 0, OBJECT = 1, VIEWPORT = 2 };

enum class HighlightingMode { NONE = 0, COMPONENT = 1, ALL = 2, SINGLE = 3 };

enum SelectionTool { NOTOOL = 0, SINGLE, MULTI, ALL, BUCKET, INVERT_TOOL,
                     POLYGON };
enum SelectionMode { NOMODE = 0, MARK, UNMARK, HOLE, INVERT };

class SelectionState
{
public:
	SelectionTool tool;
	SelectionMode mode;

public:
	SelectionState() { reset(); }
	void reset() { tool = NOTOOL; mode = NOMODE; }

	bool operator==(const SelectionTool& tool) { return this->tool == tool; }
	bool operator==(const SelectionMode& mode) { return this->mode == mode; }

	bool operator!=(const SelectionTool& tool) { return !(*this == tool); }
	bool operator!=(const SelectionMode& mode) { return !(*this == mode); }
};

class GLWidget;

class PolygonSelection
{
private:
	QList<QPointF> vertices;

private:
	QLineF lastLine() const;

public:
	void reset();
	bool add(const QPointF& p);
	bool update(const QPointF& p);

	bool isEmpty() const;
	bool isValid() const;
	bool isClosed() const;

	bool getPoint(QPointF& p) const;

	void draw(const GLWidget *glwidget) const;
};

class GradientControl
{
private:
	double _gradient;
	double _meshsize;
	int _npoints;
	const double *_pointlist;
	const double *_refinesize;

public:
	static GradientControl& getInstance();

	double gradient() const;
	double meshsize() const;
	int npoints() const;
	const double *pointlist() const;
	const double *refinesize() const;

public:
	void update(double gradient, double meshsize, int npoints,
	            const double *pointlist, const double *refinesize);
};

bool check_opengl_error(std::string& err_msg);

#endif
