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

#include <QtGui/QtGui>
#include "geometry.h"
#include "mainwindow.h"
#include "glwidget.h"

C_Model Model;
C_PLC PLC;

QDateTime startdate, enddate;

/********** Class AnimatedButton **********/
AnimatedButton::AnimatedButton(const QString &text, const QString &icon,
															 QWidget *parent) : QPushButton(text, parent)
{
	animation = new QMovie(icon);
	connect(animation, SIGNAL(frameChanged(int)),
					this, SLOT(updateAnimation(int)));

	connect(this, SIGNAL(clicked(bool)), this, SLOT(startAnimation()));
}

AnimatedButton::~AnimatedButton()
{
	if (animation)
		delete animation;
}

void AnimatedButton::updateAnimation(int)
{
	setIcon(animation->currentPixmap());
	setText("  " + text().trimmed());
}

void AnimatedButton::startAnimation()
{
	setEnabled(false);
	animation->start();
}

void AnimatedButton::stopAnimation()
{
	animation->stop();
	setIcon(QIcon());
	setText(text().trimmed());
	setEnabled(true);
}

/********** Class C_TableViewDelegate **********/
C_TableViewDelegate::C_TableViewDelegate(QObject *parent) : QItemDelegate(parent)
{
}

QWidget
		*
		C_TableViewDelegate::createEditor(QWidget *parent, const QStyleOptionViewItem & /*option*/, const QModelIndex & /*index*/) const
{
	// These parameters are hard-coded - users will need to adjust them according to their needs
	QDoubleSpinBox *editor = new QDoubleSpinBox(parent);
	editor->setMinimum(0);
	editor->setMaximum(1e10);
	editor->setDecimals(6);
	return editor;
}

/********** Class C_MainWindow **********/

// ******************** //
//	Public functions	//
// ******************** //

MainWindow::MainWindow()
{
	this->boldFont.setBold(true);

	this->threadPreMesh = new C_Thread(this);
	connect(threadPreMesh, SIGNAL(finished()), this, SLOT(threadFinishedPreMesh()));

	this->threadMesh = new C_Thread(this);
	connect(threadMesh, SIGNAL(finished()), this, SLOT(threadFinishedMesh()));

	centralWidget = new QWidget;
	setCentralWidget(centralWidget);

	glWidget = new GLWidget(this);
	connect(glWidget, SIGNAL(SelectionWasMade(unsigned char, unsigned char, unsigned char)), this, SLOT(findSelection(unsigned char, unsigned char, unsigned char)));

	glWidgetArea = new QScrollArea;
	glWidgetArea->setWidget(this->glWidget);
	glWidgetArea->setWidgetResizable(true);
	glWidgetArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
	glWidgetArea->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
	glWidgetArea->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	glWidgetArea->setMinimumSize(100, 100);

	this->createActions();
	this->createMenus();
	this->createDockWindows();

	QGridLayout *centralLayout = new QGridLayout;
	centralLayout->addWidget(glWidgetArea, 0, 0);
	centralWidget->setLayout(centralLayout);

	readSettings();
}

void MainWindow::runThread(QString Attribute)
{
	if (Attribute == "PREMESHJOB")
		this->preMeshJob();
	if (Attribute == "MESHJOB")
		this->MeshJob();
	if (Attribute == "MATERIAL_SELECTION")
		this->materialSelectionJob();
}

void MainWindow::runThreadPool(QString Attribute, int Object1, int Object2, int currentStep, int totalSteps)
{
	//	convex hull calculation
	if (Attribute == "CONVEXHULL")
	{
		emit progress_replace("   > " + QString::number(100 * currentStep / totalSteps) + "% (" + QString::number(currentStep) + "/" + QString::number(totalSteps) + ") " + Model.Surfaces[Object1].Name + " (" + Model.Surfaces[Object1].Type + ")");
		Model.Surfaces[Object1].calculate_normal_vector();
		Model.Surfaces[Object1].rotate(true);
		Model.Surfaces[Object1].calculate_convex_hull();
		Model.Surfaces[Object1].interpolation("ConvexHull", this->interpolationMethod->currentText());
		Model.Surfaces[Object1].rotate(false);
	}
	//	segmentation (coarse and fine)
	if (Attribute == "SEGMENTS")
	{
		emit progress_replace("   > " + QString::number(100 * currentStep / totalSteps) + "% (" + QString::number(currentStep) + "/" + QString::number(totalSteps) + ") " + Model.Polylines[Object1].Name + " (" + Model.Polylines[Object1].Type + ")");
		Model.Polylines[Object1].calculate_segments(false);
		Model.Polylines[Object1].Intersections.clear();
		Model.Polylines[Object1].Path.calculate_min_max();
	}
	if (Attribute == "SEGMENTS_FINE")
	{
		emit progress_replace("   > " + QString::number(100 * currentStep / totalSteps) + "% (" + QString::number(currentStep) + "/" + QString::number(totalSteps) + ") " + Model.Polylines[Object1].Name + " (" + Model.Polylines[Object1].Type + ")");
		Model.Polylines[Object1].calculate_segments(true);
		Model.Polylines[Object1].Path.calculate_min_max();
	}
	// 2D triangulation (coarse and fine)
	if (Attribute == "TRIANGLES")
	{

		emit progress_replace("   > " + QString::number(100 * currentStep / totalSteps) + "% (" + QString::number(currentStep) + "/" + QString::number(totalSteps) + ") " + Model.Surfaces[Object1].Name + " (" + Model.Surfaces[Object1].Type + ")");
		Model.Surfaces[Object1].rotate(true);
		Model.Surfaces[Object1].calculate_triangles(false, Model.preMeshGradient);
		Model.Surfaces[Object1].interpolation("Mesh", this->interpolationMethod->currentText());
		Model.Surfaces[Object1].rotate(false);
		Model.Surfaces[Object1].Intersections.clear();
		Model.Surfaces[Object1].calculate_min_max();
		for (int t = 0; t != Model.Surfaces[Object1].Ts.length(); t++)
		{
			Model.Surfaces[Object1].Ts[t].calculate_min_max();
			Model.Surfaces[Object1].Ts[t].setNormalVector();
		}
	}
	if (Attribute == "TRIANGLES_FINE")
	{
		emit progress_replace("   > " + QString::number(100 * currentStep / totalSteps) + "% (" + QString::number(currentStep) + "/" + QString::number(totalSteps) + ") " + Model.Surfaces[Object1].Name + " (" + Model.Surfaces[Object1].Type + ")");
		Model.Surfaces[Object1].calculate_normal_vector();
		Model.Surfaces[Object1].rotate(true);
		Model.Surfaces[Object1].separate_Constraints();
		Model.Surfaces[Object1].calculate_triangles(true, Model.meshGradient);
		Model.Surfaces[Object1].interpolation("Mesh", this->interpolationMethod->currentText());
		Model.Surfaces[Object1].rotate(false);
		for (int t = 0; t != Model.Surfaces[Object1].Ts.length(); t++)
		{
			Model.Surfaces[Object1].Ts[t].calculate_min_max();
			Model.Surfaces[Object1].Ts[t].setNormalVector();
		}
	}
	//	intersection - surfaces-surfaces
	if (Attribute == "INTERSECTION_MESH_MESH")
	{
		emit progress_replace("   > " + QString::number(100 * currentStep / totalSteps) + "% (" + QString::number(currentStep) + "/" + QString::number(totalSteps) + ") " + Model.Surfaces[Object1].Name + " with " + Model.Surfaces[Object2].Name);
		Model.calculate_int_polyline(Object1, Object2);
	}
	//	intersection - polylines-surfaces
	if (Attribute == "INTERSECTION_POLYLINE_MESH")
	{
		emit progress_replace("   > " + QString::number(100 * currentStep / totalSteps) + "% (" + QString::number(currentStep) + "/" + QString::number(totalSteps) + ") " + Model.Polylines[Object1].Name + " with " + Model.Surfaces[Object2].Name);
		Model.calculate_int_point(Object1, Object2);
	}
	//	intersections - triple points
	if (Attribute == "INTERSECTION_TRIPLEPOINTS")
	{
		emit progress_replace("   > " + QString::number(100 * currentStep / totalSteps) + "% (" + QString::number(currentStep) + "/" + QString::number(totalSteps) + ") ");
		Model.calculate_int_triplepoints(Object1, Object2);
	}
}

// ******************** //
//	Protected functions	//
// ******************** //

void MainWindow::closeEvent(QCloseEvent *)
{
	this->writeSettings();

	/* Stop worker threads. */
	this->threadPreMesh->quit();
	this->threadPreMesh->wait();
	this->threadMesh->quit();
	this->threadMesh->wait();
}

void MainWindow::keyPressEvent(QKeyEvent *event)
{
	if (event->key() == Qt::Key_P)
		if (glWidget->isPerspec())
			glWidget->setPerspec(false);
		else
			glWidget->setPerspec(true);

	if (event->key() == Qt::Key_0)
		glWidget->resetView();

	if (event->key() == Qt::Key_W)
		glWidget->changeViewport(GLWidget::X_ROT, 5);

	if (event->key() == Qt::Key_S)
		glWidget->changeViewport(GLWidget::X_ROT, -5);

	if (event->key() == Qt::Key_A)
		glWidget->changeViewport(GLWidget::Z_ROT, 5);

	if (event->key() == Qt::Key_D)
		glWidget->changeViewport(GLWidget::Z_ROT, -5);

	if (event->key() == Qt::Key_M)
	{
		int index = rotationCombo->findData(RotationMode::MODEL);
		if (index != -1)
			rotationCombo->setCurrentIndex(index);
	}

	if (event->key() == Qt::Key_O)
	{
		int index = rotationCombo->findData(RotationMode::OBJECT);
		if (index != -1)
			rotationCombo->setCurrentIndex(index);
	}

	if (event->key() == Qt::Key_V)
	{
		int index = rotationCombo->findData(RotationMode::VIEWPORT);
		if (index != -1)
			rotationCombo->setCurrentIndex(index);
	}
}

// ******************** //
//	Private functions	//
// ******************** //

void MainWindow::createActions()
{
	//	action - open
	openAct = new QAction(QIcon(":/images/open.png"), tr("&Open..."), this);
	openAct->setShortcuts(QKeySequence::Open);
	openAct->setStatusTip(tr("Open an existing Model"));
	connect(openAct, SIGNAL(triggered()), this, SLOT(open()));
	//	action - save
	saveAct = new QAction(QIcon(":/images/save.png"), tr("&Save"), this);
	saveAct->setShortcuts(QKeySequence::Save);
	saveAct->setStatusTip(tr("Save the existing Model"));
	connect(saveAct, SIGNAL(triggered()), this, SLOT(save()));
	//	action - save as
	saveAsAct = new QAction(QIcon(":/images/saveAs.png"), tr("Save&As..."), this);
	saveAsAct->setShortcut(tr("F11"));
	saveAsAct->setStatusTip(tr("Save the existing Model"));
	connect(saveAsAct, SIGNAL(triggered()), this, SLOT(saveAs()));
	//	action - save screen
	saveScreenAct = new QAction(QIcon(":/images/camera.png"), tr("Save Screenshot..."), this);
	saveScreenAct->setShortcut(tr("F12"));
	saveScreenAct->setStatusTip(tr("Save Screenshot"));
	connect(saveScreenAct, SIGNAL(triggered()), this, SLOT(saveScreen()));
	//	action - importGoCad
	importGoCadAct = new QAction(tr("GoCad..."), this);
	importGoCadAct->setStatusTip(tr("Import a GoCad Model"));
	connect(importGoCadAct, SIGNAL(triggered()), this, SLOT(importGoCad()));
	//	action - export 2 Feflow
	exportFeFlowAct = new QAction(QIcon(":/images/feflow.png"), tr("FeFlow..."), this);
	exportFeFlowAct->setStatusTip(tr("Export the final surface mesh to FeFlow"));
	connect(exportFeFlowAct, SIGNAL(triggered()), this, SLOT(exportFeFlow()));
	//	action - export 2 OGS + TIN
	exportOGSAct = new QAction(QIcon(":/images/ogs.png"), tr("OGS..."), this);
	exportOGSAct->setStatusTip(tr("Export the final mesh to OpenGeosys"));
	connect(exportOGSAct, SIGNAL(triggered()), this, SLOT(exportOGS()));
	exportTINAct = new QAction(QIcon(":/images/ogs.png"), tr("TIN..."), this);
	exportTINAct->setStatusTip(tr("Export TIN of selected surfaces"));
	connect(exportTINAct, SIGNAL(triggered()), this, SLOT(exportTIN()));
	//	action - export 2 Abaqus ANSYS
	exportABAQUSAct = new QAction(QIcon(":/images/abaqus.png"), tr("ABAQUS..."), this);
	exportABAQUSAct->setStatusTip(tr("Export the final mesh to ABAQUS"));
	connect(exportABAQUSAct, SIGNAL(triggered()), this, SLOT(exportABAQUS()));
	//	action - export 2 Comsol
	exportCOMSOLAct = new QAction(QIcon(":/images/comsol.png"), tr("COMSOL..."), this);
	exportCOMSOLAct->setStatusTip(tr("Export the final mesh to COMSOL"));
	connect(exportCOMSOLAct, SIGNAL(triggered()), this, SLOT(exportCOMSOL()));
	//	action - export 2 Paraview 3D and 2D
	exportVTU3DAct = new QAction(QIcon(":/images/paraview.png"), tr("PARAVIEW..."), this);
	exportVTU3DAct->setStatusTip(tr("Export the final mesh to Paraview"));
	connect(exportVTU3DAct, SIGNAL(triggered()), this, SLOT(exportVTU3D()));
	exportVTU2DAct = new QAction(QIcon(":/images/paraview.png"), tr("PARAVIEW..."), this);
	exportVTU2DAct->setStatusTip(tr("Export VTU of selected surfaces"));
	connect(exportVTU2DAct, SIGNAL(triggered()), this, SLOT(exportVTU2D()));
	//	action - export 2 Exodus
	exportEXODUSAct = new QAction(QIcon(":/images/moose.png"), tr("MOOSE..."), this);
	exportEXODUSAct->setStatusTip(tr("Export the final mesh to MOOSE"));
#ifdef NOEXODUS
	exportEXODUSAct->setDisabled(true);
#endif
	connect(exportEXODUSAct, SIGNAL(triggered()), this, SLOT(exportEXODUS()));
	//	action - add unit
	addUnitAct = new QAction(QIcon(":/images/units.png"), tr("Unit..."), this);
	addUnitAct->setStatusTip(tr("Add a Unit to the current project"));
	connect(addUnitAct, SIGNAL(triggered()), this, SLOT(addUnit()));
	//	action - add fault
	addFaultAct = new QAction(QIcon(":/images/faults.png"), tr("Fault..."), this);
	addFaultAct->setStatusTip(tr("Add a Fault to the current project"));
	connect(addFaultAct, SIGNAL(triggered()), this, SLOT(addFault()));
	//	action - add border
	addBorderAct = new QAction(QIcon(":/images/borders.png"), tr("Border..."), this);
	addBorderAct->setStatusTip(tr("Add a Border to the current project"));
	connect(addBorderAct, SIGNAL(triggered()), this, SLOT(addBorder()));
	//	action - add well
	addWellAct = new QAction(QIcon(":/images/wells.png"), tr("Well..."), this);
	addWellAct->setStatusTip(tr("Add a Well to the current project"));
	connect(addWellAct, SIGNAL(triggered()), this, SLOT(addWell()));
	//	action - delete
	deleteSurfaceAct = new QAction(QIcon(":/images/delete_icon.png"), tr("Surfaces..."), this);
	deleteSurfaceAct->setStatusTip(tr("Delete selected surfaces"));
	connect(deleteSurfaceAct, SIGNAL(triggered()), this, SLOT(deleteSurface()));
	//	action - import PLC object
	importPLCAct = new QAction(tr("PLC..."), this);
	importPLCAct->setStatusTip(tr("Import a Piecewise Linear Complex"));
	connect(importPLCAct, SIGNAL(triggered()), this, SLOT(settingCallByMenu()));
	// action - mesh PLC object
	meshPLCAct = new QAction(tr("mesh PLC"), this);
	meshPLCAct->setStatusTip(tr("Mesh the Piecewise Linear Complex"));
	connect(meshPLCAct, SIGNAL(triggered()), this, SLOT(settingCallByMenu()));
	//	action - select refinement
	editRefinement = new QAction(tr("Refinement"), this);
	editRefinement->setStatusTip(tr("Defines object specific refinement"));
	connect(editRefinement, SIGNAL(triggered()), this, SLOT(settingCallByMenu()));
	//	action - select interpolation
	editInterpolation = new QAction(tr("Interpolation"), this);
	editInterpolation->setStatusTip(tr("Defines Surface Interpolation Method"));
	connect(editInterpolation, SIGNAL(triggered()), this, SLOT(settingCallByMenu()));
	// action - change gradient
	preMeshEditGradient = new QAction(tr("Gradient"), this);
	preMeshEditGradient->setStatusTip(tr("Change gradient"));
	preMeshEditGradient->setProperty("phase", "PreMesh");
	connect(preMeshEditGradient, SIGNAL(triggered()), this, SLOT(settingCallByMenu()));
	//	action - preMesh
	editPreMesh = new QAction(tr("> Execute PreMesh <"), this);
	editPreMesh->setFont(boldFont);
	editPreMesh->setStatusTip(tr("Generates the PLC"));
	connect(editPreMesh, SIGNAL(triggered()), this, SLOT(settingCallByMenu()));
	//	action - model selection
	editSelection = new QAction(tr("Selection"), this);
	editSelection->setStatusTip(tr("Selection of constaints"));
	connect(editSelection, SIGNAL(triggered()), this, SLOT(settingCallByMenu()));
	// action - change gradient
	meshEditGradient = new QAction(tr("Gradient"), this);
	meshEditGradient->setStatusTip(tr("Change gradient"));
	meshEditGradient->setProperty("phase", "Mesh");
	connect(meshEditGradient, SIGNAL(triggered()), this, SLOT(settingCallByMenu()));
	//	action - set Material 1D-2D-3D
	editMaterial1d = new QAction(QIcon(":/images/1D.png"), tr("1D"), this);
	editMaterial1d->setStatusTip(tr("Changes the 1D materials"));
	connect(editMaterial1d, SIGNAL(triggered()), this, SLOT(settingCallByMenu()));
	editMaterial2d = new QAction(QIcon(":/images/2D.png"), tr("2D"), this);
	editMaterial2d->setStatusTip(tr("Changes the 2D materials"));
	connect(editMaterial2d, SIGNAL(triggered()), this, SLOT(settingCallByMenu()));
	editMaterial3d = new QAction(QIcon(":/images/3D.png"), tr("3D"), this);
	editMaterial3d->setStatusTip(tr("Changes the 3D materials"));
	connect(editMaterial3d, SIGNAL(triggered()), this, SLOT(settingCallByMenu()));
	//	action - tetgen options
	editTetgen = new QAction(tr("Tetgen"), this);
	editTetgen->setStatusTip(tr("Tetgen command line option"));
	connect(editTetgen, SIGNAL(triggered()), this, SLOT(settingCallByMenu()));
	//	action - Mesh
	editMesh = new QAction(tr("> Execute Mesh <"), this);
	editMesh->setFont(boldFont);
	editMesh->setStatusTip(tr("Generates the tetreahedral mesh"));
	connect(editMesh, SIGNAL(triggered()), this, SLOT(settingCallByMenu()));
	//	action - exit
	exitAct = new QAction(tr("E&xit"), this);
	exitAct->setShortcut(tr("Ctrl+Q"));
	exitAct->setStatusTip(tr("Exit the application"));
	connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));
	//	action - visualize orientation axis
	viewAxisAct = new QAction(tr("Orientation Axis"), this);
	viewAxisAct->setStatusTip(tr("Turns orientation axis on/off"));
	viewAxisAct->setCheckable(true);
	viewAxisAct->setChecked(true);
	connect(viewAxisAct, SIGNAL(triggered()), this, SLOT(viewAxis()));
	//	action - print about message*
	aboutAct = new QAction(tr("&About"), this);
	aboutAct->setStatusTip(tr("Show the application's About box"));
	connect(aboutAct, SIGNAL(triggered()), this, SLOT(about()));
}

void MainWindow::createMenus()
{
	//	First menu bar that recalls actions to open, save(As), import, export, add, delete and exit
	fileMenu = menuBar()->addMenu(tr("&File"));
	fileMenu->addAction(openAct);
	fileMenu->addAction(saveAct);
	fileMenu->addAction(saveAsAct);
	fileMenu->addAction(saveScreenAct);
	fileMenu->addSeparator();
	fileMenuImport = fileMenu->addMenu(tr("Import"));
	fileMenuImport->addAction(this->importGoCadAct);
	fileMenuExport = fileMenu->addMenu(tr("Export as ... "));
	fileMenuExport3D = fileMenuExport->addMenu(tr("3D mesh"));
	fileMenuExport3D->addAction(this->exportFeFlowAct);
	fileMenuExport3D->addAction(this->exportOGSAct);
	fileMenuExport3D->addAction(this->exportCOMSOLAct);
	fileMenuExport3D->addAction(this->exportABAQUSAct);
	fileMenuExport3D->addAction(this->exportEXODUSAct);
	fileMenuExport3D->addAction(this->exportVTU3DAct);
	fileMenuExport2D = fileMenuExport->addMenu(tr("2D mesh"));
	fileMenuExport2D->addAction(this->exportTINAct);
	fileMenuExport2D->addAction(this->exportVTU2DAct);
	fileMenu->addSeparator();
	fileMenuAdd = fileMenu->addMenu(tr("Add"));
	fileMenuAdd->addAction(this->addUnitAct);
	fileMenuAdd->addAction(this->addFaultAct);
	fileMenuAdd->addAction(this->addBorderAct);
	fileMenuAdd->addAction(this->addWellAct);
	fileMenu->addSeparator();
	fileMenuDelete = fileMenu->addMenu(tr("Delete"));
	fileMenuDelete->addAction(this->deleteSurfaceAct);
	fileMenu->addSeparator();
	fileMenu->addAction(exitAct);
	//  Third menu bar that recalls actions to import a PLC object and mesh it via direct tetgen call
	meshPLCMenu = menuBar()->addMenu(tr("&mesh PLC"));
	meshPLCMenu->addAction(this->importPLCAct);
	meshPLCMenu->addAction(this->editTetgen);
	meshPLCMenu->addAction(this->meshPLCAct);
	//	Second menu bar that recalls actions to edit material for 1D, 2D, and 3D cases
	editMenu = menuBar()->addMenu(tr("&Edit"));
	this->editMenu->addAction(this->editRefinement);
	this->editMenu->addAction(this->editInterpolation);
	editMenu->addAction(this->preMeshEditGradient);
	editMenu->addAction(this->editPreMesh);
	this->editMenuMaterial = this->editMenu->addMenu(tr("Material"));
	this->editMenuMaterial->addAction(this->editMaterial1d);
	this->editMenuMaterial->addAction(this->editMaterial2d);
	this->editMenuMaterial->addAction(this->editMaterial3d);
	editMenu->addAction(this->editSelection);
	editMenu->addAction(this->meshEditGradient);
	this->editMenu->addAction(this->editTetgen);
	this->editMenu->addAction(this->editMesh);
	//	Third menu bar that recalls actions to view and print help message
	viewMenu = menuBar()->addMenu(tr("&View"));
	viewMenu->addAction(this->viewAxisAct);
	menuBar()->addSeparator();
	helpMenu = menuBar()->addMenu(tr("&Help"));
	helpMenu->addAction(aboutAct);
}

void MainWindow::createDockWindows()
{
	this->createStatusDock();
	this->createViewDock();
	this->createMeshQualityDock();
}

void MainWindow::createStatusDock()
{
	statusDock = new QDockWidget(tr("Status"), this);
	statusDock->setObjectName(tr("StatusDock"));
	statusDock->setAllowedAreas(Qt::BottomDockWidgetArea);
	statusTE = new QTextEdit(statusDock);
	statusTE->setText(">Welcome to MeshIT");
	statusTE->setReadOnly(true);
	connect(this, SIGNAL(progress_append(QString)), statusTE, SLOT(append(QString)));
	connect(this, SIGNAL(progress_replace(QString)), this, SLOT(replace(QString)));
	connect(&Model, SIGNAL(PrintError(QString)),
					this, SLOT(printErrorMessage(QString)));
	statusDock->setWidget(this->statusTE);
	addDockWidget(Qt::BottomDockWidgetArea, statusDock);
	viewMenu->addAction(statusDock->toggleViewAction());
}

void MainWindow::createViewDock()
{
	viewDock = new QDockWidget(tr("View"), this);
	viewDock->setObjectName(tr("ViewDock"));
	viewDock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

	viewWidget = new QWidget(viewDock);
	viewVBox = new QVBoxLayout(viewWidget);
	viewTree = new QTreeWidget(viewWidget);
	viewTree->setHeaderHidden(true);

	//	add unit widget tree to the View dock
	QTreeWidgetItem *unitsTitle = new QTreeWidgetItem(viewTree);
	unitsTitle->setText(0, "Units");
	unitsTitle->setIcon(0, QIcon(":/images/units.png"));
	unitsNames = new QTreeWidgetItem(unitsTitle);
	unitsNamesCB = new QComboBox(viewTree);
	viewTree->setItemWidget(unitsNames, 0, unitsNamesCB);
	//	scattered data (input) points
	unitsScatteredDataPoints = new QTreeWidgetItem(unitsTitle);
	unitsScatteredDataPoints->setText(0, "Scattered Data Points");
	unitsScatteredDataPoints->setCheckState(0, Qt::Unchecked);
	//	convex hulls
	unitsConvexHull = new QTreeWidgetItem(unitsTitle);
	unitsConvexHull->setText(0, "Convex Hull");
	unitsConvexHull->setCheckState(0, Qt::Unchecked);
	//	surfaces - faces and edges
	unitsSurface = new QTreeWidgetItem(unitsTitle);
	unitsSurface->setText(0, "Surface");
	unitsFaces = new QTreeWidgetItem(unitsSurface);
	unitsFaces->setText(0, "Faces");
	unitsFaces->setCheckState(0, Qt::Unchecked);
	unitsEdges = new QTreeWidgetItem(unitsSurface);
	unitsEdges->setText(0, "Edges");
	unitsEdges->setCheckState(0, Qt::Unchecked);
	//	intersections - edges and vertices
	unitsIntersection = new QTreeWidgetItem(unitsTitle);
	unitsIntersection->setText(0, "Intersections");
	unitsIntersectionEdges = new QTreeWidgetItem(unitsIntersection);
	unitsIntersectionEdges->setText(0, "Edges");
	unitsIntersectionEdges->setCheckState(0, Qt::Unchecked);
	unitsIntersectionVertices = new QTreeWidgetItem(unitsIntersection);
	unitsIntersectionVertices->setText(0, "Vertices");
	unitsIntersectionVertices->setCheckState(0, Qt::Unchecked);
	connect(unitsNamesCB, SIGNAL(currentIndexChanged(QString)), this, SLOT(setUShowGBox(QString)));
	connect(viewTree, SIGNAL(itemClicked(QTreeWidgetItem *, int)), this, SLOT(setUMeshes()));
	//	add fault tree widget to the View dock
	QTreeWidgetItem *faultsTitle = new QTreeWidgetItem(viewTree);
	faultsTitle->setText(0, "Faults");
	faultsTitle->setIcon(0, QIcon(":/images/faults.png"));
	faultsNames = new QTreeWidgetItem(faultsTitle);
	faultsNamesCB = new QComboBox(viewTree);
	viewTree->setItemWidget(faultsNames, 0, faultsNamesCB);
	//	scattered data (input) points
	faultsScatteredDataPoints = new QTreeWidgetItem(faultsTitle);
	faultsScatteredDataPoints->setText(0, "Scattered Data Points");
	faultsScatteredDataPoints->setCheckState(0, Qt::Unchecked);
	//	convex hulls
	faultsConvexHull = new QTreeWidgetItem(faultsTitle);
	faultsConvexHull->setText(0, "Convex Hull");
	faultsConvexHull->setCheckState(0, Qt::Unchecked);
	//	surfaces - faces and edges
	faultsSurface = new QTreeWidgetItem(faultsTitle);
	faultsSurface->setText(0, "Surface");
	faultsFaces = new QTreeWidgetItem(faultsSurface);
	faultsFaces->setText(0, "Faces");
	faultsFaces->setCheckState(0, Qt::Unchecked);
	faultsEdges = new QTreeWidgetItem(faultsSurface);
	faultsEdges->setText(0, "Edges");
	faultsEdges->setCheckState(0, Qt::Unchecked);
	//	intersections - edges and vertices
	faultsIntersection = new QTreeWidgetItem(faultsTitle);
	faultsIntersection->setText(0, "Intersections");
	faultsIntersectionEdges = new QTreeWidgetItem(faultsIntersection);
	faultsIntersectionEdges->setText(0, "Edges");
	faultsIntersectionEdges->setCheckState(0, Qt::Unchecked);
	faultsIntersectionVertices = new QTreeWidgetItem(faultsIntersection);
	faultsIntersectionVertices->setText(0, "Vertices");
	faultsIntersectionVertices->setCheckState(0, Qt::Unchecked);
	connect(faultsNamesCB, SIGNAL(currentIndexChanged(QString)), this, SLOT(setFShowGBox(QString)));
	connect(viewTree, SIGNAL(itemClicked(QTreeWidgetItem *, int)), this, SLOT(setFMeshes()));
	//	add border tree widget to the View dock
	QTreeWidgetItem *bordersTitle = new QTreeWidgetItem(viewTree);
	bordersTitle->setText(0, "Borders");
	bordersTitle->setIcon(0, QIcon(":/images/borders.png"));
	bordersNames = new QTreeWidgetItem(bordersTitle);
	bordersNamesCB = new QComboBox(viewTree);
	viewTree->setItemWidget(bordersNames, 0, bordersNamesCB);
	//	scattered data (input) points
	bordersScatteredDataPoints = new QTreeWidgetItem(bordersTitle);
	bordersScatteredDataPoints->setText(0, "Scattered Data Points");
	bordersScatteredDataPoints->setCheckState(0, Qt::Unchecked);
	//	convex hulls
	bordersConvexHull = new QTreeWidgetItem(bordersTitle);
	bordersConvexHull->setText(0, "Convex Hull");
	bordersConvexHull->setCheckState(0, Qt::Unchecked);
	//	surfaces - faces and edges
	bordersSurface = new QTreeWidgetItem(bordersTitle);
	bordersSurface->setText(0, "Surface");
	bordersFaces = new QTreeWidgetItem(bordersSurface);
	bordersFaces->setText(0, "Faces");
	bordersFaces->setCheckState(0, Qt::Unchecked);
	bordersEdges = new QTreeWidgetItem(bordersSurface);
	bordersEdges->setText(0, "Edges");
	bordersEdges->setCheckState(0, Qt::Unchecked);
	//	intersections - edges and vertices
	bordersIntersection = new QTreeWidgetItem(bordersTitle);
	bordersIntersection->setText(0, "Intersections");
	bordersIntersectionEdges = new QTreeWidgetItem(bordersIntersection);
	bordersIntersectionEdges->setText(0, "Edges");
	bordersIntersectionEdges->setCheckState(0, Qt::Unchecked);
	bordersIntersectionVertices = new QTreeWidgetItem(bordersIntersection);
	bordersIntersectionVertices->setText(0, "Vertices");
	bordersIntersectionVertices->setCheckState(0, Qt::Unchecked);
	connect(bordersNamesCB, SIGNAL(currentIndexChanged(QString)), this, SLOT(setBShowGBox(QString)));
	connect(viewTree, SIGNAL(itemClicked(QTreeWidgetItem *, int)), this, SLOT(setBMeshes()));
	//	add well tree widget to the View dock
	QTreeWidgetItem *wellsTitle = new QTreeWidgetItem(viewTree);
	wellsTitle->setText(0, "Wells");
	wellsTitle->setIcon(0, QIcon(":/images/wells.png"));
	wellsNames = new QTreeWidgetItem(wellsTitle);
	wellsNamesCB = new QComboBox(viewTree);
	viewTree->setItemWidget(wellsNames, 0, wellsNamesCB);
	//	scattered data (input) points
	wellsScatteredDataPoints = new QTreeWidgetItem(wellsTitle);
	wellsScatteredDataPoints->setText(0, "Scattered Data Points");
	wellsScatteredDataPoints->setCheckState(0, Qt::Unchecked);
	//	paths - edges and vertices
	wellsPath = new QTreeWidgetItem(wellsTitle);
	wellsPath->setText(0, "Path");
	wellsEdges = new QTreeWidgetItem(wellsPath);
	wellsEdges->setText(0, "Edges");
	wellsEdges->setCheckState(0, Qt::Unchecked);
	wellsVertices = new QTreeWidgetItem(wellsPath);
	wellsVertices->setText(0, "Vertices");
	wellsVertices->setCheckState(0, Qt::Unchecked);
	//	intersections - vertices
	wellsIntersectionVertices = new QTreeWidgetItem(wellsTitle);
	wellsIntersectionVertices->setText(0, "Intersections");
	wellsIntersectionVertices->setCheckState(0, Qt::Unchecked);
	connect(wellsNamesCB, SIGNAL(currentIndexChanged(QString)), this, SLOT(setWShowGBox(QString)));
	connect(viewTree, SIGNAL(itemClicked(QTreeWidgetItem *, int)), this, SLOT(setWMeshes()));
	//	add 3D material tree widget to the View dock
	matsTitle = new QTreeWidgetItem(viewTree);
	matsTitle->setText(0, "Materials");
	matsTitle->setHidden(true);
	matsTitle->setIcon(0, QIcon(":/images/tetrahedron.png"));
	Model.drawTets = false;
	matsNames = new QTreeWidgetItem(matsTitle);
	matsNamesCB = new QComboBox(viewTree);
	viewTree->setItemWidget(matsNames, 0, matsNamesCB);
	//	faces
	matsFaces = new QTreeWidgetItem(matsTitle);
	matsFaces->setText(0, "Faces");
	matsFaces->setCheckState(0, Qt::Unchecked);
	//	edges
	matsEdges = new QTreeWidgetItem(matsTitle);
	matsEdges->setText(0, "Edges");
	matsEdges->setCheckState(0, Qt::Unchecked);
	connect(matsNamesCB, SIGNAL(currentIndexChanged(QString)), this, SLOT(setMShowGBox(QString)));
	connect(viewTree, SIGNAL(itemClicked(QTreeWidgetItem *, int)), this, SLOT(setMMeshes()));
	//	add meshView box to the View dock
	viewVBox->addWidget(viewTree);
	tetViewGBox = new QGroupBox(tr("MeshView"), viewWidget);
	tetViewGrid = new QGridLayout(tetViewGBox);
	//	x-direction - cut
	xCutEnable = new QCheckBox(tr("x CutPlane"), tetViewGBox);
	xCutEnable->setChecked(false);
	connect(xCutEnable, SIGNAL(clicked(bool)), this, SLOT(callMakeTets()));
	tetViewGrid->addWidget(xCutEnable, 1, 0, 1, 1);
	//	x-direction - flip
	xDirChange = new QCheckBox(tr("Flip"), tetViewGBox);
	xDirChange->setChecked(false);
	connect(xDirChange, SIGNAL(clicked(bool)), this, SLOT(callMakeTets()));
	tetViewGrid->addWidget(xDirChange, 1, 1, 1, 1);
	//	x-direction - slider
	xCutSlider = new QSlider(Qt::Horizontal);
	xCutSlider->setRange(-256, 256);
	xCutSlider->setValue(0);
	connect(xCutSlider, SIGNAL(valueChanged(int)), this, SLOT(callMakeTets()));
	tetViewGrid->addWidget(xCutSlider, 2, 0, 1, 2);
	//	y-direction - cut
	yCutEnable = new QCheckBox(tr("y CutPlane"), tetViewGBox);
	yCutEnable->setChecked(false);
	connect(yCutEnable, SIGNAL(clicked(bool)), this, SLOT(callMakeTets()));
	tetViewGrid->addWidget(yCutEnable, 3, 0, 1, 1);
	//	y-direction - flip
	yDirChange = new QCheckBox(tr("Flip"), tetViewGBox);
	yDirChange->setChecked(false);
	connect(yDirChange, SIGNAL(clicked(bool)), this, SLOT(callMakeTets()));
	tetViewGrid->addWidget(yDirChange, 3, 1, 1, 1);
	//	y-direction - slider
	yCutSlider = new QSlider(Qt::Horizontal);
	yCutSlider->setRange(-256, 256);
	yCutSlider->setValue(0);
	connect(yCutSlider, SIGNAL(valueChanged(int)), this, SLOT(callMakeTets()));
	tetViewGrid->addWidget(yCutSlider, 4, 0, 1, 2);
	//	z-direction - cut
	zCutEnable = new QCheckBox(tr("z CutPlane"), tetViewGBox);
	zCutEnable->setChecked(false);
	connect(zCutEnable, SIGNAL(clicked(bool)), this, SLOT(callMakeTets()));
	tetViewGrid->addWidget(zCutEnable, 5, 0, 1, 1);
	//	z-direction - flip
	zDirChange = new QCheckBox(tr("Flip"), tetViewGBox);
	zDirChange->setChecked(false);
	connect(zDirChange, SIGNAL(clicked(bool)), this, SLOT(callMakeTets()));
	tetViewGrid->addWidget(zDirChange, 5, 1, 1, 1);
	//	z-direction - slider
	zCutSlider = new QSlider(Qt::Horizontal);
	zCutSlider->setRange(-256, 256);
	zCutSlider->setValue(0);
	connect(zCutSlider, SIGNAL(valueChanged(int)), this, SLOT(callMakeTets()));
	tetViewGrid->addWidget(zCutSlider, 6, 0, 1, 2);
	tetViewGBox->setLayout(tetViewGrid);
	viewVBox->addWidget(tetViewGBox);

	/* Add widget used to select the rotation center. */
	QGroupBox *rotationBox = new QGroupBox(tr("Rotation center"), viewWidget);
	QGridLayout *rotationLayout = new QGridLayout(rotationBox);
	rotationCombo = new QComboBox();
	rotationCombo->insertItem(0, "Model", RotationMode::MODEL);
	rotationCombo->setItemData(0, "Sets the rotation center to the center of "
																"the entire model.\nKeyboard shortcut: M",
														 Qt::ToolTipRole);
	rotationCombo->insertItem(1, "Object", RotationMode::OBJECT);
	rotationCombo->setItemData(1, "Sets the rotation center based on the "
																"visible objects only.\nKeyboard shortcut: O",
														 Qt::ToolTipRole);
	rotationCombo->insertItem(2, "Viewport", RotationMode::VIEWPORT);
	rotationCombo->setItemData(2, "Sets the rotation center to the center of "
																"the viewport.\nKeyboard shortcut: V",
														 Qt::ToolTipRole);

	connect(rotationCombo, SIGNAL(currentIndexChanged(int)),
					this, SLOT(onRotationCenterChange(int)));

	rotationLayout->addWidget(rotationCombo);
	viewVBox->addWidget(rotationBox);

	/* Navigation sensitivity. */
	QLabel *label = new QLabel("Navigation sensitivity");
	viewVBox->addWidget(label);

	QSlider *sensitivitySlider = new QSlider(Qt::Horizontal);
	sensitivitySlider->setRange(1, 100);
	sensitivitySlider->setValue(glWidget->zoomSensitivity * 100);
	connect(sensitivitySlider, SIGNAL(valueChanged(int)),
					this, SLOT(onSensitivityChange()));
	viewVBox->addWidget(sensitivitySlider);

	//	add info box to the View Dock
	infoGBox = new QGroupBox(tr("Info"), viewWidget);
	connect(&Model, SIGNAL(ModelInfoChanged()), this, SLOT(updateInfo()));
	infoGBox->setDisabled(true);
	//	layout
	infoGrid = new QGridLayout(infoGBox);
	infoX = new QLabel("X", infoGBox);
	infoGrid->addWidget(infoX, 1, 0, 1, 1);
	infoY = new QLabel("Y", infoGBox);
	infoGrid->addWidget(infoY, 2, 0, 1, 1);
	infoZ = new QLabel("Z", infoGBox);
	infoGrid->addWidget(infoZ, 3, 0, 1, 1);
	infoMin = new QLabel("min", infoGBox);
	infoGrid->addWidget(infoMin, 0, 1, 1, 1);
	infoMax = new QLabel("max", infoGBox);
	infoGrid->addWidget(infoMax, 0, 2, 1, 1);
	infoRan = new QLabel("Range", infoGBox);
	infoGrid->addWidget(infoRan, 0, 3, 1, 1);
	//	x-dimension information
	infoMinX = new QLabel("0", infoGBox);
	infoGrid->addWidget(infoMinX, 1, 1, 1, 1);
	infoMaxX = new QLabel("0", infoGBox);
	infoGrid->addWidget(infoMaxX, 1, 2, 1, 1);
	infoRanX = new QLabel("0", infoGBox);
	infoGrid->addWidget(infoRanX, 1, 3, 1, 1);
	//	y-dimension information
	infoMinY = new QLabel("0", infoGBox);
	infoGrid->addWidget(infoMinY, 2, 1, 1, 1);
	infoMaxY = new QLabel("0", infoGBox);
	infoGrid->addWidget(infoMaxY, 2, 2, 1, 1);
	infoRanY = new QLabel("0", infoGBox);
	infoGrid->addWidget(infoRanY, 2, 3, 1, 1);
	//	z-dimension information
	infoMinZ = new QLabel("0", infoGBox);
	infoGrid->addWidget(infoMinZ, 3, 1, 1, 1);
	infoMaxZ = new QLabel("0", infoGBox);
	infoGrid->addWidget(infoMaxZ, 3, 2, 1, 1);
	infoRanZ = new QLabel("0", infoGBox);
	infoGrid->addWidget(infoRanZ, 3, 3, 1, 1);
	//	tetrahedrons information
	infoTetNumber = new QLabel("no. Tetrahedrons: " + QString::number(Model.Ts.length()), infoGBox);
	infoGrid->addWidget(infoTetNumber, 4, 0, 1, 4);
	infoGBox->setLayout(infoGrid);
	viewVBox->addWidget(infoGBox);
	viewSpacer = new QSpacerItem(0, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
	viewVBox->addItem(viewSpacer);
	viewWidget->setLayout(viewVBox);

	viewDock->setWidget(this->viewWidget);
	addDockWidget(Qt::RightDockWidgetArea, viewDock);
	viewMenu->addAction(viewDock->toggleViewAction());
}

void MainWindow::createMeshQualityDock()
{
	meshDock = new QDockWidget(tr("Mesh Quality"), this);
	meshDock->setObjectName(tr("MeshQualityDock"));
	meshDock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

	meshWidget = new QWidget(meshDock);
	meshVBox = new QVBoxLayout(meshWidget);
	meshTree = new QTreeWidget(meshWidget);
	this->meshTree->setHeaderHidden(true);

	//	add preMesh tree widget to the mesh quality dock
	this->meshPreMesh = new QTreeWidgetItem(meshTree);
	this->meshPreMesh->setText(0, "PreMesh");
	//	preMesh - refinement
	this->meshPreMeshRefinement = new QTreeWidgetItem(meshPreMesh);
	this->meshPreMeshRefinement->setText(0, "Refinement");
	//	preMesh - interpolation
	this->meshPreMeshInterpolation = new QTreeWidgetItem(meshPreMesh);
	this->meshPreMeshInterpolation->setText(0, "Interpolation");
	// preMesh - gradient
	this->meshPreMeshGradient = new QTreeWidgetItem(meshPreMesh);
	this->meshPreMeshGradient->setText(0, "Gradient");
	//	preMesh - execute
	this->meshPreMeshExecute = new QTreeWidgetItem(meshPreMesh);
	this->meshPreMeshExecute->setFont(0, boldFont);
	this->meshPreMeshExecute->setText(0, "> Execute PreMesh <");
	//	add Mesh tree widget to the mesh quality dock
	this->meshMesh = new QTreeWidgetItem(meshTree);
	this->meshMesh->setText(0, "Mesh");
	// mesh - material 1D, 2D, 3D
	this->meshMeshMaterial = new QTreeWidgetItem(meshMesh);
	this->meshMeshMaterial->setText(0, "Material");
	this->meshMeshMaterial1d = new QTreeWidgetItem(this->meshMeshMaterial);
	this->meshMeshMaterial1d->setText(0, "1D");
	this->meshMeshMaterial1d->setIcon(0, QIcon(":/images/1D.png"));
	this->meshMeshMaterial2d = new QTreeWidgetItem(this->meshMeshMaterial);
	this->meshMeshMaterial2d->setText(0, "2D");
	this->meshMeshMaterial2d->setIcon(0, QIcon(":/images/2D.png"));
	this->meshMeshMaterial3d = new QTreeWidgetItem(this->meshMeshMaterial);
	this->meshMeshMaterial3d->setText(0, "3D");
	this->meshMeshMaterial3d->setIcon(0, QIcon(":/images/3D.png"));
	// mesh - selection
	this->meshMeshSelection = new QTreeWidgetItem(meshMesh);
	this->meshMeshSelection->setText(0, "Selection");
	// mesh - gradient
	this->meshMeshGradient = new QTreeWidgetItem(meshMesh);
	this->meshMeshGradient->setText(0, "Gradient");
	//	mesh - tetgen quality options
	this->meshMeshTetgen = new QTreeWidgetItem(meshMesh);
	this->meshMeshTetgen->setText(0, "Tetgen");
	// mesh - execute
	this->meshMeshExecute = new QTreeWidgetItem(meshMesh);
	this->meshMeshExecute->setFont(0, boldFont);
	this->meshMeshExecute->setText(0, "> Execute Mesh <");

	connect(meshTree, SIGNAL(itemClicked(QTreeWidgetItem *, int)), this, SLOT(settingCallByTree(QTreeWidgetItem *, int)));

	this->meshVBox->addWidget(this->meshTree);

	// adding tetgen options group box
	this->tetgenGBox = new QGroupBox(tr("Tetgen"), this->meshWidget);
	this->tetgenGBox->setHidden(true);
	this->tetgenGrid = new QGridLayout(this->tetgenGBox);
	this->tetgenLineEdit = new QLineEdit("pq1.2AY", this->tetgenGBox);
	this->tetgenGrid->addWidget(this->tetgenLineEdit, 0, 0, 1, 1);
	this->tetgenGBox->setLayout(this->tetgenGrid);
	this->meshVBox->addWidget(this->tetgenGBox);
	// adding refinement options group box
	refinementGBox = new QGroupBox(tr("Refinement"), meshWidget);
	refinementGBox->setHidden(true);
	refinementGrid = new QGridLayout(refinementGBox);
	refinementTableWidget = new QTableWidget(0, 2, refinementGBox);
	refinementTableWidget->verticalHeader()->hide();
	refinementTableWidget->horizontalHeader()->setStretchLastSection(true);
	refinementTableWidget->setHorizontalHeaderItem(0, new QTableWidgetItem("Object"));
	refinementTableWidget->setHorizontalHeaderItem(1, new QTableWidgetItem("Refinement"));
	refinementGrid->addWidget(refinementTableWidget, 0, 0, 1, 1);
	refinementGBox->setLayout(refinementGrid);
	meshVBox->addWidget(refinementGBox);
	// adding gradient change group box
	meshGradientGBox = new QGroupBox(tr("Radius"), meshWidget);
	meshGradientGBox->setHidden(true);
	meshGradientGrid = new QGridLayout(meshGradientGBox);
	meshGradientValue = new QDoubleSpinBox();
	meshGradientValue->setMaximum(9999.99);
	meshGradientValue->setMinimum(0);
	connect(meshGradientValue, SIGNAL(valueChanged(double)),
					this, SLOT(meshGradientUpdate(double)));
	meshGradientValue->setValue(Model.meshGradient);
	meshGradientValue->setSingleStep(0.25);
	meshGradientGrid->addWidget(meshGradientValue, 0, 0, 1, 1);
	meshGradientGBox->setLayout(meshGradientGrid);
	meshVBox->addWidget(meshGradientGBox);
	// adding interpolation options group box
	interpolationGBox = new QGroupBox(tr("Interpolation"), meshWidget);
	interpolationGBox->setHidden(true);
	interpolationGrid = new QGridLayout(interpolationGBox);
	interpolationMethod = new QComboBox(interpolationGBox);
	interpolationMethod->addItem("IDW");
	interpolationMethod->addItem("SPLINE");
	interpolationMethod->addItem("KRIGING");
	interpolationMethod->setCurrentIndex(0);
	connect(interpolationMethod, SIGNAL(currentIndexChanged(QString)), this, SLOT(interpolationSetMethod(QString)));
	interpolationGrid->addWidget(interpolationMethod, 0, 0, 1, 1);
	interpolationGBox->setLayout(interpolationGrid);
	meshVBox->addWidget(interpolationGBox);

	/* Gradient! */
	preMeshGradientGBox = new QGroupBox(tr("Radius"), meshWidget);
	preMeshGradientGBox->setHidden(true);
	preMeshGradientGrid = new QGridLayout(preMeshGradientGBox);
	preMeshGradientValue = new QDoubleSpinBox();
	preMeshGradientValue->setMaximum(9999.99);
	preMeshGradientValue->setMinimum(0);
	connect(preMeshGradientValue, SIGNAL(valueChanged(double)),
					this, SLOT(preMeshGradientUpdate(double)));
	preMeshGradientValue->setValue(Model.preMeshGradient);
	preMeshGradientValue->setSingleStep(0.25);
	preMeshGradientGrid->addWidget(preMeshGradientValue, 0, 0, 1, 1);
	preMeshGradientGBox->setLayout(preMeshGradientGrid);
	meshVBox->addWidget(preMeshGradientGBox);

	// adding selection options group box
	this->selectionGBox = new QGroupBox(tr("Selection"), this->meshWidget);
	this->selectionGBox->setHidden(true);
	this->selectionGrid = new QGridLayout(this->selectionGBox);
	this->selectionComponent = new QComboBox(this->selectionGBox);
	connect(this->selectionComponent, SIGNAL(currentIndexChanged(QString)), this, SLOT(setSShowGBox(QString)));
	connect(this->selectionComponent, SIGNAL(currentIndexChanged(QString)),
					this, SLOT(onHighlightingChange()));
	this->selectionGrid->addWidget(this->selectionComponent, 0, 0, 1, 4);

	this->selectionSingle = new QLabel(tr("Single"));
	this->selectionSingle->setAlignment(Qt::AlignHCenter);
	this->selectionGrid->addWidget(this->selectionSingle, 1, 1, 1, 1);

	this->selectionMulti = new QLabel(tr("Multi"));
	this->selectionMulti->setAlignment(Qt::AlignHCenter);
	this->selectionGrid->addWidget(this->selectionMulti, 1, 2, 1, 1);

	this->selectionAll = new QLabel(tr("All"));
	this->selectionAll->setAlignment(Qt::AlignHCenter);
	this->selectionGrid->addWidget(this->selectionAll, 1, 3, 1, 1);

	this->selectionMark = new QLabel(tr("Mark"));
	this->selectionGrid->addWidget(this->selectionMark, 2, 0, 1, 1);

	this->selectionMarkSingle = new QPushButton(QIcon(":/images/markSingle.png"), "", this->selectionGBox);
	this->selectionMarkSingle->setCheckable(true);
	this->selectionMarkSingle->setProperty("type", SINGLE);
	this->selectionMarkSingle->setProperty("mode", MARK);
	connect(this->selectionMarkSingle, SIGNAL(clicked(bool)), this, SLOT(selectionSetFlags(bool)));
	this->selectionGrid->addWidget(this->selectionMarkSingle, 2, 1, 1, 1);
	selectionButtons.push_back(this->selectionMarkSingle);

	this->selectionMarkMulti = new QPushButton(QIcon(":/images/markMulti.png"), "", this->selectionGBox);
	this->selectionMarkMulti->setCheckable(true);
	this->selectionMarkMulti->setProperty("type", MULTI);
	this->selectionMarkMulti->setProperty("mode", MARK);
	connect(this->selectionMarkMulti, SIGNAL(clicked(bool)), this, SLOT(selectionSetFlags(bool)));
	this->selectionGrid->addWidget(this->selectionMarkMulti, 2, 2, 1, 1);
	selectionButtons.push_back(this->selectionMarkMulti);

	this->selectionMarkAll = new QPushButton(QIcon(":/images/markAll.png"), "", this->selectionGBox);
	this->selectionMarkAll->setCheckable(true);
	this->selectionMarkAll->setProperty("type", ALL);
	this->selectionMarkAll->setProperty("mode", MARK);
	connect(this->selectionMarkAll, SIGNAL(clicked(bool)), this, SLOT(selectionSetFlags(bool)));
	this->selectionGrid->addWidget(this->selectionMarkAll, 2, 3, 1, 1);
	selectionButtons.push_back(this->selectionMarkAll);

	this->selectionUnmark = new QLabel(tr("Unmark"));
	this->selectionGrid->addWidget(this->selectionUnmark, 3, 0, 1, 1);

	this->selectionUnmarkSingle = new QPushButton(QIcon(":/images/unmarkSingle.png"), "", this->selectionGBox);
	this->selectionUnmarkSingle->setCheckable(true);
	this->selectionUnmarkSingle->setProperty("type", SINGLE);
	this->selectionUnmarkSingle->setProperty("mode", UNMARK);
	connect(this->selectionUnmarkSingle, SIGNAL(clicked(bool)), this, SLOT(selectionSetFlags(bool)));
	this->selectionGrid->addWidget(this->selectionUnmarkSingle, 3, 1, 1, 1);
	selectionButtons.push_back(this->selectionUnmarkSingle);

	this->selectionUnmarkMulti = new QPushButton(QIcon(":/images/unmarkMulti.png"), "", this->selectionGBox);
	this->selectionUnmarkMulti->setCheckable(true);
	this->selectionUnmarkMulti->setProperty("type", MULTI);
	this->selectionUnmarkMulti->setProperty("mode", UNMARK);
	connect(this->selectionUnmarkMulti, SIGNAL(clicked(bool)), this, SLOT(selectionSetFlags(bool)));
	this->selectionGrid->addWidget(this->selectionUnmarkMulti, 3, 2, 1, 1);
	selectionButtons.push_back(this->selectionUnmarkMulti);

	this->selectionUnmarkAll = new QPushButton(QIcon(":/images/unmarkAll.png"), "", this->selectionGBox);
	this->selectionUnmarkAll->setCheckable(true);
	this->selectionUnmarkAll->setProperty("type", ALL);
	this->selectionUnmarkAll->setProperty("mode", UNMARK);
	connect(this->selectionUnmarkAll, SIGNAL(clicked(bool)), this, SLOT(selectionSetFlags(bool)));
	this->selectionGrid->addWidget(this->selectionUnmarkAll, 3, 3, 1, 1);
	selectionButtons.push_back(this->selectionUnmarkAll);

	this->selectionHole = new QLabel(tr("Hole"));
	this->selectionGrid->addWidget(this->selectionHole, 4, 0, 1, 1);
	this->selectionHoleSingle = new QPushButton(QIcon(":/images/holeSingle.png"), "", this->selectionGBox);
	this->selectionHoleSingle->setCheckable(true);
	this->selectionHoleSingle->setProperty("type", SINGLE);
	this->selectionHoleSingle->setProperty("mode", HOLE);
	connect(this->selectionHoleSingle, SIGNAL(clicked(bool)), this, SLOT(selectionSetFlags(bool)));
	this->selectionGrid->addWidget(this->selectionHoleSingle, 4, 1, 1, 1);
	selectionButtons.push_back(this->selectionHoleSingle);

	this->selectionHoleMulti = new QPushButton(QIcon(":/images/holeMulti.png"), "", this->selectionGBox);
	this->selectionHoleMulti->setCheckable(true);
	this->selectionHoleMulti->setProperty("type", MULTI);
	this->selectionHoleMulti->setProperty("mode", HOLE);
	connect(this->selectionHoleMulti, SIGNAL(clicked(bool)), this, SLOT(selectionSetFlags(bool)));
	this->selectionGrid->addWidget(this->selectionHoleMulti, 4, 2, 1, 1);
	selectionButtons.push_back(this->selectionHoleMulti);

	this->selectionHoleAll = new QPushButton(QIcon(":/images/holeAll.png"), "", this->selectionGBox);
	this->selectionHoleAll->setCheckable(true);
	this->selectionHoleAll->setProperty("type", ALL);
	this->selectionHoleAll->setProperty("mode", HOLE);
	connect(this->selectionHoleAll, SIGNAL(clicked(bool)), this, SLOT(selectionSetFlags(bool)));
	this->selectionGrid->addWidget(this->selectionHoleAll, 4, 3, 1, 1);
	selectionButtons.push_back(this->selectionHoleAll);

	QLabel *selectionPaintBucket = new QLabel(tr("Bucket"));
	selectionPaintBucket->setAlignment(Qt::AlignHCenter);
	this->selectionGrid->addWidget(selectionPaintBucket, 6, 1, 1, 1);

	QLabel *selectionPolygon = new QLabel(tr("Polygon"));
	selectionPolygon->setAlignment(Qt::AlignHCenter);
	this->selectionGrid->addWidget(selectionPolygon, 6, 2, 1, 1);

	this->selectionGrid->addWidget(new QLabel(tr("Mark")), 7, 0, 1, 1);

	this->selectionMarkBucket = new QPushButton(QIcon(":/images/markBucket.png"), "", this->selectionGBox);
	this->selectionMarkBucket->setCheckable(true);
	this->selectionMarkBucket->setProperty("type", BUCKET);
	this->selectionMarkBucket->setProperty("mode", MARK);
	this->selectionGrid->addWidget(this->selectionMarkBucket, 7, 1, 1, 1);
	connect(this->selectionMarkBucket, SIGNAL(clicked(bool)), this, SLOT(selectionSetFlags(bool)));
	selectionButtons.push_back(this->selectionMarkBucket);

	this->selectionMarkPolygon = new QPushButton(QIcon(":/images/markPolygon.png"), "", this->selectionGBox);
	this->selectionMarkPolygon->setCheckable(true);
	this->selectionMarkPolygon->setProperty("type", POLYGON);
	this->selectionMarkPolygon->setProperty("mode", MARK);
	this->selectionGrid->addWidget(this->selectionMarkPolygon, 7, 2, 1, 1);
	connect(this->selectionMarkPolygon, SIGNAL(clicked(bool)), this, SLOT(selectionSetFlags(bool)));
	selectionButtons.push_back(this->selectionMarkPolygon);

	this->selectionGrid->addWidget(new QLabel(tr("Unmark")), 8, 0, 1, 1);

	this->selectionUnmarkBucket = new QPushButton(QIcon(":/images/unmarkBucket.png"), "", this->selectionGBox);
	this->selectionUnmarkBucket->setCheckable(true);
	this->selectionUnmarkBucket->setProperty("type", BUCKET);
	this->selectionUnmarkBucket->setProperty("mode", UNMARK);
	this->selectionGrid->addWidget(this->selectionUnmarkBucket, 8, 1, 1, 1);
	connect(this->selectionUnmarkBucket, SIGNAL(clicked(bool)), this, SLOT(selectionSetFlags(bool)));
	selectionButtons.push_back(this->selectionUnmarkBucket);

	this->selectionUnmarkPolygon = new QPushButton(QIcon(":/images/unmarkPolygon.png"), "", this->selectionGBox);
	this->selectionUnmarkPolygon->setCheckable(true);
	this->selectionUnmarkPolygon->setProperty("type", POLYGON);
	this->selectionUnmarkPolygon->setProperty("mode", UNMARK);
	this->selectionGrid->addWidget(this->selectionUnmarkPolygon, 8, 2, 1, 1);
	connect(this->selectionUnmarkPolygon, SIGNAL(clicked(bool)), this, SLOT(selectionSetFlags(bool)));
	selectionButtons.push_back(this->selectionUnmarkPolygon);

	this->selectionGrid->addWidget(new QLabel(tr("Hole")), 9, 0, 1, 1);

	this->selectionHoleBucket = new QPushButton(QIcon(":/images/holeBucket.png"), "", this->selectionGBox);
	this->selectionHoleBucket->setCheckable(true);
	this->selectionHoleBucket->setProperty("type", BUCKET);
	this->selectionHoleBucket->setProperty("mode", HOLE);
	this->selectionGrid->addWidget(this->selectionHoleBucket, 9, 1, 1, 1);
	connect(this->selectionHoleBucket, SIGNAL(clicked(bool)), this, SLOT(selectionSetFlags(bool)));
	selectionButtons.push_back(this->selectionHoleBucket);

	this->selectionHolePolygon = new QPushButton(QIcon(":/images/holePolygon.png"), "", this->selectionGBox);
	this->selectionHolePolygon->setCheckable(true);
	this->selectionHolePolygon->setProperty("type", POLYGON);
	this->selectionHolePolygon->setProperty("mode", HOLE);
	this->selectionGrid->addWidget(this->selectionHolePolygon, 9, 2, 1, 1);
	connect(this->selectionHolePolygon, SIGNAL(clicked(bool)), this, SLOT(selectionSetFlags(bool)));
	selectionButtons.push_back(this->selectionHolePolygon);

	this->selectionInvert = new QPushButton("Invert Selection", this->selectionGBox);
	this->selectionInvert->setCheckable(true);
	this->selectionInvert->setProperty("type", INVERT_TOOL);
	this->selectionInvert->setProperty("mode", INVERT);
	connect(this->selectionInvert, SIGNAL(clicked(bool)), this, SLOT(selectionSetFlags(bool)));
	this->selectionGrid->addWidget(this->selectionInvert, 10, 0, 1, 4);
	selectionButtons.push_back(this->selectionInvert);

	this->selectionMaterial = new AnimatedButton(
			"Material Selection", ":/images/loading.gif", selectionGBox);
	connect(this->selectionMaterial, SIGNAL(clicked(bool)),
					this, SLOT(applyMaterialSelection()));
	this->selectionGrid->addWidget(this->selectionMaterial, 11, 0, 1, 4);

	/* Add dropdown menu used to control the highlighting mode. */
	highlightingCombo = new QComboBox(selectionGBox);
	selectionGrid->addWidget(highlightingCombo, 12, 0, 1, 4);
	highlightingCombo->addItem("No highlighting",
														 int(HighlightingMode::NONE));
	highlightingCombo->addItem("Highlight unmarked",
														 int(HighlightingMode::COMPONENT));
	highlightingCombo->addItem("Highlight all unmarked",
														 int(HighlightingMode::ALL));
	highlightingCombo->addItem("Highlight individual",
														 int(HighlightingMode::SINGLE));

	connect(highlightingCombo, SIGNAL(currentIndexChanged(int)),
					this, SLOT(onHighlightingChange()));

	/* Add table used to list missing selections. */
	selectionTable = new QTableWidget(0, 3, selectionGBox);
	selectionGrid->addWidget(selectionTable, 13, 0, 1, 4);
	connect(selectionTable, SIGNAL(itemActivated(QTableWidgetItem *)),
					this, SLOT(tableClicked(QTableWidgetItem *)));

	selectionTable->verticalHeader()->hide();
	selectionTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
	selectionTable->setSelectionBehavior(QAbstractItemView::SelectRows);
	selectionTable->setSelectionMode(QAbstractItemView::SingleSelection);
	selectionTable->hide();

	this->selectionGBox->setLayout(this->selectionGrid);
	this->meshVBox->addWidget(this->selectionGBox);

	// adding material options 1D group box
	this->material1dGBox = new QGroupBox(tr("Material 1D"), meshWidget);
	this->material1dGBox->setHidden(true);
	this->material1dGrid = new QGridLayout(material1dGBox);
	this->material1dList = new QListWidget(material1dGBox);
	this->material1dGrid->addWidget(material1dList, 0, 0, 1, 1);
	this->material1dGBox->setLayout(material1dGrid);
	this->meshVBox->addWidget(material1dGBox);
	// adding material options 2D group box
	this->material2dGBox = new QGroupBox(tr("Material 2D"), meshWidget);
	this->material2dGBox->setHidden(true);
	this->material2dGrid = new QGridLayout(material2dGBox);
	this->material2dList = new QListWidget(material2dGBox);
	this->material2dGrid->addWidget(material2dList, 0, 0, 1, 1);
	this->material2dGBox->setLayout(material2dGrid);
	this->meshVBox->addWidget(material2dGBox);

	// adding material options 3D group box
	this->material3dGBox = new QGroupBox(tr("Material 3D"), meshWidget);
	this->material3dGBox->setHidden(true);
	this->material3dGrid = new QGridLayout(this->material3dGBox);
	this->material3dList = new QListWidget(this->material3dGBox);
	connect(this->material3dList, SIGNAL(currentRowChanged(int)), this, SLOT(material3dFillLocation(int)));
	this->material3dGrid->addWidget(this->material3dList, 0, 0, 1, 4);
	this->material3dAddRemoveMaterialLabel = new QLabel(tr("Add/Remove Material"));
	this->material3dGrid->addWidget(this->material3dAddRemoveMaterialLabel, 1, 0, 1, 2);
	this->material3dAddMaterialPushButton = new QPushButton("+", this->material3dGBox);
	connect(this->material3dAddMaterialPushButton, SIGNAL(clicked(bool)), this, SLOT(material3dAdd()));
	this->material3dGrid->addWidget(this->material3dAddMaterialPushButton, 1, 2, 1, 1);
	this->material3dRemoveMaterialPushButton = new QPushButton("-", this->material3dGBox);
	connect(this->material3dRemoveMaterialPushButton, SIGNAL(clicked(bool)), this, SLOT(material3dRemove()));
	this->material3dGrid->addWidget(this->material3dRemoveMaterialPushButton, 1, 3, 1, 1);
	this->material3dListLocation = new QListWidget(this->material3dGBox);
	connect(this->material3dListLocation, SIGNAL(currentRowChanged(int)), this, SLOT(material3dFillValue(int)));
	this->material3dGrid->addWidget(this->material3dListLocation, 2, 0, 1, 4);
	this->material3dAddRemoveLocationLabel = new QLabel(tr("Add/Remove Location"));
	this->material3dGrid->addWidget(this->material3dAddRemoveLocationLabel, 3, 0, 1, 2);
	this->material3dAddLocationPushButton = new QPushButton("+", this->material3dGBox);
	connect(this->material3dAddLocationPushButton, SIGNAL(clicked(bool)), this, SLOT(material3dLocationAdd()));
	this->material3dGrid->addWidget(this->material3dAddLocationPushButton, 3, 2, 1, 1);
	this->material3dRemoveLocationPushButton = new QPushButton("-", this->material3dGBox);
	connect(this->material3dRemoveLocationPushButton, SIGNAL(clicked(bool)), this, SLOT(material3dLocationRemove()));
	this->material3dGrid->addWidget(this->material3dRemoveLocationPushButton, 3, 3, 1, 1);
	this->material3dXLabel = new QLabel(tr("X:"));
	this->material3dGrid->addWidget(this->material3dXLabel, 4, 0, 1, 2);
	this->material3dXValue = new QDoubleSpinBox();
	this->material3dXValue->setProperty("location", "X");
	connect(this->material3dXValue, SIGNAL(valueChanged(double)), this, SLOT(material3dSetLocationFromDSpinBox(double)));
	this->material3dGrid->addWidget(this->material3dXValue, 4, 2, 1, 2);
	this->material3dYLabel = new QLabel(tr("Y:"));
	this->material3dGrid->addWidget(this->material3dYLabel, 6, 0, 1, 2);
	this->material3dYValue = new QDoubleSpinBox();
	this->material3dYValue->setProperty("location", "Y");
	connect(this->material3dYValue, SIGNAL(valueChanged(double)), this, SLOT(material3dSetLocationFromDSpinBox(double)));
	this->material3dGrid->addWidget(this->material3dYValue, 6, 2, 1, 2);
	this->material3dZLabel = new QLabel(tr("Z:"));
	this->material3dGrid->addWidget(this->material3dZLabel, 8, 0, 1, 2);
	this->material3dZValue = new QDoubleSpinBox();
	this->material3dZValue->setProperty("location", "Z");
	connect(this->material3dZValue, SIGNAL(valueChanged(double)), this, SLOT(material3dSetLocationFromDSpinBox(double)));
	this->material3dGrid->addWidget(this->material3dZValue, 8, 2, 1, 2);
	this->material3dXSlider = new QSlider(Qt::Horizontal);
	this->material3dXSlider->setProperty("location", "X");
	connect(this->material3dXSlider, SIGNAL(sliderMoved(int)), this, SLOT(material3dSetLocationFromSlider(int)));
	this->material3dGrid->addWidget(this->material3dXSlider, 5, 0, 1, 4);
	this->material3dYSlider = new QSlider(Qt::Horizontal);
	this->material3dYSlider->setProperty("location", "Y");
	connect(this->material3dYSlider, SIGNAL(sliderMoved(int)), this, SLOT(material3dSetLocationFromSlider(int)));
	this->material3dGrid->addWidget(this->material3dYSlider, 7, 0, 1, 4);
	this->material3dZSlider = new QSlider(Qt::Horizontal);
	this->material3dZSlider->setProperty("location", "Z");
	connect(this->material3dZSlider, SIGNAL(sliderMoved(int)), this, SLOT(material3dSetLocationFromSlider(int)));
	this->material3dGrid->addWidget(this->material3dZSlider, 9, 0, 1, 4);
	this->material3dGBox->setLayout(material3dGrid);
	meshVBox->addWidget(material3dGBox);

	/* Error widget. */
	errorWidget = new QTreeWidgetItem(meshTree);
	errorWidget->setText(0, "Errors");
	errorWidget->setHidden(true);

	errorGBox = new QGroupBox(tr("Errors"), meshWidget);
	errorLayout = new QVBoxLayout(errorGBox);

	errorGBox->setHidden(true);
	errorGBox->setLayout(errorLayout);
	meshVBox->addWidget(errorGBox);

	errorLabel = new QLabel(errorGBox);
	errorLayout->addWidget(errorLabel);
	errorLayout->addWidget(new QLabel("List of involved triangles:",
																		errorGBox));

	/* Add table used to list self intersections. */
	errorTable = new QTableWidget(0, 1, errorGBox);
	errorLayout->addWidget(errorTable);
	connect(errorTable, SIGNAL(itemActivated(QTableWidgetItem *)),
					this, SLOT(errorTableClicked(QTableWidgetItem *)));

	errorTable->verticalHeader()->hide();
	errorTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
	errorTable->setSelectionBehavior(QAbstractItemView::SelectRows);
	errorTable->setSelectionMode(QAbstractItemView::SingleSelection);

	connect(&Model, SIGNAL(ErrorInfoChanged(QString)),
					this, SLOT(updateErrorInfo(QString)));

	QPushButton *errorClearButton = new QPushButton("Clear markers", errorGBox);
	errorLayout->addWidget(errorClearButton);

	connect(errorClearButton, SIGNAL(clicked(bool)),
					this, SLOT(clearErrorMarkers()));

	meshSpacer = new QSpacerItem(0, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
	meshVBox->addItem(meshSpacer);
	meshWidget->setLayout(meshVBox);
	meshDock->setWidget(meshWidget);
	addDockWidget(Qt::LeftDockWidgetArea, meshDock);
	viewMenu->addAction(meshDock->toggleViewAction());
}

void MainWindow::readSettings()
{
	QSettings settings("GFZ", "MeshIT");
	restoreGeometry(settings.value("geometry").toByteArray());
	restoreState(settings.value("windowState").toByteArray());
}

void MainWindow::writeSettings()
{
	QSettings settings("GFZ", "MeshIT");
	settings.setValue("geometry", saveGeometry());
	settings.setValue("windowState", saveState());
}

void MainWindow::preMeshJob()
{
	clearErrorInfo();

	int currentStep, totalSteps;
	//	start
	startdate = QDateTime::currentDateTime();
	emit progress_append(">Start Time: " + startdate.toString() + "\n");
	//	convex hull
	emit progress_append(">Start calculating convexhull...\n");
	currentStep = 0;
	totalSteps = Model.Surfaces.length();
	for (int s = 0; s != Model.Surfaces.length(); s++)
	{
		C_Task *task = new C_Task(this, "CONVEXHULL", s, 0, ++currentStep, totalSteps);
		QThreadPool::globalInstance()->start(task);
	}
	QThreadPool::globalInstance()->waitForDone();
	emit progress_append(">...finished");
	//	segmentation coarse
	emit progress_append(">Start coarse segmentation...\n");
	currentStep = 0;
	totalSteps = Model.Polylines.length();
	for (int p = 0; p != Model.Polylines.length(); p++)
	{
		C_Task *task = new C_Task(this, "SEGMENTS", p, 0, ++currentStep, totalSteps);
		QThreadPool::globalInstance()->start(task);
	}
	QThreadPool::globalInstance()->waitForDone();
	emit progress_append(">...finished");
	//	2D triangulation coarse
	emit progress_append(">Start coarse triangulation...\n");
	currentStep = 0;
	totalSteps = Model.Surfaces.length();
	for (int s = 0; s != Model.Surfaces.length(); s++)
	{
		C_Task *task = new C_Task(this, "TRIANGLES", s, 0, ++currentStep, totalSteps);
		QThreadPool::globalInstance()->start(task);
	}
	QThreadPool::globalInstance()->waitForDone();
	emit progress_append(">...finished");
	//	intersection: surface-surface
	emit progress_append(">Start calculating surface-surface intersections...\n");
	Model.Intersections.clear();
	currentStep = 0;
	totalSteps = Model.Surfaces.length() * (Model.Surfaces.length() - 1) / 2;
	if (totalSteps > 0)
	{
		for (int s1 = 0; s1 != Model.Surfaces.length() - 1; s1++)
			for (int s2 = s1 + 1; s2 != Model.Surfaces.length(); s2++)
			{
				C_Task *task = new C_Task(this, "INTERSECTION_MESH_MESH", s1, s2, ++currentStep, totalSteps);
				QThreadPool::globalInstance()->start(task);
			}
		QThreadPool::globalInstance()->waitForDone();
	}
	emit progress_append(">...finished");
	//	intersection: polyline-surface
	emit progress_append(">Start calculating polyline-surface intersections...\n");
	currentStep = 0;
	totalSteps = Model.Polylines.length() * Model.Surfaces.length();
	if (totalSteps > 0)
	{
		for (int p = 0; p != Model.Polylines.length(); p++)
			for (int s = 0; s != Model.Surfaces.length(); s++)
			{
				C_Task *task = new C_Task(this, "INTERSECTION_POLYLINE_MESH", p, s, ++currentStep, totalSteps);
				QThreadPool::globalInstance()->start(task);
			}
		QThreadPool::globalInstance()->waitForDone();
	}
	emit progress_append(">...finished");
	//	intersection: calculate size
	Model.calculate_size_of_intersections();
	//	intersection: triple points
	emit progress_append(">Start calculating intersection triplepoints...\n");
	Model.TPs.clear();
	currentStep = 0;
	totalSteps = Model.Intersections.length() * (Model.Intersections.length() - 1) / 2;
	if (totalSteps > 0)
	{
		for (int i1 = 0; i1 != Model.Intersections.length() - 1; i1++)
			for (int i2 = i1 + 1; i2 != Model.Intersections.length(); i2++)
			{
				C_Task *task = new C_Task(this, "INTERSECTION_TRIPLEPOINTS", i1, i2, ++currentStep, totalSteps);
				QThreadPool::globalInstance()->start(task);
			}
		QThreadPool::globalInstance()->waitForDone();
	}
	Model.insert_int_triplepoints();
	emit progress_append(">...finished");

	emit progress_append(">Start aligning Convex Hulls to Intersections...\n");
	for (int s = 0; s != Model.Surfaces.length(); s++)
	{
		emit progress_replace("   >(" + QString::number(s + 1) + "/" + QString::number(Model.Surfaces.length()) + ") " + Model.Surfaces[s].Name + " (" + Model.Surfaces[s].Type + ")");
		Model.Surfaces[s].alignIntersectionsToConvexHull();
	}
	emit progress_append(">...finished");
	//	model constraints
	emit progress_append(">Start calculating constraints...\n");
	for (int s = 0; s != Model.Surfaces.length(); s++)
		Model.Surfaces[s].calculate_Constraints();
	for (int p = 0; p != Model.Polylines.length(); p++)
		Model.Polylines[p].calculate_Constraints();
	emit progress_append(">...finished");
	Model.calculate_size_of_constraints();
	//	end
	enddate = QDateTime::currentDateTime();
	emit progress_append(">End Time: " + enddate.toString() + "\n");
	emit progress_append(">elapsed Time: " + QString::number(startdate.msecsTo(enddate)) + "\n");
}

void MainWindow::MeshJob()
{
	clearErrorInfo();

	int currentStep, totalSteps;
	//	start
	startdate = QDateTime::currentDateTime();
	emit progress_append(">Start Time: " + startdate.toString() + "\n");

	//	segmentation fine
	emit progress_append(">Start fine segmentation...\n");
	currentStep = 0;
	totalSteps = Model.Polylines.length();
	for (int p = 0; p != Model.Polylines.length(); p++)
	{
		C_Task *task = new C_Task(this, "SEGMENTS_FINE", p, 0, ++currentStep, totalSteps);
		QThreadPool::globalInstance()->start(task);
	}
	QThreadPool::globalInstance()->waitForDone();
	emit progress_append(">...finished");
	//	2D triangulation fine
	emit progress_append(">Start fine triangulation...\n");
	currentStep = 0;
	totalSteps = Model.Surfaces.length();
	for (int s = 0; s != Model.Surfaces.length(); s++)
	{
		C_Task *task = new C_Task(this, "TRIANGLES_FINE", s, 0, ++currentStep, totalSteps);
		QThreadPool::globalInstance()->start(task);
	}
	QThreadPool::globalInstance()->waitForDone();
	emit progress_append(">...finished");
	//	3D tetrahedralization
	emit progress_append(">Start tetrahedralization...");
	Model.calculate_tets(this->tetgenLineEdit->text());
	emit progress_append(">...finished");

	emit progress_append(">Start verification...");
	Model.verify_materials();
	emit progress_append(">...finished");

	//	end
	enddate = QDateTime::currentDateTime();
	emit progress_append(">End Time: " + enddate.toString() + "\n");
	emit progress_append(">elapsed Time: " + QString::number(startdate.msecsTo(enddate)) + "\n");
}

/* Generate selections based on the assigned materials. This invokes tetgen a
 * first time (without refinement switch -q) to find adjacent constraints. */
void MainWindow::materialSelectionJob()
{
	clearErrorInfo();

	int currentStep, totalSteps;

	emit progress_append("> Start material selection\n");

	Model.select_all_constraints();

	currentStep = 0;
	totalSteps = Model.Polylines.length();
	for (int p = 0; p != Model.Polylines.length(); p++)
	{
		C_Task *task = new C_Task(this, "SEGMENTS_FINE", p, 0,
															++currentStep, totalSteps);
		QThreadPool::globalInstance()->start(task);
	}
	QThreadPool::globalInstance()->waitForDone();

	currentStep = 0;
	totalSteps = Model.Surfaces.length();
	for (int s = 0; s != Model.Surfaces.length(); s++)
	{
		C_Task *task = new C_Task(this, "TRIANGLES_FINE", s, 0,
															++currentStep, totalSteps);
		QThreadPool::globalInstance()->start(task);
	}
	QThreadPool::globalInstance()->waitForDone();

	/* Save previous model. */
	C_Mesh3D *prevModel = Model.Mesh;
	Model.Mesh = 0;

	/* Refinement switch -q is not required here since the result is only used
	 * temporarily to find selections based on the defined materials. */
	Model.calculate_tets("pAY");

	if (!Model.Mesh)
	{
		Model.deselect_all_constraints();
		emit Model.PrintError("> Material selection failed\n");
		return;
	}

	Model.material_selections();

	/* Resotre previous model. */
	delete Model.Mesh;
	Model.Mesh = prevModel;

	emit progress_append("> Material selection completed\n");
	emit Model.ModelInfoChanged();
}

void MainWindow::clearMesh()
{
	if (Model.Mesh)
	{
		delete Model.Mesh;
		Model.Mesh = 0;
		Model.listTets = NULL;
		Model.drawTets = false;
		this->matsTitle->setHidden(true);
	}
}

// ******************** //
//	Private slots		//
// ******************** //

void MainWindow::open()
{
	Model.FileNameTmp = QFileDialog::getOpenFileName(this, tr("Open MeshIt File"), Model.FilePath, tr("MeshIt File (*.pvd)"));
	QApplication::processEvents();
	if (!Model.FileNameTmp.isEmpty())
	{
		this->reset();
		Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
		Model.FileNameModel = Model.FileNameTmp;
		emit progress_append(">Start opening " + Model.FileNameTmp + "...");
		QApplication::processEvents();
		Model.Open();
		emit progress_append(">...finished");
	}
	this->FinishedRead();
}

void MainWindow::save()
{
	if (Model.FileNameModel.isEmpty())
	{
		Model.FileNameTmp = QFileDialog::getSaveFileName(this, tr("Save As MeshIt File"), Model.FilePath, tr("MeshIt File (*.pvd)"));
		if (!Model.FileNameTmp.isEmpty())
		{
			Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
			Model.FileNameModel = Model.FileNameTmp;
		}
	}
	if (!Model.FileNameModel.isEmpty())
	{
		emit progress_append(">Start saving " + Model.FileNameModel + "...");
		QApplication::processEvents();
		Model.Save();
		emit progress_append(">...finished");
	}
}

void MainWindow::saveAs()
{
	Model.FileNameTmp = QFileDialog::getSaveFileName(this, tr("Save As MeshIt File"), Model.FilePath, tr("MeshIt File (*.pvd)"));
	if (!Model.FileNameTmp.isEmpty())
	{
		Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
		Model.FileNameModel = Model.FileNameTmp;
		emit progress_append(">Start saving " + Model.FileNameTmp + "...");
		QApplication::processEvents();
		Model.Save();
		emit progress_append(">...finished");
	}
}

void MainWindow::saveScreen()
{
	Model.FileNameTmp = QFileDialog::getSaveFileName(this, tr("Save Screenshot"), Model.FilePath, tr("Portable Network Graphics (*.png);;JPEG (*.jpg);;Bitmap (*.bmp);;Tagged Image File Format (*.tif)"));
	if (!Model.FileNameTmp.isEmpty())
	{
		Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
		emit progress_append(">Start saving " + Model.FileNameTmp + "...");
		QApplication::processEvents();
		glWidget->saveScreen(Model.FileNameTmp);
		emit progress_append(">...finished");
	}
}

void MainWindow::importGoCad()
{
	Model.FileNameTmp = QFileDialog::getOpenFileName(this, tr("Open File"), Model.FilePath, tr("GoCad (*.gp)"));
	if (!Model.FileNameTmp.isEmpty())
	{
		Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
		this->reset();
		emit progress_append(">Start importing GOCAD " + Model.FileNameTmp + "...");
		QApplication::processEvents();
		Model.ReadGocadFile();
		emit progress_append(">...finished");
	}
	this->FinishedRead();
}

void MainWindow::exportVTU3D()
{
	Model.FileNameTmp = QFileDialog::getSaveFileName(this, tr("Export File"), Model.FilePath, tr("Paraview 3D Mesh (*.vtu)"));
	if (!Model.FileNameTmp.isEmpty())
	{
		Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
		emit progress_append(">Start exporting " + Model.FileNameTmp + "...");
		QApplication::processEvents();
		Model.ExportVTU3D();
		emit progress_append(">...finished");
	}
}

void MainWindow::exportVTU2D()
{
	QDialog *dialog = new QDialog(this);
	dialog->setModal(true);
	QGridLayout *layout = new QGridLayout(dialog);
	QListWidget *surfaces = new QListWidget(dialog);
	QListWidgetItem *item;
	QRadioButton *sep;
	QRadioButton *com;
	for (int s = 0; s != Model.Surfaces.length(); s++)
	{
		if (Model.Surfaces[s].Type == "UNIT")
			item = new QListWidgetItem(QIcon(":/images/units.png"), Model.Surfaces[s].Name);
		if (Model.Surfaces[s].Type == "FAULT")
			item = new QListWidgetItem(QIcon(":/images/faults.png"), Model.Surfaces[s].Name);
		if (Model.Surfaces[s].Type == "BORDER")
			item = new QListWidgetItem(QIcon(":/images/borders.png"), Model.Surfaces[s].Name);
		surfaces->addItem(item);
	}
	surfaces->setSelectionMode(QAbstractItemView::ExtendedSelection);
	layout->addWidget(surfaces, 0, 0, 1, 2);
	sep = new QRadioButton("separate");
	sep->setChecked(true);
	layout->addWidget(sep, 1, 0, 1, 1);
	com = new QRadioButton("combined");
	layout->addWidget(com, 1, 1, 1, 1);
	QHBoxLayout *buttonLayout = new QHBoxLayout;
	layout->addLayout(buttonLayout, 2, 0, 1, 2);
	buttonLayout->addStretch();
	QPushButton *cancelButton = new QPushButton(tr("Cancel"));
	connect(cancelButton, SIGNAL(clicked()), dialog, SLOT(reject()));
	buttonLayout->addWidget(cancelButton);
	QPushButton *okButton = new QPushButton(tr("Ok"));
	connect(okButton, SIGNAL(clicked()), dialog, SLOT(accept()));
	buttonLayout->addWidget(okButton);
	okButton->setDefault(true);
	int ret = dialog->exec();
	if (ret == QDialog::Rejected)
		return;
	if (sep->isChecked())
	{
		for (int m = 0; m != surfaces->count(); m++)
			if (surfaces->item(m)->isSelected())
			{
				Model.FileNameTmp = QFileDialog::getSaveFileName(this, tr("Export File"), Model.FilePath + "/" + surfaces->item(m)->text(), tr("Paraview 2D surface (*.vtu)"));
				if (!Model.FileNameTmp.isEmpty())
				{
					Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
					emit progress_append(">Start exporting " + Model.FileNameTmp + "...");
					QApplication::processEvents();
					Model.ExportVTU2D(QString::number(m));
					emit progress_append(">...finished");
				}
			}
	}
	else
	{
		QString meshIDs = "";
		for (int m = 0; m != surfaces->count(); m++)
		{
			if (surfaces->item(m)->isSelected())
				if (meshIDs.length() == 0)
					meshIDs += QString::number(m);
				else
					meshIDs += "," + QString::number(m);
		}
		Model.FileNameTmp = QFileDialog::getSaveFileName(this, tr("Export File"), Model.FilePath + "/all", tr("Paraview 2D surface (*.vtu)"));
		if (!Model.FileNameTmp.isEmpty())
		{
			Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
			emit progress_append(">Start exporting " + Model.FileNameTmp + "...");
			QApplication::processEvents();
			Model.ExportVTU2D(meshIDs);
			emit progress_append(">...finished");
		}
	}
}

void MainWindow::exportFeFlow()
{
	Model.FileNameTmp = QFileDialog::getSaveFileName(this, tr("Export File"), Model.FilePath, tr("ASCII Interchange FEM format (*.fem)"));
	if (!Model.FileNameTmp.isEmpty())
	{
		Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
		emit progress_append(">Start exporting " + Model.FileNameTmp + "...");
		QApplication::processEvents();
		Model.ExportFeFlow();
		emit progress_append(">...finished");
	}
}

void MainWindow::exportOGS()
{
	Model.FileNameTmp = QFileDialog::getSaveFileName(this, tr("Export File"), Model.FilePath, tr("OpenGeoSys Mesh (*.msh)"));
	if (!Model.FileNameTmp.isEmpty())
	{
		Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
		emit progress_append(">Start exporting " + Model.FileNameTmp + "...");
		QApplication::processEvents();
		Model.ExportOGS();
		emit progress_append(">...finished");
	}
}

void MainWindow::exportTIN()
{
	QDialog *dialog = new QDialog(this);
	dialog->setModal(true);
	QGridLayout *layout = new QGridLayout(dialog);
	QListWidget *surfaces = new QListWidget(dialog);
	QListWidgetItem *item;
	QRadioButton *sep;
	QRadioButton *com;
	QCheckBox *all;
	for (int s = 0; s != Model.Surfaces.length(); s++)
	{
		if (Model.Surfaces[s].Type == "UNIT")
			item = new QListWidgetItem(QIcon(":/images/units.png"), Model.Surfaces[s].Name);
		if (Model.Surfaces[s].Type == "FAULT")
			item = new QListWidgetItem(QIcon(":/images/faults.png"), Model.Surfaces[s].Name);
		if (Model.Surfaces[s].Type == "BORDER")
			item = new QListWidgetItem(QIcon(":/images/borders.png"), Model.Surfaces[s].Name);
		surfaces->addItem(item);
	}
	surfaces->setSelectionMode(QAbstractItemView::ExtendedSelection);
	layout->addWidget(surfaces, 0, 0, 1, 2);
	sep = new QRadioButton("separate");
	sep->setChecked(true);
	layout->addWidget(sep, 1, 0, 1, 1);
	com = new QRadioButton("combined");
	layout->addWidget(com, 1, 1, 1, 1);
	QHBoxLayout *buttonLayout = new QHBoxLayout;
	layout->addLayout(buttonLayout, 2, 0, 1, 2);
	buttonLayout->addStretch();
	QPushButton *cancelButton = new QPushButton(tr("Cancel"));
	connect(cancelButton, SIGNAL(clicked()), dialog, SLOT(reject()));
	buttonLayout->addWidget(cancelButton);
	QPushButton *okButton = new QPushButton(tr("Ok"));
	connect(okButton, SIGNAL(clicked()), dialog, SLOT(accept()));
	buttonLayout->addWidget(okButton);
	okButton->setDefault(true);
	int ret = dialog->exec();
	if (ret == QDialog::Rejected)
		return;
	if (sep->isChecked())
	{
		for (int m = 0; m != surfaces->count(); m++)
			if (surfaces->item(m)->isSelected())
			{
				Model.FileNameTmp = QFileDialog::getSaveFileName(this, tr("Export File"), Model.FilePath + "/" + surfaces->item(m)->text(), tr("TIN surface (*.tin)"));
				if (!Model.FileNameTmp.isEmpty())
				{
					Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
					emit progress_append(">Start exporting " + Model.FileNameTmp + "...");
					QApplication::processEvents();
					Model.ExportTIN(QString::number(m));
					emit progress_append(">...finished");
				}
			}
	}
	else
	{
		QString meshIDs = "";
		for (int m = 0; m != surfaces->count(); m++)
		{
			if (surfaces->item(m)->isSelected())
				if (meshIDs.length() == 0)
					meshIDs += QString::number(m);
				else
					meshIDs += "," + QString::number(m);
		}
		Model.FileNameTmp = QFileDialog::getSaveFileName(this, tr("Export File"), Model.FilePath + "/all", tr("TIN surface (*.tin)"));
		if (!Model.FileNameTmp.isEmpty())
		{
			Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
			emit progress_append(">Start exporting " + Model.FileNameTmp + "...");
			QApplication::processEvents();
			Model.ExportTIN(meshIDs);
			emit progress_append(">...finished");
		}
	}
}

void MainWindow::exportCOMSOL()
{
	Model.FileNameTmp = QFileDialog::getSaveFileName(this, tr("Export File"), Model.FilePath, tr("COMSOL Mesh (*.mphtxt)"));
	if (!Model.FileNameTmp.isEmpty())
	{
		Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
		emit progress_append(">Start exporting " + Model.FileNameTmp + "...");
		QApplication::processEvents();
		Model.ExportCOMSOL();
		emit progress_append(">...finished");
	}
}

void MainWindow::exportABAQUS()
{
	QDialog *dialog = new QDialog(this);
	dialog->setModal(true);
	QGridLayout *layout = new QGridLayout(dialog);
	QListWidget *surfaces = new QListWidget(dialog);
	QListWidgetItem *item;
	QRadioButton *sep;
	QRadioButton *com;
	for (int s = 0; s != Model.Surfaces.length(); s++)
	{
		if (Model.Surfaces[s].Type == "UNIT")
			item = new QListWidgetItem(QIcon(":/images/units.png"), Model.Surfaces[s].Name);
		if (Model.Surfaces[s].Type == "FAULT")
			item = new QListWidgetItem(QIcon(":/images/faults.png"), Model.Surfaces[s].Name);
		if (Model.Surfaces[s].Type == "BORDER")
			item = new QListWidgetItem(QIcon(":/images/borders.png"), Model.Surfaces[s].Name);
		item->setCheckState(Qt::Unchecked);
		surfaces->addItem(item);
	}
	layout->addWidget(surfaces, 0, 0, 1, 2);
	sep = new QRadioButton("separate");
	sep->setChecked(true);
	QList<int> borders;
	QHBoxLayout *buttonLayout = new QHBoxLayout;
	layout->addLayout(buttonLayout, 2, 0, 1, 2);
	buttonLayout->addStretch();
	QPushButton *cancelButton = new QPushButton(tr("Cancel"));
	connect(cancelButton, SIGNAL(clicked()), dialog, SLOT(reject()));
	buttonLayout->addWidget(cancelButton);
	QPushButton *okButton = new QPushButton(tr("Ok"));
	connect(okButton, SIGNAL(clicked()), dialog, SLOT(accept()));
	buttonLayout->addWidget(okButton);
	okButton->setDefault(true);
	int ret = dialog->exec();
	if (ret == QDialog::Rejected)
		return;
	QString meshIDs = "";
	if (sep->isChecked())
		for (int m = 0; m != surfaces->count(); m++)
		{
			if (surfaces->item(m)->checkState() == Qt::Checked)
				if (meshIDs.length() == 0)
					meshIDs += QString::number(m);
				else
					meshIDs += "," + QString::number(m);
		}
	Model.FileNameTmp = QFileDialog::getSaveFileName(this, tr("Export File"), Model.FilePath, tr("Abaqus 3D Mesh (*.inp)"));
	Model.FileNameTmp = QDir::toNativeSeparators(Model.FileNameTmp);
	if (!Model.FileNameTmp.isEmpty())
	{
		Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
		emit progress_append(">Start exporting " + Model.FileNameTmp + "...");
		QApplication::processEvents();
		Model.ExportABAQUS(meshIDs);
		emit progress_append(">...finished");
	}
}

void MainWindow::exportEXODUS()
{
#ifndef NOEXODUS
	QDialog *dialog = new QDialog(this);
	dialog->setModal(true);
	QGridLayout *layout = new QGridLayout(dialog);
	QListWidget *surfaces = new QListWidget(dialog);
	QListWidgetItem *item;
	QRadioButton *com;
	for (int s = 0; s != Model.Surfaces.length(); s++)
	{
		if (Model.Surfaces[s].Type == "UNIT")
			item = new QListWidgetItem(QIcon(":/images/units.png"), Model.Surfaces[s].Name);
		if (Model.Surfaces[s].Type == "FAULT")
			item = new QListWidgetItem(QIcon(":/images/faults.png"), Model.Surfaces[s].Name);
		if (Model.Surfaces[s].Type == "BORDER")
			item = new QListWidgetItem(QIcon(":/images/borders.png"), Model.Surfaces[s].Name);
		item->setCheckState(Qt::Unchecked);
		surfaces->addItem(item);
	}
	layout->addWidget(surfaces, 0, 0, 1, 2);
	QHBoxLayout *angleLayout = new QHBoxLayout;
	layout->addLayout(angleLayout, 1, 0, 1, 2);
	angleLayout->addStretch();
	QLabel *angleLabel = new QLabel("Rotate mesh counterclockwise by: ");
	angleLayout->addWidget(angleLabel);
	QDoubleSpinBox *angelValue = new QDoubleSpinBox();
	angelValue->setValue(Model.ExportRotationAngle);
	angelValue->setRange(0, 360);
	angelValue->setDecimals(1);
	connect(angelValue, SIGNAL(valueChanged(double)), this, SLOT(ExportRotationAngelUpdate(double)));
	angleLayout->addWidget(angelValue);
	QList<int> borders;
	QHBoxLayout *buttonLayout = new QHBoxLayout;
	layout->addLayout(buttonLayout, 2, 0, 1, 2);
	buttonLayout->addStretch();
	QPushButton *cancelButton = new QPushButton(tr("Cancel"));
	connect(cancelButton, SIGNAL(clicked()), dialog, SLOT(reject()));
	buttonLayout->addWidget(cancelButton);
	QPushButton *okButton = new QPushButton(tr("Ok"));
	connect(okButton, SIGNAL(clicked()), dialog, SLOT(accept()));
	buttonLayout->addWidget(okButton);
	okButton->setDefault(true);
	int ret = dialog->exec();
	if (ret == QDialog::Rejected)
		return;
	QString meshIDs = "";
	for (int m = 0; m != surfaces->count(); m++)
	{
		if (surfaces->item(m)->checkState() == Qt::Checked)
			if (meshIDs.length() == 0)
				meshIDs += QString::number(m);
			else
				meshIDs += "," + QString::number(m);
	}
	Model.FileNameTmp = QFileDialog::getSaveFileName(this, tr("Export File"), Model.FilePath, tr("Exodus 3D Mesh (*.e)"));
	Model.FileNameTmp = QDir::toNativeSeparators(Model.FileNameTmp);
	if (!Model.FileNameTmp.isEmpty())
	{
		Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
		emit progress_append(">Start exporting " + Model.FileNameTmp + "...");
		QApplication::processEvents();
		Model.ExportEXODUS(meshIDs);
		emit progress_append(">...finished");
	}
#endif
}

void MainWindow::addUnit()
{
	Model.FileNamesTmp = QFileDialog::getOpenFileNames(this, tr("Select one or more Unit Data files to open"), Model.FilePath, tr("All Supported Files (*.vtu *.dat *.txt *.csv);;VTK UnstructuredGrid Files (*.vtu);;XYZ Coordinate Files (*.dat *.txt *.csv *.xyz)"));
	if (!Model.FileNamesTmp.isEmpty())
		for (int f = 0; f != Model.FileNamesTmp.length(); f++)
		{
			Model.FileNameTmp = Model.FileNamesTmp[f];
			Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
			emit progress_append(">Start adding unit " + Model.FileNameTmp + "...");
			QApplication::processEvents();
			Model.AddSurface("UNIT");
			emit progress_append(">...finished");
		}
	this->FinishedRead();
}

void MainWindow::addFault()
{
	Model.FileNamesTmp = QFileDialog::getOpenFileNames(this, tr("Select one or more Fault Data files to open"), Model.FilePath, tr("All Supported Files (*.vtu *.dat *.txt *.csv);;VTK UnstructuredGrid Files (*.vtu);;XYZ Coordinate Files (*.dat *.txt *.csv)"));
	if (!Model.FileNamesTmp.isEmpty())
		for (int f = 0; f != Model.FileNamesTmp.length(); f++)
		{
			Model.FileNameTmp = Model.FileNamesTmp[f];
			Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
			emit progress_append(">Start adding fault " + Model.FileNameTmp + "...");
			QApplication::processEvents();
			Model.AddSurface("FAULT");
			emit progress_append(">...finished");
		}
	this->FinishedRead();
}

void MainWindow::addBorder()
{
	Model.FileNamesTmp = QFileDialog::getOpenFileNames(this, tr("Select one or more Border Data files to open"), Model.FilePath, tr("All Supported Files (*.vtu *.dat *.txt *.csv);;VTK UnstructuredGrid Files (*.vtu);;XYZ Coordinate Files (*.dat *.txt *.csv)"));
	if (!Model.FileNamesTmp.isEmpty())
		for (int f = 0; f != Model.FileNamesTmp.length(); f++)
		{
			Model.FileNameTmp = Model.FileNamesTmp[f];
			Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
			emit progress_append(">Start adding border " + Model.FileNameTmp + "...");
			QApplication::processEvents();
			Model.AddSurface("BORDER");
			emit progress_append(">...finished");
		}
	this->FinishedRead();
}

void MainWindow::addWell()
{
	Model.FileNamesTmp = QFileDialog::getOpenFileNames(this, tr("Select one or more Well Data files to open"), Model.FilePath, tr("All Supported Files (*.vtu *.dat *.txt *.csv);;VTK UnstructuredGrid Files (*.vtu);;XYZ Coordinate Files (*.dat *.txt *.csv)"));
	if (!Model.FileNamesTmp.isEmpty())
		for (int f = 0; f != Model.FileNamesTmp.length(); f++)
		{
			Model.FileNameTmp = Model.FileNamesTmp[f];
			Model.FilePath = Model.FileNameTmp.section("/", 0, -2);
			emit progress_append(">Start adding well " + Model.FileNameTmp + "...");
			QApplication::processEvents();
			Model.AddPolyline("WELL");
			emit progress_append(">...finished");
		}
	this->FinishedRead();
}

void MainWindow::deleteSurface()
{
	QDialog *dialog = new QDialog(this);
	dialog->setModal(true);
	QGridLayout *layout = new QGridLayout(dialog);
	QListWidget *surfaces_and_wells = new QListWidget(dialog);
	QListWidgetItem *item;
	for (int s = 0; s != Model.Surfaces.length(); s++)
	{
		if (Model.Surfaces[s].Type == "UNIT")
			item = new QListWidgetItem(QIcon(":/images/units.png"), Model.Surfaces[s].Name);
		if (Model.Surfaces[s].Type == "FAULT")
			item = new QListWidgetItem(QIcon(":/images/faults.png"), Model.Surfaces[s].Name);
		if (Model.Surfaces[s].Type == "BORDER")
			item = new QListWidgetItem(QIcon(":/images/borders.png"), Model.Surfaces[s].Name);
		item->setCheckState(Qt::Unchecked);
		surfaces_and_wells->addItem(item);
	}
	for (int p = 0; p != Model.Polylines.length(); p++)
	{
		if (Model.Polylines[p].Type == "WELL")
			item = new QListWidgetItem(QIcon(":/images/wells.png"), Model.Polylines[p].Name);
		item->setCheckState(Qt::Unchecked);
		surfaces_and_wells->addItem(item);
	}
	layout->addWidget(surfaces_and_wells, 0, 0, 1, 2);
	QHBoxLayout *buttonLayout = new QHBoxLayout;
	layout->addLayout(buttonLayout, 2, 0, 1, 2);
	buttonLayout->addStretch();
	QPushButton *cancelButton = new QPushButton(tr("Cancel"));
	connect(cancelButton, SIGNAL(clicked()), dialog, SLOT(reject()));
	buttonLayout->addWidget(cancelButton);
	QPushButton *okButton = new QPushButton(tr("Ok"));
	connect(okButton, SIGNAL(clicked()), dialog, SLOT(accept()));
	buttonLayout->addWidget(okButton);
	okButton->setDefault(true);
	int ret = dialog->exec();
	if (ret == QDialog::Rejected)
		return;
	for (int s = 0; s != surfaces_and_wells->count(); s++)
		if (surfaces_and_wells->item(s)->checkState() == Qt::Checked)
		{
			emit progress_append(">Start deleting geometric feature " + QString::number(s) + "...");
			QApplication::processEvents();
			Model.DeleteSurface(surfaces_and_wells->item(s)->text());
			emit progress_append(">...finished");
		}
}

void MainWindow::about()
{
	QMessageBox::about(this, tr("About MeshIT"),
										 tr("The tool <b>MeshIT</b> uses <a href='http://www.cs.cmu.edu/~quake/triangle.html' style='color: inherit; text-decoration: none;'><i>Triangle</i></a> "
												"and <a href='http://tetgen.berlios.de' style='color: inherit; text-decoration: none;'><i>Tetgen</i></a> "
												"to generate a quality tetrahedral mesh based on structural geological informations. "
												"This procedure is fully automized and needs at least scattered data point as input!<br><br>"
												"Main developers: <b>Mauro Cacace</b> and <b>Guido Bl&ouml;cher</b>.<br><br>"
												"Some extensions were added by <a href='https://www.perfacct.eu' style='color: inherit; text-decoration: none;'><b>PERFACCT</b></a>."));
}

void MainWindow::reset()
{
	Model.Intersections.clear();
	Model.Surfaces.clear();
	Model.Polylines.clear();
	Model.Mats.clear();
	if (Model.Mesh)
	{
		delete Model.Mesh;
		Model.Mesh = 0;
	}
	this->FillNameCombos();
}

void MainWindow::viewAxis()
{
	glWidget->drawAxis = this->viewAxisAct->isChecked();
	glWidget->updateGL();
}

void MainWindow::FillNameCombos()
{
	//	fill interpolation algotihm
	if (Model.intAlgorythm == "IDW")
		this->interpolationMethod->setCurrentIndex(0);
	if (Model.intAlgorythm == "SPLINE")
		this->interpolationMethod->setCurrentIndex(1);
	if (Model.intAlgorythm == "KRIGING")
		this->interpolationMethod->setCurrentIndex(2);
	//	clear old - units, faults, borders and wells
	this->unitsNamesCB->clear();
	this->faultsNamesCB->clear();
	this->bordersNamesCB->clear();
	this->wellsNamesCB->clear();
	this->matsNamesCB->clear();
	this->selectionComponent->clear();
	//	(re)fill - units, faults, borders
	for (int m = 0; m != Model.Mats.length(); m++)
		this->matsNamesCB->addItem(QIcon(":/images/tetrahedron.png"), "Material " + QString::number(m));
	for (int s = 0; s != Model.Surfaces.length(); s++)
	{
		if (Model.Surfaces[s].Type == "UNIT")
		{
			this->unitsNamesCB->addItem(Model.Surfaces[s].Name);
			this->selectionComponent->addItem(QIcon(":/images/units.png"), Model.Surfaces[s].Name);
		}
		if (Model.Surfaces[s].Type == "FAULT")
		{
			this->faultsNamesCB->addItem(Model.Surfaces[s].Name);
			this->selectionComponent->addItem(QIcon(":/images/faults.png"), Model.Surfaces[s].Name);
		}
		if (Model.Surfaces[s].Type == "BORDER")
		{
			this->bordersNamesCB->addItem(Model.Surfaces[s].Name);
			this->selectionComponent->addItem(QIcon(":/images/borders.png"), Model.Surfaces[s].Name);
		}
		if (Model.Surfaces[s].MaterialID >= 0)
			this->matsNamesCB->addItem(QIcon(":/images/triangle.png"), "Material " + QString::number(Model.Surfaces[s].MaterialID) + " (" + Model.Surfaces[s].Name + ")");
	}
	//	(re)fill - wells
	for (int p = 0; p != Model.Polylines.length(); p++)
		if (Model.Polylines[p].MaterialID >= 0)
			this->matsNamesCB->addItem(QIcon(":/images/line.png"), "Material " + QString::number(Model.Polylines[p].MaterialID) + " (" + Model.Polylines[p].Name + ")");
	for (int p = 0; p != Model.Polylines.length(); p++)
		if (Model.Polylines[p].Type == "WELL")
		{
			this->wellsNamesCB->addItem(Model.Polylines[p].Name);
			this->selectionComponent->addItem(QIcon(":/images/wells.png"), Model.Polylines[p].Name);
		}
	if (this->unitsNamesCB->count() > 1)
		this->unitsNamesCB->insertItem(0, "all");
	if (this->faultsNamesCB->count() > 1)
		this->faultsNamesCB->insertItem(0, "all");
	if (this->bordersNamesCB->count() > 1)
		this->bordersNamesCB->insertItem(0, "all");
	if (this->wellsNamesCB->count() > 1)
		this->wellsNamesCB->insertItem(0, "all");
	if (this->matsNamesCB->count() > 1)
		this->matsNamesCB->insertItem(0, "all");
}

void MainWindow::setUShowGBox(QString Name)
{
	if (Name == "all")
	{
		int count = 0, dSD = 0, dCH = 0, dF = 0, dE = 0, dIE = 0, dIV = 0;
		for (int s = 0; s != Model.Surfaces.length(); s++)
			if (Model.Surfaces[s].Type == "UNIT")
			{
				count++;
				if (Model.Surfaces[s].drawScatteredData)
					dSD++;
				if (Model.Surfaces[s].drawConvexHull)
					dCH++;
				if (Model.Surfaces[s].drawFaces)
					dF++;
				if (Model.Surfaces[s].drawEdges)
					dE++;
				if (Model.Surfaces[s].drawIntEdges)
					dIE++;
				if (Model.Surfaces[s].drawIntVertices)
					dIV++;
			}
		//	scattered data points
		if (dSD == 0)
			this->unitsScatteredDataPoints->setCheckState(0, Qt::Unchecked);
		else if (dSD == count)
			this->unitsScatteredDataPoints->setCheckState(0, Qt::Checked);
		else
			this->unitsScatteredDataPoints->setCheckState(0, Qt::PartiallyChecked);
		//	convex hull
		if (dCH == 0)
			this->unitsConvexHull->setCheckState(0, Qt::Unchecked);
		else if (dCH == count)
			this->unitsConvexHull->setCheckState(0, Qt::Checked);
		else
			this->unitsConvexHull->setCheckState(0, Qt::PartiallyChecked);
		//	faces
		if (dF == 0)
			this->unitsFaces->setCheckState(0, Qt::Unchecked);
		else if (dF == count)
			this->unitsFaces->setCheckState(0, Qt::Checked);
		else
			this->unitsFaces->setCheckState(0, Qt::PartiallyChecked);
		//	edges
		if (dE == 0)
			this->unitsEdges->setCheckState(0, Qt::Unchecked);
		else if (dE == count)
			this->unitsEdges->setCheckState(0, Qt::Checked);
		else
			this->unitsEdges->setCheckState(0, Qt::PartiallyChecked);
		//	intersection edges
		if (dIE == 0)
			this->unitsIntersectionEdges->setCheckState(0, Qt::Unchecked);
		else if (dIE == count)
			this->unitsIntersectionEdges->setCheckState(0, Qt::Checked);
		else
			this->unitsIntersectionEdges->setCheckState(0, Qt::PartiallyChecked);
		//	intersection vertices
		if (dIV == 0)
			this->unitsIntersectionVertices->setCheckState(0, Qt::Unchecked);
		else if (dIV == count)
			this->unitsIntersectionVertices->setCheckState(0, Qt::Checked);
		else
			this->unitsIntersectionVertices->setCheckState(0, Qt::PartiallyChecked);
	}
	else
	{
		int i = -1;
		for (int s = 0; s != Model.Surfaces.length(); s++)
			if (Model.Surfaces[s].Name == Name && Model.Surfaces[s].Type == "UNIT")
				i = s;

		if (i > -1)
		{
			// initialize all
			this->unitsScatteredDataPoints->setCheckState(0, Qt::Unchecked);
			this->unitsConvexHull->setCheckState(0, Qt::Unchecked);
			this->unitsFaces->setCheckState(0, Qt::Unchecked);
			this->unitsEdges->setCheckState(0, Qt::Unchecked);
			this->unitsIntersectionEdges->setCheckState(0, Qt::Unchecked);
			this->unitsIntersectionVertices->setCheckState(0, Qt::Unchecked);
			//	scattered data points
			if (Model.Surfaces[i].drawScatteredData)
				this->unitsScatteredDataPoints->setCheckState(0, Qt::Checked);
			//	convex hull
			if (Model.Surfaces[i].drawConvexHull)
				this->unitsConvexHull->setCheckState(0, Qt::Checked);
			// faces
			if (Model.Surfaces[i].drawFaces)
				this->unitsFaces->setCheckState(0, Qt::Checked);
			// edges
			if (Model.Surfaces[i].drawEdges)
				this->unitsEdges->setCheckState(0, Qt::Checked);
			// intersection edges
			if (Model.Surfaces[i].drawIntEdges)
				this->unitsIntersectionEdges->setCheckState(0, Qt::Checked);
			// intersection vertices
			if (Model.Surfaces[i].drawIntVertices)
				this->unitsIntersectionVertices->setCheckState(0, Qt::Checked);
		}
	}
}

void MainWindow::setUMeshes()
{
	if (this->unitsNamesCB->currentText() == "all")
	{
		for (int s = 0; s != Model.Surfaces.length(); s++)
			if (Model.Surfaces[s].Type == "UNIT")
			{
				//	scattered data points
				if (this->unitsScatteredDataPoints->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawScatteredData = true;
				if (this->unitsScatteredDataPoints->checkState(0) == Qt::Unchecked)
					Model.Surfaces[s].drawScatteredData = false;
				// convex hull
				if (this->unitsConvexHull->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawConvexHull = true;
				if (this->unitsConvexHull->checkState(0) == Qt::Unchecked)
					Model.Surfaces[s].drawConvexHull = false;
				//	faces
				if (this->unitsFaces->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawFaces = true;
				if (this->unitsFaces->checkState(0) == Qt::Unchecked)
					Model.Surfaces[s].drawFaces = false;
				// edges
				if (this->unitsEdges->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawEdges = true;
				if (this->unitsEdges->checkState(0) == Qt::Unchecked)
					Model.Surfaces[s].drawEdges = false;
				// intersection edges
				if (this->unitsIntersectionEdges->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawIntEdges = true;
				if (this->unitsIntersectionEdges->checkState(0) == Qt::Unchecked)
					Model.Surfaces[s].drawIntEdges = false;
				// intersection vertices
				if (this->unitsIntersectionVertices->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawIntVertices = true;
				if (this->unitsIntersectionVertices->checkState(0) == Qt::Unchecked)
					Model.Surfaces[s].drawIntVertices = false;
			}
	}
	else
	{
		int i = -1;
		for (int s = 0; s != Model.Surfaces.length(); s++)
			if (Model.Surfaces[s].Name == unitsNamesCB->currentText())
				i = s;
		if (i > -1)
		{
			// initialize all
			Model.Surfaces[i].drawScatteredData = false;
			Model.Surfaces[i].drawConvexHull = false;
			Model.Surfaces[i].drawFaces = false;
			Model.Surfaces[i].drawEdges = false;
			Model.Surfaces[i].drawIntEdges = false;
			Model.Surfaces[i].drawIntVertices = false;
			// scattered data points
			if (this->unitsScatteredDataPoints->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawScatteredData = true;
			// convex hull
			if (this->unitsConvexHull->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawConvexHull = true;
			// faces
			if (this->unitsFaces->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawFaces = true;
			// edges
			if (this->unitsEdges->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawEdges = true;
			// intersection edges
			if (this->unitsIntersectionEdges->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawIntEdges = true;
			// intersection vertices
			if (this->unitsIntersectionVertices->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawIntVertices = true;
		}
	}
	glWidget->update();
}

void MainWindow::setFShowGBox(QString Name)
{
	if (Name == "all")
	{
		int count = 0, dSD = 0, dCH = 0, dF = 0, dE = 0, dIE = 0, dIV = 0;
		for (int s = 0; s != Model.Surfaces.length(); s++)
			if (Model.Surfaces[s].Type == "FAULT")
			{
				count++;
				if (Model.Surfaces[s].drawScatteredData)
					dSD++;
				if (Model.Surfaces[s].drawConvexHull)
					dCH++;
				if (Model.Surfaces[s].drawFaces)
					dF++;
				if (Model.Surfaces[s].drawEdges)
					dE++;
				if (Model.Surfaces[s].drawIntEdges)
					dIE++;
				if (Model.Surfaces[s].drawIntVertices)
					dIV++;
			}
		// scattered data input points
		if (dSD == 0)
			this->faultsScatteredDataPoints->setCheckState(0, Qt::Unchecked);
		else if (dSD == count)
			this->faultsScatteredDataPoints->setCheckState(0, Qt::Checked);
		else
			this->faultsScatteredDataPoints->setCheckState(0, Qt::PartiallyChecked);
		// convex hull
		if (dCH == 0)
			this->faultsConvexHull->setCheckState(0, Qt::Unchecked);
		else if (dCH == count)
			this->faultsConvexHull->setCheckState(0, Qt::Checked);
		else
			this->faultsConvexHull->setCheckState(0, Qt::PartiallyChecked);
		// faces
		if (dF == 0)
			this->faultsFaces->setCheckState(0, Qt::Unchecked);
		else if (dF == count)
			this->faultsFaces->setCheckState(0, Qt::Checked);
		else
			this->faultsFaces->setCheckState(0, Qt::PartiallyChecked);
		// edges
		if (dE == 0)
			this->faultsEdges->setCheckState(0, Qt::Unchecked);
		else if (dE == count)
			this->faultsEdges->setCheckState(0, Qt::Checked);
		else
			this->faultsEdges->setCheckState(0, Qt::PartiallyChecked);
		// intersection edges
		if (dIE == 0)
			this->faultsIntersectionEdges->setCheckState(0, Qt::Unchecked);
		else if (dIE == count)
			this->faultsIntersectionEdges->setCheckState(0, Qt::Checked);
		else
			this->faultsIntersectionEdges->setCheckState(0, Qt::PartiallyChecked);
		// intersection vertices
		if (dIV == 0)
			this->faultsIntersectionVertices->setCheckState(0, Qt::Unchecked);
		else if (dIV == count)
			this->faultsIntersectionVertices->setCheckState(0, Qt::Checked);
		else
			this->faultsIntersectionVertices->setCheckState(0, Qt::PartiallyChecked);
	}
	else
	{
		int i = -1;
		for (int s = 0; s != Model.Surfaces.length(); s++)
			if (Model.Surfaces[s].Name == Name && Model.Surfaces[s].Type == "FAULT")
				i = s;
		if (i > -1)
		{
			// initialize all
			this->faultsScatteredDataPoints->setCheckState(0, Qt::Unchecked);
			this->faultsConvexHull->setCheckState(0, Qt::Unchecked);
			this->faultsFaces->setCheckState(0, Qt::Unchecked);
			this->faultsEdges->setCheckState(0, Qt::Unchecked);
			this->faultsIntersectionEdges->setCheckState(0, Qt::Unchecked);
			this->faultsIntersectionVertices->setCheckState(0, Qt::Unchecked);
			// scattered data points
			if (Model.Surfaces[i].drawScatteredData)
				this->faultsScatteredDataPoints->setCheckState(0, Qt::Checked);
			// convex hull
			if (Model.Surfaces[i].drawConvexHull)
				this->faultsConvexHull->setCheckState(0, Qt::Checked);
			// faces
			if (Model.Surfaces[i].drawFaces)
				this->faultsFaces->setCheckState(0, Qt::Checked);
			// edges
			if (Model.Surfaces[i].drawEdges)
				this->faultsEdges->setCheckState(0, Qt::Checked);
			// intersection edges
			if (Model.Surfaces[i].drawIntEdges)
				this->faultsIntersectionEdges->setCheckState(0, Qt::Checked);
			// intersection vertices
			if (Model.Surfaces[i].drawIntVertices)
				this->faultsIntersectionVertices->setCheckState(0, Qt::Checked);
		}
	}
}

void MainWindow::setFMeshes()
{
	if (this->faultsNamesCB->currentText() == "all")
	{
		for (int s = 0; s != Model.Surfaces.length(); s++)
			if (Model.Surfaces[s].Type == "FAULT")
			{
				// initialize all
				Model.Surfaces[s].drawScatteredData = false;
				Model.Surfaces[s].drawConvexHull = false;
				Model.Surfaces[s].drawFaces = false;
				Model.Surfaces[s].drawEdges = false;
				Model.Surfaces[s].drawIntEdges = false;
				Model.Surfaces[s].drawIntVertices = false;
				// scattered data points
				if (this->faultsScatteredDataPoints->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawScatteredData = true;
				// convex hull
				if (this->faultsConvexHull->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawConvexHull = true;
				// faces
				if (this->faultsFaces->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawFaces = true;
				// edges
				if (this->faultsEdges->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawEdges = true;
				// intersection edges
				if (this->faultsIntersectionEdges->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawIntEdges = true;
				// intersection vertices
				if (this->faultsIntersectionVertices->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawIntVertices = true;
			}
	}
	else
	{
		int i = -1;
		for (int s = 0; s != Model.Surfaces.length(); s++)
			if (Model.Surfaces[s].Name == faultsNamesCB->currentText())
				i = s;
		if (i > -1)
		{
			// initialize all
			Model.Surfaces[i].drawScatteredData = false;
			Model.Surfaces[i].drawConvexHull = false;
			Model.Surfaces[i].drawFaces = false;
			Model.Surfaces[i].drawEdges = false;
			Model.Surfaces[i].drawIntEdges = false;
			Model.Surfaces[i].drawIntVertices = false;
			// scattered data points
			if (this->faultsScatteredDataPoints->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawScatteredData = true;
			// convex hull
			if (this->faultsConvexHull->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawConvexHull = true;
			// faces
			if (this->faultsFaces->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawFaces = true;
			// edges
			if (this->faultsEdges->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawEdges = true;
			// intersection edges
			if (this->faultsIntersectionEdges->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawIntEdges = true;
			// intersection vertices
			if (this->faultsIntersectionVertices->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawIntVertices = true;
		}
	}
	glWidget->update();
}

void MainWindow::setBShowGBox(QString Name)
{
	if (Name == "all")
	{
		int count = 0, dSD = 0, dCH = 0, dF = 0, dE = 0, dIE = 0, dIV = 0;
		for (int s = 0; s != Model.Surfaces.length(); s++)
			if (Model.Surfaces[s].Type == "BORDER")
			{
				count++;
				if (Model.Surfaces[s].drawScatteredData)
					dSD++;
				if (Model.Surfaces[s].drawConvexHull)
					dCH++;
				if (Model.Surfaces[s].drawFaces)
					dF++;
				if (Model.Surfaces[s].drawEdges)
					dE++;
				if (Model.Surfaces[s].drawIntEdges)
					dIE++;
				if (Model.Surfaces[s].drawIntVertices)
					dIV++;
			}
		// scattered data points
		if (dSD == 0)
			this->bordersScatteredDataPoints->setCheckState(0, Qt::Unchecked);
		else if (dSD == count)
			this->bordersScatteredDataPoints->setCheckState(0, Qt::Checked);
		else
			this->bordersScatteredDataPoints->setCheckState(0, Qt::PartiallyChecked);
		// convex hull
		if (dCH == 0)
			this->bordersConvexHull->setCheckState(0, Qt::Unchecked);
		else if (dCH == count)
			this->bordersConvexHull->setCheckState(0, Qt::Checked);
		else
			this->bordersConvexHull->setCheckState(0, Qt::PartiallyChecked);
		// faces
		if (dF == 0)
			this->bordersFaces->setCheckState(0, Qt::Unchecked);
		else if (dF == count)
			this->bordersFaces->setCheckState(0, Qt::Checked);
		else
			this->bordersFaces->setCheckState(0, Qt::PartiallyChecked);
		// edges
		if (dE == 0)
			this->bordersEdges->setCheckState(0, Qt::Unchecked);
		else if (dE == count)
			this->bordersEdges->setCheckState(0, Qt::Checked);
		else
			this->bordersEdges->setCheckState(0, Qt::PartiallyChecked);
		// intersection edges
		if (dIE == 0)
			this->bordersIntersectionEdges->setCheckState(0, Qt::Unchecked);
		else if (dIE == count)
			this->bordersIntersectionEdges->setCheckState(0, Qt::Checked);
		else
			this->bordersIntersectionEdges->setCheckState(0, Qt::PartiallyChecked);
		// intersection vertices
		if (dIV == 0)
			this->bordersIntersectionVertices->setCheckState(0, Qt::Unchecked);
		else if (dIV == count)
			this->bordersIntersectionVertices->setCheckState(0, Qt::Checked);
		else
			this->bordersIntersectionVertices->setCheckState(0, Qt::PartiallyChecked);
	}
	else
	{
		int i = -1;
		for (int s = 0; s != Model.Surfaces.length(); s++)
			if (Model.Surfaces[s].Name == Name && Model.Surfaces[s].Type == "BORDER")
				i = s;
		if (i > -1)
		{
			// initialize all
			this->bordersScatteredDataPoints->setCheckState(0, Qt::Unchecked);
			this->bordersConvexHull->setCheckState(0, Qt::Unchecked);
			this->bordersFaces->setCheckState(0, Qt::Unchecked);
			this->bordersEdges->setCheckState(0, Qt::Unchecked);
			this->bordersIntersectionEdges->setCheckState(0, Qt::Unchecked);
			this->bordersIntersectionVertices->setCheckState(0, Qt::Unchecked);
			// scattered data points
			if (Model.Surfaces[i].drawScatteredData)
				this->bordersScatteredDataPoints->setCheckState(0, Qt::Checked);
			// convex hull
			if (Model.Surfaces[i].drawConvexHull)
				this->bordersConvexHull->setCheckState(0, Qt::Checked);
			// faces
			if (Model.Surfaces[i].drawFaces)
				this->bordersFaces->setCheckState(0, Qt::Checked);
			// edges
			if (Model.Surfaces[i].drawEdges)
				this->bordersEdges->setCheckState(0, Qt::Checked);
			// intersection edges
			if (Model.Surfaces[i].drawIntEdges)
				this->bordersIntersectionEdges->setCheckState(0, Qt::Checked);
			// intersection vertices
			if (Model.Surfaces[i].drawIntVertices)
				this->bordersIntersectionVertices->setCheckState(0, Qt::Checked);
		}
	}
}

void MainWindow::setBMeshes()
{
	if (this->bordersNamesCB->currentText() == "all")
	{
		for (int s = 0; s != Model.Surfaces.length(); s++)
			if (Model.Surfaces[s].Type == "BORDER")
			{
				// initialize all
				Model.Surfaces[s].drawScatteredData = false;
				Model.Surfaces[s].drawConvexHull = false;
				Model.Surfaces[s].drawFaces = false;
				Model.Surfaces[s].drawEdges = false;
				Model.Surfaces[s].drawIntEdges = false;
				Model.Surfaces[s].drawIntVertices = false;
				// scattered data points
				if (this->bordersScatteredDataPoints->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawScatteredData = true;
				// convex hull
				if (this->bordersConvexHull->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawConvexHull = true;
				// faces
				if (this->bordersFaces->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawFaces = true;
				// edges
				if (this->bordersEdges->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawEdges = true;
				// intersection edges
				if (this->bordersIntersectionEdges->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawIntEdges = true;
				// intersection vertices
				if (this->bordersIntersectionVertices->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawIntVertices = true;
			}
	}
	else
	{
		int i = -1;
		for (int s = 0; s != Model.Surfaces.length(); s++)
			if (Model.Surfaces[s].Name == bordersNamesCB->currentText())
				i = s;
		if (i > -1)
		{
			// initialize all
			Model.Surfaces[i].drawScatteredData = false;
			Model.Surfaces[i].drawConvexHull = false;
			Model.Surfaces[i].drawFaces = false;
			Model.Surfaces[i].drawEdges = false;
			Model.Surfaces[i].drawIntEdges = false;
			Model.Surfaces[i].drawIntVertices = false;
			// scattered data points
			if (this->bordersScatteredDataPoints->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawScatteredData = true;
			// convex hull
			if (this->bordersConvexHull->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawConvexHull = true;
			// faces
			if (this->bordersFaces->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawFaces = true;
			// edges
			if (this->bordersEdges->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawEdges = true;
			// intersection edges
			if (this->bordersIntersectionEdges->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawIntEdges = true;
			// intersection vertices
			if (this->bordersIntersectionVertices->checkState(0) == Qt::Checked)
				Model.Surfaces[i].drawIntVertices = true;
		}
	}
	glWidget->update();
}

void MainWindow::setWShowGBox(QString Name)
{
	if (Name == "all")
	{
		int count = 0, dSD = 0, dE = 0, dV = 0, dIV = 0;
		for (int p = 0; p != Model.Polylines.length(); p++)
		{
			count++;
			if (Model.Polylines[p].drawScatteredData)
				dSD++;
			if (Model.Polylines[p].drawEdges)
				dE++;
			if (Model.Polylines[p].drawVertices)
				dV++;
			if (Model.Polylines[p].drawIntVertices)
				dIV++;
		}
		// scattered data points
		if (dSD == 0)
			this->wellsScatteredDataPoints->setCheckState(0, Qt::Unchecked);
		else if (dSD == count)
			this->wellsScatteredDataPoints->setCheckState(0, Qt::Checked);
		else
			this->wellsScatteredDataPoints->setCheckState(0, Qt::PartiallyChecked);
		// edges
		if (dE == 0)
			this->wellsEdges->setCheckState(0, Qt::Unchecked);
		else if (dE == count)
			this->wellsEdges->setCheckState(0, Qt::Checked);
		else
			this->wellsEdges->setCheckState(0, Qt::PartiallyChecked);
		// vertices
		if (dV == 0)
			this->wellsVertices->setCheckState(0, Qt::Unchecked);
		else if (dV == count)
			this->wellsVertices->setCheckState(0, Qt::Checked);
		else
			this->wellsVertices->setCheckState(0, Qt::PartiallyChecked);
		// intersection vertices
		if (dIV == 0)
			this->wellsIntersectionVertices->setCheckState(0, Qt::Unchecked);
		else if (dIV == count)
			this->wellsIntersectionVertices->setCheckState(0, Qt::Checked);
		else
			this->wellsIntersectionVertices->setCheckState(0, Qt::PartiallyChecked);
	}
	else
	{
		int i = -1;
		for (int p = 0; p != Model.Polylines.length(); p++)
			if (Model.Polylines[p].Name == Name && Model.Polylines[p].Type == "WELL")
				i = p;
		if (i > -1)
		{
			// initialize all
			this->wellsScatteredDataPoints->setCheckState(0, Qt::Unchecked);
			this->wellsEdges->setCheckState(0, Qt::Unchecked);
			this->wellsVertices->setCheckState(0, Qt::Unchecked);
			this->wellsIntersectionVertices->setCheckState(0, Qt::Unchecked);
			// scattered data points
			if (Model.Polylines[i].drawScatteredData)
				this->wellsScatteredDataPoints->setCheckState(0, Qt::Checked);
			// edges
			if (Model.Polylines[i].drawEdges)
				this->wellsEdges->setCheckState(0, Qt::Checked);
			// vertices
			if (Model.Polylines[i].drawVertices)
				this->wellsVertices->setCheckState(0, Qt::Checked);
			// intersection vertices
			if (Model.Polylines[i].drawIntVertices)
				this->wellsIntersectionVertices->setCheckState(0, Qt::Checked);
		}
	}
}

void MainWindow::setWMeshes()
{
	if (this->wellsNamesCB->currentText() == "all")
	{
		for (int p = 0; p != Model.Polylines.length(); p++)
			if (Model.Polylines[p].Type == "WELL")
			{
				// initialize all
				Model.Polylines[p].drawScatteredData = false;
				Model.Polylines[p].drawEdges = false;
				Model.Polylines[p].drawVertices = false;
				Model.Polylines[p].drawIntVertices = false;
				// scattered data points
				if (this->wellsScatteredDataPoints->checkState(0) == Qt::Checked)
					Model.Polylines[p].drawScatteredData = true;
				// edges
				if (this->wellsEdges->checkState(0) == Qt::Checked)
					Model.Polylines[p].drawEdges = true;
				// vertices
				if (this->wellsVertices->checkState(0) == Qt::Checked)
					Model.Polylines[p].drawVertices = true;
				// intersection vertices
				if (this->wellsIntersectionVertices->checkState(0) == Qt::Checked)
					Model.Polylines[p].drawIntVertices = true;
			}
	}
	else
	{
		int i = -1;
		for (int p = 0; p != Model.Polylines.length(); p++)
			if (Model.Polylines[p].Name == wellsNamesCB->currentText())
				i = p;
		if (i > -1)
		{
			// initialize all
			Model.Polylines[i].drawScatteredData = false;
			Model.Polylines[i].drawEdges = false;
			Model.Polylines[i].drawVertices = false;
			Model.Polylines[i].drawIntVertices = false;
			// scattered data points
			if (this->wellsScatteredDataPoints->checkState(0) == Qt::Checked)
				Model.Polylines[i].drawScatteredData = true;
			// edges
			if (this->wellsEdges->checkState(0) == Qt::Checked)
				Model.Polylines[i].drawEdges = true;
			// vertices
			if (this->wellsVertices->checkState(0) == Qt::Checked)
				Model.Polylines[i].drawVertices = true;
			// intersection vertices
			if (this->wellsIntersectionVertices->checkState(0) == Qt::Checked)
				Model.Polylines[i].drawIntVertices = true;
		}
	}
	glWidget->update();
}

void MainWindow::setMShowGBox(QString Name)
{
	if (Name == "all")
	{
		int count = 0, dF = 0, dE = 0;
		for (int m = 0; m != Model.Mats.length(); m++)
		{
			count++;
			if (Model.Mats[m].drawMatFaces)
				dF++;
			if (Model.Mats[m].drawMatEdges)
				dE++;
		}
		for (int s = 0; s != Model.Surfaces.length(); s++)
			if (Model.Surfaces[s].MaterialID >= 0)
			{
				count++;
				if (Model.Surfaces[s].drawMatFaces)
					dF++;
				if (Model.Surfaces[s].drawMatEdges)
					dE++;
			}
		if (dF == 0)
			this->matsFaces->setCheckState(0, Qt::Unchecked);
		else if (dF == count)
			this->matsFaces->setCheckState(0, Qt::Checked);
		else
			this->matsFaces->setCheckState(0, Qt::PartiallyChecked);
		if (dE == 0)
			this->matsEdges->setCheckState(0, Qt::Unchecked);
		else if (dE == count)
			this->matsEdges->setCheckState(0, Qt::Checked);
		else
			this->matsEdges->setCheckState(0, Qt::PartiallyChecked);
	}
	else
	{
		int i = Name.section(" ", 1, 1).toInt();
		if (i < Model.Mats.length())
		{
			// initialize all materials
			this->matsFaces->setCheckState(0, Qt::Unchecked);
			this->matsEdges->setCheckState(0, Qt::Unchecked);
			// material faces
			if (Model.Mats[i].drawMatFaces)
				this->matsFaces->setCheckState(0, Qt::Checked);
			// material edges
			if (Model.Mats[i].drawMatEdges)
				this->matsEdges->setCheckState(0, Qt::Checked);
		}
		else
		{
			i = -1;
			for (int s = 0; s != Model.Surfaces.length(); s++)
				if (Model.Surfaces[s].Name == Name.section("(", 1, 1).section(")", 0, 0) && Model.Surfaces[s].MaterialID >= 0)
					i = s;
			if (i > -1)
			{
				// initialize all materials
				this->matsFaces->setCheckState(0, Qt::Unchecked);
				this->matsEdges->setCheckState(0, Qt::Unchecked);
				// material faces
				if (Model.Surfaces[i].drawMatFaces)
					this->matsFaces->setCheckState(0, Qt::Checked);
				// material edges
				if (Model.Surfaces[i].drawMatEdges)
					this->matsEdges->setCheckState(0, Qt::Checked);
			}
		}
	}
}

void MainWindow::setMMeshes()
{
	if (this->matsNamesCB->currentText() == "all")
	{
		for (int m = 0; m != Model.Mats.length(); m++)
		{
			// initial all materials
			Model.Mats[m].drawMatFaces = false;
			Model.Mats[m].drawMatEdges = false;
			// material faces
			if (this->matsFaces->checkState(0) == Qt::Checked)
				Model.Mats[m].drawMatFaces = true;
			// material edges
			if (this->matsEdges->checkState(0) == Qt::Checked)
				Model.Mats[m].drawMatEdges = true;
		}
		for (int s = 0; s != Model.Surfaces.length(); s++)
			if (Model.Surfaces[s].MaterialID >= 0)
			{
				// initialize all materials
				Model.Surfaces[s].drawMatFaces = false;
				Model.Surfaces[s].drawMatEdges = false;
				// material faces
				if (this->matsFaces->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawMatFaces = true;
				// material edges
				if (this->matsEdges->checkState(0) == Qt::Checked)
					Model.Surfaces[s].drawMatEdges = true;
			}
	}
	else
	{
		int i = matsNamesCB->currentText().section(" ", 1, 1).toInt();
		if (i < Model.Mats.length())
		{
			// initial all materials
			Model.Mats[i].drawMatFaces = false;
			Model.Mats[i].drawMatEdges = false;
			// material faces
			if (this->matsFaces->checkState(0) == Qt::Checked)
				Model.Mats[i].drawMatFaces = true;
			// material edges
			if (this->matsEdges->checkState(0) == Qt::Checked)
				Model.Mats[i].drawMatEdges = true;
		}
		else
		{
			i = -1;
			for (int s = 0; s != Model.Surfaces.length(); s++)
				if (Model.Surfaces[s].Name == matsNamesCB->currentText().section("(", 1, 1).section(")", 0, 0))
					i = s;
			if (i > -1)
			{
				// initialize all materials
				Model.Surfaces[i].drawMatFaces = false;
				Model.Surfaces[i].drawMatEdges = false;
				// material faces
				if (this->matsFaces->checkState(0) == Qt::Checked)
					Model.Surfaces[i].drawMatFaces = true;
				// material edges
				if (this->matsEdges->checkState(0) == Qt::Checked)
					Model.Surfaces[i].drawMatEdges = true;
			}
		}
	}
	this->callMakeTets();
}

void MainWindow::setSShowGBox(QString MeshName)
{
	// prevent showing of selection if the box has not focus
	QComboBox *senderCB = qobject_cast<QComboBox *>(sender());
	if (senderCB && !senderCB->isVisible())
		MeshName = "DoNotShowAnything";

	for (int s = 0; s != Model.Surfaces.length(); s++)
		if (Model.Surfaces[s].Name == MeshName)
			Model.Surfaces[s].drawConstraints = true;
		else
			Model.Surfaces[s].drawConstraints = false;
	for (int p = 0; p != Model.Polylines.length(); p++)
		if (Model.Polylines[p].Name == MeshName)
			Model.Polylines[p].drawConstraints = true;
		else
			Model.Polylines[p].drawConstraints = false;
	glWidget->update();
}

void MainWindow::interpolationSetMethod(QString method)
{
	Model.intAlgorythm = method;
}

void MainWindow::material3dFillValue(int Location)
{

	if (Location != -1)
	{
		this->material3dXSlider->blockSignals(true);
		this->material3dYSlider->blockSignals(true);
		this->material3dZSlider->blockSignals(true);
		this->material3dXValue->blockSignals(true);
		this->material3dYValue->blockSignals(true);
		this->material3dZValue->blockSignals(true);
		this->material3dXSlider->setRange(256 * (Model.min.x() - Model.shift.x()) * Model.scale, 256 * (Model.max.x() - Model.shift.x()) * Model.scale);
		this->material3dYSlider->setRange(256 * (Model.min.y() - Model.shift.y()) * Model.scale, 256 * (Model.max.y() - Model.shift.y()) * Model.scale);
		this->material3dZSlider->setRange(256 * (Model.min.z() - Model.shift.z()) * Model.scale, 256 * (Model.max.z() - Model.shift.z()) * Model.scale);
		this->material3dXSlider->setSingleStep(1);
		this->material3dYSlider->setSingleStep(1);
		this->material3dZSlider->setSingleStep(1);
		this->material3dXValue->setRange(Model.min.x(), Model.max.x());
		this->material3dYValue->setRange(Model.min.y(), Model.max.y());
		this->material3dZValue->setRange(Model.min.z(), Model.max.z());
		this->material3dXValue->setSingleStep((Model.max.x() - Model.min.x()) / 512.0);
		this->material3dYValue->setSingleStep((Model.max.y() - Model.min.y()) / 512.0);
		this->material3dZValue->setSingleStep((Model.max.z() - Model.min.z()) / 512.0);
		this->material3dXSlider->blockSignals(false);
		this->material3dYSlider->blockSignals(false);
		this->material3dZSlider->blockSignals(false);
		this->material3dXValue->blockSignals(false);
		this->material3dYValue->blockSignals(false);
		this->material3dZValue->blockSignals(false);
		this->material3dXValue->setValue(Model.Mats[this->material3dList->currentRow()].Locations[Location].x() / Model.scale + Model.shift.x());
		this->material3dYValue->setValue(Model.Mats[this->material3dList->currentRow()].Locations[Location].y() / Model.scale + Model.shift.y());
		this->material3dZValue->setValue(Model.Mats[this->material3dList->currentRow()].Locations[Location].z() / Model.scale + Model.shift.z());
	}
	this->callMakeMats();
}

void MainWindow::material3dSetLocationFromDSpinBox(double value)
{
	if (this->material3dListLocation->currentRow() == -1)
	{
		QMessageBox msgBox;
		msgBox.setText("One material and one loacation must be selected to change the position!");
		msgBox.exec();
		return;
	}
	QDoubleSpinBox *senderDSB = qobject_cast<QDoubleSpinBox *>(sender());
	if (senderDSB)
	{
		if (senderDSB->property("location") == "X")
		{
			double scaledValue = (value - Model.shift.x()) * Model.scale;
			Model.Mats[this->material3dList->currentRow()].Locations[this->material3dListLocation->currentRow()].setX(scaledValue);
			this->material3dXSlider->blockSignals(true);
			this->material3dXSlider->setValue(256 * scaledValue);
			this->material3dXSlider->blockSignals(false);
		}
		else if (senderDSB->property("location") == "Y")
		{
			double scaledValue = (value - Model.shift.y()) * Model.scale;
			Model.Mats[this->material3dList->currentRow()].Locations[this->material3dListLocation->currentRow()].setY(scaledValue);
			this->material3dYSlider->blockSignals(true);
			this->material3dYSlider->setValue(256 * scaledValue);
			this->material3dYSlider->blockSignals(false);
		}
		else if (senderDSB->property("location") == "Z")
		{
			double scaledValue = (value - Model.shift.z()) * Model.scale;
			Model.Mats[this->material3dList->currentRow()].Locations[this->material3dListLocation->currentRow()].setZ(scaledValue);
			this->material3dZSlider->blockSignals(true);
			this->material3dZSlider->setValue(256 * scaledValue);
			this->material3dZSlider->blockSignals(false);
		}
	}
	this->callMakeMats();
}

void MainWindow::material3dSetLocationFromSlider(int value)
{
	if (this->material3dListLocation->currentRow() == -1)
	{
		QMessageBox msgBox;
		msgBox.setText("One material and one loacation must be selected to change the position!");
		msgBox.exec();
		return;
	}
	QSlider *senderS = qobject_cast<QSlider *>(sender());
	if (senderS)
	{
		if (senderS->property("location") == "X")
		{
			double scaledValue = value / 256.0;
			Model.Mats[this->material3dList->currentRow()].Locations[this->material3dListLocation->currentRow()].setX(scaledValue);
			this->material3dXValue->blockSignals(true);
			this->material3dXValue->setValue(scaledValue / Model.scale + Model.shift.x());
			this->material3dXValue->blockSignals(false);
		}
		else if (senderS->property("location") == "Y")
		{
			double scaledValue = value / 256.0;
			Model.Mats[this->material3dList->currentRow()].Locations[this->material3dListLocation->currentRow()].setY(scaledValue);
			this->material3dYValue->blockSignals(true);
			this->material3dYValue->setValue(scaledValue / Model.scale + Model.shift.y());
			this->material3dYValue->blockSignals(false);
		}
		else if (senderS->property("location") == "Z")
		{
			double scaledValue = value / 256.0;
			Model.Mats[this->material3dList->currentRow()].Locations[this->material3dListLocation->currentRow()].setZ(scaledValue);
			this->material3dZValue->blockSignals(true);
			this->material3dZValue->setValue(scaledValue / Model.scale + Model.shift.z());
			this->material3dZValue->blockSignals(false);
		}
	}
	this->callMakeMats();
}

void MainWindow::interpolationFill()
{
}

void MainWindow::ExportRotationAngelUpdate(double angel)
{
	Model.ExportRotationAngle = angel;
}

void MainWindow::refinementFill()
{
	this->refinementTableWidget->clearContents();
	this->refinementTableWidget->setColumnCount(2);
	this->refinementTableWidget->setRowCount(Model.Surfaces.length() + Model.Polylines.length());
	// surfaces
	for (int s = 0; s != Model.Surfaces.length(); s++)
	{
		if (Model.Surfaces[s].Type == "UNIT")
			this->refinementTableWidgetItemObject = new QTableWidgetItem(QIcon(":/images/units.png"), Model.Surfaces[s].Name, 0);
		if (Model.Surfaces[s].Type == "FAULT")
			this->refinementTableWidgetItemObject = new QTableWidgetItem(QIcon(":/images/faults.png"), Model.Surfaces[s].Name, 0);
		if (Model.Surfaces[s].Type == "BORDER")
			this->refinementTableWidgetItemObject = new QTableWidgetItem(QIcon(":/images/borders.png"), Model.Surfaces[s].Name, 0);
		this->refinementTableWidget->setItem(s, 0, this->refinementTableWidgetItemObject);
		this->refinementTableWidgetItemValue = new QDoubleSpinBox();
		this->refinementTableWidgetItemValue->setMaximum(1e7);
		this->refinementTableWidgetItemValue->setMinimum(0);
		this->refinementTableWidgetItemValue->setProperty("type", "Surface");
		this->refinementTableWidgetItemValue->setProperty("number", s);

		connect(refinementTableWidgetItemValue, SIGNAL(valueChanged(double)), this, SLOT(refinementUpdate(double)));

		this->refinementTableWidgetItemValue->setDecimals(3);
		this->refinementTableWidgetItemValue->setValue(Model.Surfaces[s].size / Model.scale);
		this->refinementTableWidget->setCellWidget(s, 1, this->refinementTableWidgetItemValue);
	}
	// polylines
	for (int p = 0; p != Model.Polylines.length(); p++)
	{
		this->refinementTableWidgetItemObject = new QTableWidgetItem(QIcon(":/images/wells.png"), Model.Polylines[p].Name, 2);
		this->refinementTableWidgetItemObject->setFlags(this->refinementTableWidgetItemObject->flags() & ~Qt::ItemIsEditable);
		this->refinementTableWidget->setItem(p + Model.Surfaces.length(), 0, this->refinementTableWidgetItemObject);
		this->refinementTableWidgetItemValue = new QDoubleSpinBox();
		this->refinementTableWidgetItemValue->setProperty("type", "Polyline");
		this->refinementTableWidgetItemValue->setProperty("number", p);

		connect(refinementTableWidgetItemValue, SIGNAL(valueChanged(double)), this, SLOT(refinementUpdate(double)));

		this->refinementTableWidgetItemValue->setDecimals(3);
		this->refinementTableWidgetItemValue->setValue(Model.Polylines[p].size / Model.scale);
		this->refinementTableWidget->setCellWidget(p + Model.Surfaces.length(), 1, this->refinementTableWidgetItemValue);
	}
}

void MainWindow::refinementUpdate(double value)
{
	QDoubleSpinBox *senderDSB = qobject_cast<QDoubleSpinBox *>(sender());
	if (senderDSB)
	{
		if (senderDSB->property("type") == "Surface")
			Model.Surfaces[senderDSB->property("number").toInt()].size = value * Model.scale;
		if (senderDSB->property("type") == "Polyline")
			Model.Polylines[senderDSB->property("number").toInt()].size = value * Model.scale;
	}
}

void MainWindow::preMeshGradientUpdate(double value)
{
	Model.preMeshGradient = value;
}

void MainWindow::meshGradientUpdate(double value)
{
	Model.meshGradient = value;
}

void MainWindow::material1dFill()
{
	this->material1dList->clear();

	for (int p = 0; p != Model.Polylines.length(); p++)
	{
		if (Model.Polylines[p].Type == "WELL")
			material1dListItem = new QListWidgetItem(QIcon(":/images/wells.png"), Model.Polylines[p].Name);
		this->material1dListItem->setCheckState(Qt::Unchecked);
		if (Model.Polylines[p].MaterialID >= 0)
			this->material1dListItem->setCheckState(Qt::Checked);

		connect(material1dList, SIGNAL(itemChanged(QListWidgetItem *)), this, SLOT(material1dUpdate(QListWidgetItem *)));

		this->material1dList->addItem(this->material1dListItem);
	}
}

void MainWindow::material1dUpdate(QListWidgetItem *item)
{
	Model.Polylines[item->listWidget()->row(item)].MaterialID = item->checkState() - 1;
	this->material3dIndexing();
}

void MainWindow::material2dFill()
{
	this->material2dList->clear();

	for (int s = 0; s != Model.Surfaces.length(); s++)
	{
		if (Model.Surfaces[s].Type == "UNIT")
			material2dListItem = new QListWidgetItem(QIcon(":/images/units.png"), Model.Surfaces[s].Name);
		if (Model.Surfaces[s].Type == "FAULT")
			material2dListItem = new QListWidgetItem(QIcon(":/images/faults.png"), Model.Surfaces[s].Name);
		if (Model.Surfaces[s].Type == "BORDER")
			material2dListItem = new QListWidgetItem(QIcon(":/images/borders.png"), Model.Surfaces[s].Name);
		this->material2dListItem->setCheckState(Qt::Unchecked);
		if (Model.Surfaces[s].MaterialID >= 0)
			this->material2dListItem->setCheckState(Qt::Checked);

		connect(material2dList, SIGNAL(itemChanged(QListWidgetItem *)), this, SLOT(material2dUpdate(QListWidgetItem *)));

		this->material2dList->addItem(this->material2dListItem);
	}
}

void MainWindow::material2dUpdate(QListWidgetItem *item)
{
	Model.Surfaces[item->listWidget()->row(item)].MaterialID = item->checkState() - 1;
	this->material3dIndexing();
}

void MainWindow::material3dFill()
{
	if (this->material3dList->currentRow() == -1)
		this->material3dList->blockSignals(true);
	Model.drawMats = true;
	this->material3dList->clear();

	for (int m = 0; m != Model.Mats.length(); m++)
	{
		this->material3dListItem = new QListWidgetItem(QIcon(":/images/tetrahedron.png"), "Material: " + QVariant(m).toString());
		this->material3dList->addItem(this->material3dListItem);
	}
	this->callMakeMats();
	this->material3dList->blockSignals(false);
}

void MainWindow::material3dFillLocation(int Material)
{
	if (this->material3dListLocation->currentRow() == -1)
		this->material3dListLocation->blockSignals(true);

	Model.drawMats = true;
	this->material3dListLocation->clear();

	if (Material < 0)
		return;
	for (int l = 0; l != Model.Mats[Material].Locations.length(); l++)
	{
		this->material3dListLocationItem = new QListWidgetItem(QIcon(":/images/tetrahedron.png"), "Location: " + QVariant(l).toString());
		this->material3dListLocation->addItem(this->material3dListLocationItem);
	}
	this->callMakeMats();
	this->material3dListLocation->blockSignals(false);
}

void MainWindow::material3dAdd()
{
	C_Material *Material = new C_Material;
	Model.Mats.append(*Material);
	this->material3dIndexing();
	this->material3dFill();
	this->material3dList->setCurrentRow(Model.Mats.length() - 1);
}

void MainWindow::material3dRemove()
{
	int Material = this->material3dList->currentRow();
	if (Material < 0)
		return;
	Model.Mats.removeAt(Material);
	this->material3dIndexing();
	this->material3dFill();
	if (Material == Model.Mats.length())
		this->material3dList->setCurrentRow(Material - 1);
	else
		this->material3dList->setCurrentRow(Material);
}

void MainWindow::material3dIndexing()
{
	int Material = Model.Mats.length();
	for (int s = 0; s != Model.Surfaces.length(); s++)
		if (Model.Surfaces[s].MaterialID >= 0)
			Model.Surfaces[s].MaterialID = Material++;
	for (int p = 0; p != Model.Polylines.length(); p++)
		if (Model.Polylines[p].MaterialID >= 0)
			Model.Polylines[p].MaterialID = Material++;
}

void MainWindow::material3dLocationAdd()
{
	int Material = this->material3dList->currentRow();
	if (Material < 0)
		return;
	Model.Mats[Material].Locations.append(C_Vector3D(0.0, 0.0, 0.0));
	this->material3dFillLocation(Material);
	this->material3dListLocation->setCurrentRow(Model.Mats[Material].Locations.length() - 1);
}

void MainWindow::material3dLocationRemove()
{
	int Material = this->material3dList->currentRow();
	int Location = this->material3dListLocation->currentRow();
	if (Material < 0)
		return;
	Model.Mats[Material].Locations.removeAt(Location);
	this->material3dFillLocation(Material);
	if (Location == Model.Mats[Material].Locations.length())
		this->material3dListLocation->setCurrentRow(Location - 1);
	else
		this->material3dListLocation->setCurrentRow(Location);
}

void MainWindow::selectionSetFlags(bool state)
{
	QPushButton *senderPB = qobject_cast<QPushButton *>(sender());

	const QVariant &type = senderPB->property("type");
	const QVariant &mode = senderPB->property("mode");

	if (state)
	{
		selectionState.tool = static_cast<SelectionTool>(type.toInt());
		selectionState.mode = static_cast<SelectionMode>(mode.toInt());
	}
	else
	{
		selectionState.reset();
	}

	glWidget->selectionState = selectionState;

	std::list<QPushButton *>::const_iterator it;
	for (it = selectionButtons.begin(); it != selectionButtons.end(); ++it)
		(*it)->setChecked(false);

	senderPB->setChecked(state);

	if (selectionState != POLYGON)
	{
		glWidget->polygonSelection.reset();
		glWidget->updateGL();
	}
}

void MainWindow::selection(bool possible)
{
	if (possible)
	{
		selectionState = selectionStatePrev;
		this->setSShowGBox(this->selectionComponent->currentText());
	}
	else
	{
		/* Save current selection state to be able to restore it once the
		 * selection component gets reactivated. */
		selectionStatePrev = selectionState;
		selectionState.reset();
		this->setSShowGBox("");
	}

	glWidget->selectionState = selectionState;
}

void MainWindow::findSelection(unsigned char R, unsigned char G, unsigned char B)
{
	for (int s1 = 0; s1 != Model.Surfaces.length(); s1++)
	{
		if (Model.Surfaces[s1].Name == this->selectionComponent->currentText())
		{
			for (int c1 = 0; c1 != Model.Surfaces[s1].Constraints.length(); c1++)
			{
				if (Model.Surfaces[s1].Constraints[c1].RGB[0] == R && Model.Surfaces[s1].Constraints[c1].RGB[1] == G && Model.Surfaces[s1].Constraints[c1].RGB[2] == B)
				{
					C_Line &cs = Model.Surfaces[s1].Constraints[c1];

					if (selectionState == MARK)
						Model.Surfaces[s1].Constraints[c1].Type = "SEGMENTS";
					else if (selectionState == UNMARK)
						Model.Surfaces[s1].Constraints[c1].Type = "UNDEFINED";
					else if (selectionState == HOLE)
						Model.Surfaces[s1].Constraints[c1].Type = "HOLES";
					else if (selectionState == INVERT)
					{
						if (cs.Type == "SEGMENTS" || cs.Type == "HOLES")
							cs.Type = "UNDEFINED";
						else if (cs.Type == "UNDEFINED")
							cs.Type = "SEGMENTS";
					}

					Model.Surfaces[s1].makeConstraints();

					for (int p = 0; p != Model.Polylines.length(); p++)
						for (int c2 = 0; c2 != Model.Polylines[p].Constraints.length(); c2++)
						{
							if (Model.Surfaces[s1].Constraints[c1].IsIdenticallyWith(Model.Polylines[p].Constraints[c2]))
							{
								Model.Polylines[p].Constraints[c2].Type = Model.Surfaces[s1].Constraints[c1].Type;
								Model.Polylines[p].makeConstraints();
							}
						}
					if (Model.Surfaces[s1].Constraints[c1].Type != "HOLES")
					{
						for (int s2 = 0; s2 != Model.Surfaces.length(); s2++)
							if (s1 != s2)
							{
								for (int c2 = 0; c2 != Model.Surfaces[s2].Constraints.length(); c2++)
								{
									if (((Model.Surfaces[s1].Constraints[c1].Type == "SEGMENTS" && Model.Surfaces[s2].Constraints[c2].Type != "HOLES") ||
											 Model.Surfaces[s1].Constraints[c1].Type == "UNDEFINED") &&
											Model.Surfaces[s1].Constraints[c1].IsIdenticallyWith(Model.Surfaces[s2].Constraints[c2]))
									{
										Model.Surfaces[s2].Constraints[c2].Type = Model.Surfaces[s1].Constraints[c1].Type;
										Model.Surfaces[s2].makeConstraints();
									}
								}
							}
					}
					glWidget->updateGL();
				}
			}
		}
	}
	for (int p = 0; p != Model.Polylines.length(); p++)
	{
		if (Model.Polylines[p].Name == this->selectionComponent->currentText())
		{
			for (int c1 = 0; c1 != Model.Polylines[p].Constraints.length(); c1++)
			{
				if (Model.Polylines[p].Constraints[c1].RGB[0] == R && Model.Polylines[p].Constraints[c1].RGB[1] == G && Model.Polylines[p].Constraints[c1].RGB[2] == B)
				{
					C_Line &cs = Model.Polylines[p].Constraints[c1];

					if (selectionState == MARK)
						Model.Polylines[p].Constraints[c1].Type = "SEGMENTS";
					else if (selectionState == UNMARK)
						Model.Polylines[p].Constraints[c1].Type = "UNDEFINED";
					else if (selectionState == INVERT)
					{
						if (cs.Type == "SEGMENTS" || cs.Type == "HOLES")
							cs.Type = "UNDEFINED";
						else if (cs.Type == "UNDEFINED")
							cs.Type = "SEGMENTS";
					}

					Model.Polylines[p].makeConstraints();

					for (int s = 0; s != Model.Surfaces.length(); s++)
						for (int c2 = 0; c2 != Model.Surfaces[s].Constraints.length(); c2++)
						{
							if (Model.Polylines[p].Constraints[c1].IsIdenticallyWith(Model.Surfaces[s].Constraints[c2]))
							{
								Model.Surfaces[s].Constraints[c2].Type = Model.Polylines[p].Constraints[c1].Type;
								Model.Surfaces[s].makeConstraints();
							}
						}
					glWidget->updateGL();
				}
			}
		}
	}
}

void MainWindow::settingCallByMenu()
{
	QAction *senderAction = qobject_cast<QAction *>(sender());
	QString phase = senderAction->property("phase").toString();
	this->settingShow(phase, senderAction->text());
}

void MainWindow::settingCallByTree(QTreeWidgetItem *item, int column)
{
	QString phase;
	QTreeWidgetItem *parent = item->parent();
	if (parent)
		phase = parent->text(0);
	this->settingShow(phase, item->text(column));
}

void MainWindow::settingShow(QString phase, QString setting)
{
	QRect treeWidgetGeometry;

	if (setting != "")
	{
		this->tetgenGBox->setHidden(true);
		this->refinementGBox->setHidden(true);
		this->interpolationGBox->setHidden(true);
		this->preMeshGradientGBox->setHidden(true);

		this->selectionGBox->setHidden(true);
		this->meshGradientGBox->setHidden(true);
		this->selection(false);

		this->material1dGBox->setHidden(true);
		this->material2dGBox->setHidden(true);
		this->material3dGBox->setHidden(true);

		this->errorGBox->setHidden(true);

		Model.drawMats = false;
		this->callMakeMats();
		if (setting == "Tetgen")
			this->tetgenGBox->setHidden(false);
		if (setting == "Refinement")
		{
			this->refinementGBox->setHidden(false);
			this->refinementFill();
		}
		if (setting == "Interpolation")
		{
			this->interpolationGBox->setHidden(false);
			this->interpolationFill();
		}

		if (setting == "Gradient" && phase == "PreMesh")
			this->preMeshGradientGBox->setHidden(false);

		bool showSelection = (setting == "Selection");
		this->selectionGBox->setHidden(!showSelection);
		this->selection(showSelection);

		if (setting == "Gradient" && phase == "Mesh")
			this->meshGradientGBox->setHidden(false);

		if (setting == "1D")
		{
			this->material1dGBox->setHidden(false);
			this->material1dFill();
		}
		if (setting == "2D")
		{
			this->material2dGBox->setHidden(false);
			this->material2dFill();
		}
		if (setting == "3D")
		{
			this->material3dGBox->setHidden(false);
			this->material3dFill();
		}

		if (setting == "Errors")
			this->errorGBox->setHidden(false);

		if (setting == "> Execute PreMesh <")
			this->preMesh();
		if (setting == "> Execute Mesh <")
			this->Mesh();

		if (setting == "PLC...")
			this->readPLC();

		if (setting == "mesh PLC")
			this->meshPLC();
	}
}

void MainWindow::preMesh()
{
	if (!this->threadPreMesh->isRunning())
	{
		clearMesh();
		this->threadPreMesh->setAttribute("PREMESHJOB");
		this->threadPreMesh->start();
	}
}

void MainWindow::Mesh()
{
	if (!this->threadMesh->isRunning())
	{
		clearMesh();
		this->threadMesh->setAttribute("MESHJOB");
		this->threadMesh->start();
	}
}

void MainWindow::readPLC()
{
	PLC.fileName = QFileDialog::getOpenFileName(this, tr("Select the PLC file to open"), PLC.filePath, tr("All Supported Files (*.obj)"));
	QApplication::processEvents();
	if (!PLC.fileName.isEmpty())
	{
		PLC.filePath = PLC.fileName.section("/", 0, -2);
		emit progress_append(">Start reading PLC file " + PLC.fileName + "...");
		QApplication::processEvents();
		PLC.readPLCFile();
		emit progress_append(">...finished");
	}
}

void MainWindow::meshPLC()
{
	emit progress_append(">Start tetrahedralization...");
	PLC.meshPLC(this->tetgenLineEdit->text());
	emit progress_append(">...finished");
}

void MainWindow::applyMaterialSelection()
{
	if (!this->threadMesh->isRunning())
	{
		this->threadMesh->setAttribute("MATERIAL_SELECTION");
		this->threadMesh->start();
	}
}

void MainWindow::fillTable()
{
	selectionTable->clear();
	selectionTable->setRowCount(0);
	selectionTable->setSortingEnabled(false);

	QStringList headerLabels;
	headerLabels << "Component"
							 << "Length"
							 << "Selected\nneighbors";
	selectionTable->setHorizontalHeaderLabels(headerLabels);
	selectionTable->horizontalHeader()->setStretchLastSection(true);
	selectionTable->horizontalHeader()->setMinimumSectionSize(125);

	int row = 0;
	for (int s = 0; s != Model.Surfaces.length(); s++)
	{
		C_Surface &surface = Model.Surfaces[s];
		for (int i = 0; i < surface.Constraints.length(); i++)
		{
			C_Line &line = surface.Constraints[i];
			if (line.Type != "UNDEFINED")
				continue;

			double pathLength = 0;
			for (int j = 0; j < line.Ns.size() - 1; ++j)
				pathLength += length(line.Ns[j + 1] - line.Ns[j]);

			int numNeighbors = 0;
			for (int k = 0; k < surface.Constraints.length(); k++)
			{
				double T = 1e-10;
				C_Line &other = surface.Constraints[k];
				if ((line.Ns.first().equals(other.Ns.first(), T) ||
						 line.Ns.first().equals(other.Ns.last(), T) ||
						 line.Ns.last().equals(other.Ns.first(), T) ||
						 line.Ns.last().equals(other.Ns.last(), T)) &&
						other.Type != "UNDEFINED")
					numNeighbors++;
			}

			selectionTable->insertRow(row);

			QTableWidgetItem *item;
			QVariant data;

			item = new QTableWidgetItem(surface.Name);
			data.setValue(static_cast<void *>(&line));
			item->setData(Qt::UserRole, data);
			selectionTable->setItem(row, 0, item);

			item = new QTableWidgetItem();
			item->setData(Qt::DisplayRole, pathLength);
			item->setTextAlignment(Qt::AlignRight);
			selectionTable->setItem(row, 1, item);

			item = new QTableWidgetItem();
			item->setData(Qt::DisplayRole, numNeighbors);
			item->setTextAlignment(Qt::AlignRight);
			selectionTable->setItem(row, 2, item);

			row++;
		}
	}

	selectionTable->setSortingEnabled(true);
	selectionTable->sortItems(1, Qt::AscendingOrder);
}

void MainWindow::onHighlightingChange()
{
	bool isValid;
	HighlightingMode mode;

	int index = highlightingCombo->currentIndex();
	QVariant data = highlightingCombo->itemData(index);
	mode = static_cast<HighlightingMode>(data.toInt(&isValid));

	if (!isValid)
		return;

	bool refillTable = !selectionTable->isVisible();

	Model.listMarkers = NULL;

	if (mode != HighlightingMode::SINGLE)
		selectionTable->hide();

	switch (mode)
	{
	case HighlightingMode::NONE:
		break;
	case HighlightingMode::COMPONENT:
		Model.makeMarkers(NULL, selectionComponent->currentText());
		break;
	case HighlightingMode::ALL:
		Model.makeMarkers();
		break;
	case HighlightingMode::SINGLE:
		if (refillTable)
			fillTable();
		selectionTable->show();
		break;
	}

	glWidget->updateGL();
}

void MainWindow::tableClicked(QTableWidgetItem *item)
{
	int row = item->row();

	QVariant vline = selectionTable->item(row, 0)->data(Qt::UserRole);
	C_Line *line = static_cast<C_Line *>(vline.value<void *>());

	QString sname = selectionTable->item(row, 0)->data(Qt::DisplayRole).toString();
	double length = selectionTable->item(row, 1)->data(Qt::DisplayRole).toDouble();

	int idx = selectionComponent->findText(sname);
	selectionComponent->setCurrentIndex(idx);

	C_Surface *surface = Model.findSurface(sname);
	if (!surface)
		return;

	for (int s = 0; s < Model.Surfaces.length(); ++s)
		Model.Surfaces[s].makeConstraints();

	surface->calculate_normal_vector();

	C_Vector3D center(0, 0, 0);
	for (int i = 0; i < line->Ns.size(); ++i)
		center += line->Ns[i];

	center /= line->Ns.size();
	center *= -1;

	Model.makeMarkers(line);

	/* Move to center of marker. Pass surface->normal_vector as the last
	 * argument to moveViewport() in order to establish a plan view. */
	double dist = std::min(glWidget->currentDist(), std::pow(2, length));
	glWidget->moveViewport(center, dist);
}

void MainWindow::clearErrorInfo()
{
	emit Model.ErrorInfoChanged("");
}

void MainWindow::updateErrorInfo(QString message)
{
	bool show = !message.isEmpty();
	errorWidget->setHidden(!show);

	clearErrorMarkers();

	if (show)
	{
		errorWidget->treeWidget()->selectionModel()->clearSelection();
		errorWidget->treeWidget()->collapseAll();

		errorLabel->setText(message);
		errorLabel->setStyleSheet("QLabel { color: red; }");
		fillErrorTable();

		errorWidget->setSelected(true);
		settingShow("Mesh", "Errors");
	}
	else
	{
		this->errorGBox->setHidden(true);
	}
}

void MainWindow::fillErrorTable()
{
	QStringList headerLabels;
	headerLabels << "Component";

	errorTable->clear();
	errorTable->setRowCount(0);
	errorTable->setColumnCount(1);
	errorTable->setSortingEnabled(false);
	errorTable->setHorizontalHeaderLabels(headerLabels);
	errorTable->horizontalHeader()->setStretchLastSection(true);
	errorTable->horizontalHeader()->setMinimumSectionSize(125);

	int row = 0;

	QList<SelfIntersection>::iterator it;
	for (it = Model.selfIntersections.begin();
			 it != Model.selfIntersections.end(); ++it)
	{
		QSet<const C_Surface *>::const_iterator sit;
		for (sit = it->second.cbegin(); sit != it->second.cend(); ++sit)
		{
			/* Add table entry. */
			errorTable->insertRow(row);

			QVariant data;
			QTableWidgetItem *item = new QTableWidgetItem((*sit)->Name);
			data.setValue(static_cast<void *>(&it->first));
			item->setData(Qt::UserRole, data);
			errorTable->setItem(row, 0, item);

			row++;
		}
	}

	errorTable->setSortingEnabled(true);
	errorTable->sortItems(0, Qt::AscendingOrder);
}

void MainWindow::errorTableClicked(QTableWidgetItem *item)
{
	int row = item->row();

	QVariant vtri = errorTable->item(row, 0)->data(Qt::UserRole);
	C_Triangle *tri = static_cast<C_Triangle *>(vtri.value<void *>());

	QString sname = errorTable->item(row, 0)->data(Qt::DisplayRole).toString();

	int idx = selectionComponent->findText(sname);
	selectionComponent->setCurrentIndex(idx);

	C_Surface *surface = Model.findSurface(sname);
	if (!surface)
		return;

	for (int s = 0; s < Model.Surfaces.length(); ++s)
	{
		Model.Surfaces[s].makeConstraints();
		Model.Surfaces[s].drawConstraints = false;
	}

	surface->calculate_normal_vector();
	surface->drawConstraints = true;

	C_Line line;
	line.Ns.append(*tri->Ns[0]);
	line.Ns.append(*tri->Ns[1]);
	line.Ns.append(*tri->Ns[2]);
	line.Ns.append(*tri->Ns[0]);

	C_Vector3D center(0, 0, 0);
	for (int i = 0; i < line.Ns.size(); ++i)
		center += line.Ns[i];

	center /= line.Ns.size();
	center *= -1;

	double lineLength = 0;
	for (int j = 0; j < line.Ns.size() - 1; ++j)
		lineLength += length(line.Ns[j + 1] - line.Ns[j]);

	Model.makeErrorMarkers(line);

	/* Move to center of marker. Pass surface->normal_vector as the last
	 * argument to moveViewport() in order to establish a plan view. */
	double dist = std::min(glWidget->currentDist(), std::pow(2, lineLength));
	glWidget->moveViewport(center, dist);
}

void MainWindow::clearErrorMarkers()
{
	for (int i = 0; i < Model.Surfaces.length(); ++i)
		Model.Surfaces[i].drawConstraints = false;

	Model.listErrorMarkers = NULL;
	glWidget->updateGL();
}

void MainWindow::FinishedRead()
{
	glWidget->resetView();
	// surfaces
	for (int s = 0; s != Model.Surfaces.length(); s++)
	{
		Model.Surfaces[s].clearScatteredData();
		Model.Surfaces[s].makeScatteredData();
		Model.Surfaces[s].makeConvexHull();
		Model.Surfaces[s].makeFaces();
		Model.Surfaces[s].makeEdges();
		Model.Surfaces[s].makeIntEdges();
		Model.Surfaces[s].makeIntVertices();
		Model.Surfaces[s].makeConstraints();
		glWidget->Model = &Model;
	}
	// polylines
	for (int p = 0; p != Model.Polylines.length(); p++)
	{
		Model.Polylines[p].makeScatteredData();
		Model.Polylines[p].makeEdges();
		Model.Polylines[p].makeVertices();
		Model.Polylines[p].makeIntVertices();
		Model.Polylines[p].makeConstraints();
		glWidget->Model = &Model;
	}
	/* Fill radius values. */
	preMeshGradientValue->setValue(Model.preMeshGradient);
	meshGradientValue->setValue(Model.meshGradient);
	this->FillNameCombos();
	if (Model.Mesh)
	{
		this->callMakeTets();
		this->matsTitle->setHidden(false);
		Model.drawTets = true;
	}
	else
	{
		this->matsTitle->setHidden(true);
		Model.drawTets = false;
	}

	clearErrorInfo();

	glWidget->updateGL();
}

void MainWindow::threadFinishedPreMesh()
{
	// surfaces
	for (int s = 0; s != Model.Surfaces.length(); s++)
	{
		Model.Surfaces[s].makeScatteredData();
		Model.Surfaces[s].makeConvexHull();
		Model.Surfaces[s].makeFaces();
		Model.Surfaces[s].makeEdges();
		Model.Surfaces[s].makeIntEdges();
		Model.Surfaces[s].makeIntVertices();
		Model.Surfaces[s].makeConstraints();
	}
	// polylines
	for (int p = 0; p != Model.Polylines.length(); p++)
	{
		Model.Polylines[p].makeScatteredData();
		Model.Polylines[p].makeEdges();
		Model.Polylines[p].makeVertices();
		Model.Polylines[p].makeIntVertices();
		Model.Polylines[p].makeConstraints();
	}
	glWidget->updateGL();
}

void MainWindow::threadFinishedMesh()
{
	// surfaces
	for (int s = 0; s != Model.Surfaces.length(); s++)
	{
		Model.Surfaces[s].makeFaces();
		Model.Surfaces[s].makeEdges();
		Model.Surfaces[s].makeConstraints();
	}
	// polylines
	for (int p = 0; p != Model.Polylines.length(); p++)
	{
		Model.Polylines[p].makeEdges();
		Model.Polylines[p].makeVertices();
		Model.Polylines[p].makeConstraints();
	}
	this->callMakeTets();
	this->FillNameCombos();

	if (Model.Mesh)
	{
		this->matsTitle->setHidden(false);
		Model.drawTets = true;
	}

	selectionMaterial->stopAnimation();

	glWidget->updateGL();
}

void MainWindow::callMakeMats()
{
	Model.makeMats(this->material3dList->currentRow(), this->material3dListLocation->currentRow());
	glWidget->update();
}

void MainWindow::callMakeTets()
{
	if (Model.Mesh)
		Model.makeTets(this->xCutEnable->isChecked(), this->xCutSlider->value(), this->xDirChange->isChecked(), this->yCutEnable->isChecked(), this->yCutSlider->value(), this->yDirChange->isChecked(), this->zCutEnable->isChecked(), this->zCutSlider->value(), this->zDirChange->isChecked());
	glWidget->update();
}

void MainWindow::updateInfo()
{
	infoMinX->setText(QString::number(Model.min.x(), 'f', 0));
	infoMaxX->setText(QString::number(Model.max.x(), 'f', 0));
	infoRanX->setText(QString::number((Model.max.x() - Model.min.x()), 'f', 0));
	infoMinY->setText(QString::number(Model.min.y(), 'f', 0));
	infoMaxY->setText(QString::number(Model.max.y(), 'f', 0));
	infoRanY->setText(QString::number((Model.max.y() - Model.min.y()), 'f', 0));
	infoMinZ->setText(QString::number(Model.min.z(), 'f', 0));
	infoMaxZ->setText(QString::number(Model.max.z(), 'f', 0));
	infoRanZ->setText(QString::number((Model.max.z() - Model.min.z()), 'f', 0));
	if (Model.Mesh)
		infoTetNumber->setText("no. Tetrahedrons: " + QString::number(Model.Mesh->numberoftetrahedra));
	else
		infoTetNumber->setText("no. Tetrahedrons: " + QString::number(0));
}

void MainWindow::replace(QString text)
{
	statusTE->moveCursor(QTextCursor::End, QTextCursor::MoveAnchor);
	statusTE->moveCursor(QTextCursor::StartOfLine, QTextCursor::KeepAnchor);
	statusTE->textCursor().insertText(text);
}

void MainWindow::printErrorMessage(QString message)
{
	const QColor &prevColor = statusTE->palette().text().color();
	statusTE->setTextColor("red");
	statusTE->append(message);
	statusTE->setTextColor(prevColor);
}

void MainWindow::onRotationCenterChange(int index)
{
	bool isValid;
	RotationMode mode;

	QComboBox *comboBox = qobject_cast<QComboBox *>(sender());
	mode = static_cast<RotationMode>(comboBox->itemData(index).toInt(&isValid));

	if (isValid)
	{
		QString itemText = comboBox->itemText(index);
		QString msgText = "Rotation mode changed to <i>" + itemText + "</i>";
		emit progress_append(msgText);

		glWidget->rotationMode = mode;
		glWidget->updateGL();
	}
}

void MainWindow::onSensitivityChange()
{
	QSlider *slider = qobject_cast<QSlider *>(sender());
	glWidget->zoomSensitivity = (double)slider->value() / 100.;
}
