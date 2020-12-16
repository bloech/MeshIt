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

#ifndef _MAINWINDOW_H_
#define _MAINWINDOW_H_

#include <list>

#include <QtWidgets/QtWidgets>
#include "core.h"

QT_BEGIN_NAMESPACE
class QAction;
class QMenu;
class QPlainTextEdit;
class QScrollArea;
QT_END_NAMESPACE

class GLWidget;
class C_Thread;

class AnimatedButton: public QPushButton
{
	Q_OBJECT

private:
	QMovie *animation;

private slots:
	void updateAnimation(int);

public slots:
	void startAnimation();
	void stopAnimation();

public:
	explicit AnimatedButton(const QString &text, const QString &icon,
	                        QWidget *parent = Q_NULLPTR);
	~AnimatedButton();
};

class C_TableViewDelegate : public QItemDelegate
{
	Q_OBJECT
public:
	C_TableViewDelegate(QObject *parent = 0);
	QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const;
};

class MainWindow : public QMainWindow
{
	Q_OBJECT
signals:
	void progress_append(QString);
	void progress_replace(QString);
	void material3DLocationChanged(QString, double);

public:
	MainWindow();
	void runThread(QString);
	void runThreadPool(QString, int, int, int, int);

protected:
	void closeEvent(QCloseEvent *);
	void keyPressEvent(QKeyEvent *);

private slots:
	/*slots to open and save routines*/
	void open();
	void save();
	void saveAs();
	void saveScreen();
	/*slots to import/export routines*/
	void importGoCad();
	void exportVTU3D();
	void exportVTU2D();
	void exportFeFlow();
	void exportOGS();
	void exportTIN();
	void exportCOMSOL();
	void exportABAQUS();
	void exportEXODUS();
	/*slots to basic routines - add geometric objects and deleting*/
	void addUnit();
	void addFault();
	void addBorder();
	void addWell();
	void deleteSurface();
	/*additional slots - visualization, combo boxes etc ...*/
	void about();
	void reset();
	void viewAxis();
	// visualize combo boxes
	void FillNameCombos();
	// combo boxes units
	void setUShowGBox(QString);
	void setUMeshes();
	// combo boxes faults
	void setFShowGBox(QString);
	void setFMeshes();
	// combo boxes borders
	void setBShowGBox(QString);
	void setBMeshes();
	// combo boxes wells
	void setWShowGBox(QString);
	void setWMeshes();
	// combo boxes materials
	void setMShowGBox(QString);
	void setMMeshes();
	// combo boxes selections
	void setSShowGBox(QString);
	// set interpolation method
	void interpolationSetMethod(QString);
	void interpolationFill();
	// fill material sliders and cut (x-, y-, z-) location
	void material3dFillValue(int);
	void material3dSetLocationFromDSpinBox(double);
	void material3dSetLocationFromSlider(int);
	// export angel
	void ExportRotationAngelUpdate(double);
	// refinement
	void refinementFill();
	void refinementUpdate(double);
	// materials
	void material1dFill();
	void material1dUpdate(QListWidgetItem *);
	void material2dFill();
	void material2dUpdate(QListWidgetItem *);
	void material3dFill();
	void material3dFillLocation(int Material);
	void material3dAdd();
	void material3dRemove();
	void material3dIndexing();
	void material3dLocationAdd();
	void material3dLocationRemove();
	// selections
	void selectionSetFlags(bool);
	void selection(bool);
	void findSelection(unsigned char, unsigned char, unsigned char);
	// gradient
	void preMeshGradientUpdate(double);
	void meshGradientUpdate(double);
	// setting
	void settingCallByMenu();
	void settingCallByTree(QTreeWidgetItem *, int);
	void settingShow(QString phase, QString setting);
	// slots to basic geometric operations
	void preMesh();
	void Mesh();
	void applyMaterialSelection();
	void fillTable();
	void onHighlightingChange();
	void tableClicked(QTableWidgetItem *item);
	void clearErrorInfo();
	void updateErrorInfo(QString message);
	void fillErrorTable();
	void errorTableClicked(QTableWidgetItem *item);
	void clearErrorMarkers();
	// wait and stop the multithreading
	void FinishedRead();
	void threadFinishedPreMesh();
	void threadFinishedMesh();
	// fill materials and tets
	void callMakeMats();
	void callMakeTets();
	// update info box
	void updateInfo();
	// updates ranges and stepping vor Sliders and Doublespinboxes
	void replace(QString);
	// Error output.
	void printErrorMessage(QString message);
	// update rotation center
	void onRotationCenterChange(int index);
	// update zoom sensitivity
	void onSensitivityChange();

private:
	void createActions();
	void createMenus();
	void createDockWindows();
	void createStatusDock();
	void createViewDock();
	void createMeshQualityDock();
	void readSettings();
	void writeSettings();
	void preMeshJob();
	void MeshJob();
	void materialSelectionJob();
	void clearMesh();

	QWidget *centralWidget;
	QScrollArea *glWidgetArea;
	GLWidget *glWidget;
	/*list of QMenus*/
	QMenu *fileMenu;
	QMenu *fileMenuImport;
	QMenu *fileMenuExport;
	QMenu *fileMenuExport3D;
	QMenu *fileMenuExport2D;
	QMenu *fileMenuAdd;
	QMenu *fileMenuDelete;
	QMenu *editMenu;
	QMenu *editMenuMaterial;
	QMenu *viewMenu;
	QMenu *helpMenu;
	/*list of QActions*/
	QAction *openAct;
	QAction *saveAct;
	QAction *saveAsAct;
	QAction *saveScreenAct;
	QAction *importGoCadAct;
	QAction *exportFeFlowAct;
	QAction *exportOGSAct;
	QAction *exportCOMSOLAct;
	QAction *exportABAQUSAct;
	QAction *exportEXODUSAct;
	QAction *exportVTU3DAct;
	QAction *exportTINAct;
	QAction *exportVTU2DAct;
	QAction *addUnitAct;
	QAction *addFaultAct;
	QAction *addBorderAct;
	QAction *addWellAct;
	QAction *deleteSurfaceAct;
	QAction *editRefinement;
	QAction *editInterpolation;
	QAction *preMeshEditGradient;
	QAction *editPreMesh;
	QAction *editSelection;
	QAction *meshEditGradient;
	QAction *editMaterial1d;
	QAction *editMaterial2d;
	QAction *editMaterial3d;
	QAction *editTetgen;
	QAction *editMesh;
	QAction *exitAct;
	QAction *viewAxisAct;
	QAction *aboutAct;
	/*status Dock*/
	QDockWidget *statusDock;
	QTextEdit *statusTE;
	/*view Dock*/
	QDockWidget *viewDock;
	QWidget *viewWidget;
	QVBoxLayout *viewVBox;
	QTreeWidget *viewTree;
	QTreeWidgetItem *unitsTitle;
	QTreeWidgetItem *unitsNames;
	QComboBox *unitsNamesCB;
	QTreeWidgetItem *unitsScatteredDataPoints;
	QTreeWidgetItem *unitsConvexHull;
	QTreeWidgetItem *unitsSurface;
	QTreeWidgetItem *unitsFaces;
	QTreeWidgetItem *unitsEdges;
	QTreeWidgetItem *unitsIntersection;
	QTreeWidgetItem *unitsIntersectionEdges;
	QTreeWidgetItem *unitsIntersectionVertices;
	QTreeWidgetItem *faultsTitle;
	QTreeWidgetItem *faultsNames;
	QComboBox *faultsNamesCB;
	QTreeWidgetItem *faultsScatteredDataPoints;
	QTreeWidgetItem *faultsConvexHull;
	QTreeWidgetItem *faultsSurface;
	QTreeWidgetItem *faultsFaces;
	QTreeWidgetItem *faultsEdges;
	QTreeWidgetItem *faultsIntersection;
	QTreeWidgetItem *faultsIntersectionEdges;
	QTreeWidgetItem *faultsIntersectionVertices;
	QTreeWidgetItem *bordersTitle;
	QTreeWidgetItem *bordersNames;
	QComboBox *bordersNamesCB;
	QTreeWidgetItem *bordersScatteredDataPoints;
	QTreeWidgetItem *bordersConvexHull;
	QTreeWidgetItem *bordersSurface;
	QTreeWidgetItem *bordersFaces;
	QTreeWidgetItem *bordersEdges;
	QTreeWidgetItem *bordersIntersection;
	QTreeWidgetItem *bordersIntersectionEdges;
	QTreeWidgetItem *bordersIntersectionVertices;
	QTreeWidgetItem *wellsTitle;
	QTreeWidgetItem *wellsNames;
	QComboBox *wellsNamesCB;
	QTreeWidgetItem *wellsScatteredDataPoints;
	QTreeWidgetItem *wellsPath;
	QTreeWidgetItem *wellsEdges;
	QTreeWidgetItem *wellsVertices;
	QTreeWidgetItem *wellsIntersectionVertices;
	QTreeWidgetItem *matsTitle;
	QTreeWidgetItem *matsNames;
	QComboBox *matsNamesCB;
	QTreeWidgetItem *matsFaces;
	QTreeWidgetItem *matsEdges;
	/*MeshVies Dock*/
	QGroupBox *tetViewGBox;
	QGridLayout *tetViewGrid;
	QCheckBox *xCutEnable;
	QCheckBox *yCutEnable;
	QCheckBox *zCutEnable;
	QCheckBox *xDirChange;
	QCheckBox *yDirChange;
	QCheckBox *zDirChange;
	QSlider *xCutSlider;
	QSlider *yCutSlider;
	QSlider *zCutSlider;
	/* Rotation center view. */
	QComboBox *rotationCombo;
	/*InfoVies Dock*/
	QGroupBox *infoGBox;
	QGridLayout *infoGrid;
	QLabel *infoX;
	QLabel *infoY;
	QLabel *infoZ;
	QLabel *infoMin;
	QLabel *infoMax;
	QLabel *infoRan;
	QLabel *infoMinX;
	QLabel *infoMaxX;
	QLabel *infoRanX;
	QLabel *infoMinY;
	QLabel *infoMaxY;
	QLabel *infoRanY;
	QLabel *infoMinZ;
	QLabel *infoMaxZ;
	QLabel *infoRanZ;
	QLabel *infoTetNumber;
	QSpacerItem *viewSpacer;
	/*MeshQuality Dock*/
	QDockWidget *meshDock;
	QWidget *meshWidget;
	QVBoxLayout *meshVBox;
	QTreeWidget *meshTree;
	QTreeWidgetItem *meshPreMesh;
	QTreeWidgetItem *meshPreMeshRefinement;
	QTreeWidgetItem *meshPreMeshInterpolation;
	QTreeWidgetItem *meshPreMeshGradient;
	QTreeWidgetItem *meshPreMeshExecute;
	QTreeWidgetItem *meshMesh;
	QTreeWidgetItem *meshMeshSelection;
	QTreeWidgetItem *meshMeshGradient;
	QTreeWidgetItem *meshMeshMaterial;
	QTreeWidgetItem *meshMeshMaterial1d;
	QTreeWidgetItem *meshMeshMaterial2d;
	QTreeWidgetItem *meshMeshMaterial3d;
	QTreeWidgetItem *meshMeshTetgen;
	QTreeWidgetItem *meshMeshExecute;
	QSpacerItem *meshSpacer;
	/*Tetgen Dock*/
	QLineEdit *tetgenLineEdit;
	QGroupBox *tetgenGBox;
	QGridLayout *tetgenGrid;
	/*Refine Dock*/
	QGroupBox *refinementGBox;
	QGridLayout *refinementGrid;
	QTableWidget *refinementTableWidget;
	QTableWidgetItem *refinementTableWidgetItemObject;
	QDoubleSpinBox *refinementTableWidgetItemValue;
	/*Interpolation Dock*/
	QGroupBox *interpolationGBox;
	QGridLayout *interpolationGrid;
	QComboBox *interpolationMethod;

	/* Gradient */
	QGroupBox *preMeshGradientGBox;
	QGridLayout *preMeshGradientGrid;
	QDoubleSpinBox *preMeshGradientValue;

	QGroupBox *meshGradientGBox;
	QGridLayout *meshGradientGrid;
	QDoubleSpinBox *meshGradientValue;

	/*Selection Dock*/
	QGroupBox *selectionGBox;
	QGridLayout *selectionGrid;
	QComboBox *selectionComponent;
	QLabel *selectionSingle;
	QLabel *selectionMulti;
	QLabel *selectionAll;
	QLabel *selectionMark;
	QLabel *selectionUnmark;
	QLabel *selectionHole;

	std::list<QPushButton*> selectionButtons;
	QPushButton *selectionMarkSingle;
	QPushButton *selectionMarkMulti;
	QPushButton *selectionMarkAll;
	QPushButton *selectionUnmarkSingle;
	QPushButton *selectionUnmarkMulti;
	QPushButton *selectionUnmarkAll;
	QPushButton *selectionHoleSingle;
	QPushButton *selectionHoleMulti;
	QPushButton *selectionHoleAll;

	/* Bucket selection. */
	QPushButton *selectionMarkBucket;
	QPushButton *selectionUnmarkBucket;
	QPushButton *selectionHoleBucket;

	/* Polygon selection. */
	QPushButton *selectionMarkPolygon;
	QPushButton *selectionUnmarkPolygon;
	QPushButton *selectionHolePolygon;

	/* Inversion tool. */
	QPushButton *selectionInvert;

	/* Material based selection tool. */
	AnimatedButton *selectionMaterial;

	QComboBox *highlightingCombo;
	QTableWidget *selectionTable;

	/*Material1D Dock*/
	QGroupBox *material1dGBox;
	QGridLayout *material1dGrid;
	QListWidget *material1dList;
	QListWidgetItem  *material1dListItem;
	/*Material2D Dock*/
	QGroupBox *material2dGBox;
	QGridLayout *material2dGrid;
	QListWidget *material2dList;
	QListWidgetItem  *material2dListItem;
	/*Material3D Dock*/
	QGroupBox *material3dGBox;
	QGridLayout *material3dGrid;
	QListWidget *material3dList;
	QListWidgetItem  *material3dListItem;
	QListWidget *material3dListLocation;
	QListWidgetItem *material3dListLocationItem;
	QLabel *material3dAddRemoveMaterialLabel;
	QLabel *material3dAddRemoveLocationLabel;
	QLabel *material3dXLabel;
	QLabel *material3dYLabel;
	QLabel *material3dZLabel;
	QPushButton *material3dAddMaterialPushButton;
	QPushButton *material3dRemoveMaterialPushButton;
	QPushButton *material3dAddLocationPushButton;
	QPushButton *material3dRemoveLocationPushButton;
	QSlider *material3dXSlider;
	QSlider *material3dYSlider;
	QSlider *material3dZSlider;
	QDoubleSpinBox *material3dXValue;
	QDoubleSpinBox *material3dYValue;
	QDoubleSpinBox *material3dZValue;

	/* Error reporting. */
	QTreeWidgetItem *errorWidget;
	QGroupBox *errorGBox;
	QVBoxLayout *errorLayout;

	QLabel *errorLabel;
	QTableWidget *errorTable;

	/*C_Threads*/
	C_Thread *threadPreMesh;
	C_Thread *threadMesh;
	QFont boldFont;

	/* State */
	SelectionState selectionState;
	SelectionState selectionStatePrev;
};

class C_Task : public QRunnable
{
public:
	C_Task(MainWindow *mainWindow, const QString &Attribute, const int &Object1, const int &Object2, const int &currentStep, const int &totalSteps) :
		mainWindow(mainWindow), Attribute_(Attribute), Object1_(Object1), Object2_(Object2), currentStep_(currentStep), totalSteps_(totalSteps)
	{};

protected:
	MainWindow *mainWindow;
	void run()
	{
		mainWindow->runThreadPool(Attribute_, Object1_, Object2_, currentStep_, totalSteps_);
	}

private:
	QString Attribute_;
	int Object1_;
	int Object2_;
	int currentStep_;
	int totalSteps_;
};

class C_Thread : public QThread
{
	Q_OBJECT
public:

	C_Thread(MainWindow *mainWindow) : mainWindow(mainWindow)
	{};
	virtual ~C_Thread()
	{};
	void run()
	{
		mainWindow->runThread(Attribute_);
	}
	void setAttribute(const QString & Attribute)
	{
		Attribute_ = Attribute;
	}

protected:
	MainWindow *mainWindow;

private:
	QString Type_;
	QString Attribute_;
};

#endif	// _MAINWINDOW_H_


