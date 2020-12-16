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

#include "mainwindow.h"
#include "commandline.h"


int
main(int argc, char *argv[])
{
	Q_INIT_RESOURCE(MeshIT);

	QApplication app(argc, argv);
	QApplication::setApplicationName("MeshIt");
	QApplication::setApplicationVersion("1.0");
	
	if (argc > 1)
	{
		QCommandLineParser parser;
		parser.setApplicationDescription("MeshIt - a 3D mesh generator for fractured reservoirs.");
		parser.addHelpOption();
		parser.addVersionOption();

		QCommandLineOption preMeshOption("p", QApplication::translate("main", "performs the premeshing and selection of all calculated intersection"));
		parser.addOption(preMeshOption);

		QCommandLineOption meshOption("m", QApplication::translate("main", "performs the meshing"));
		parser.addOption(meshOption);

		QCommandLineOption inputDirectoryOption(QStringList() << "i" << "input",
			QApplication::translate("main", "reads MeshIt input state file <directory>."),
			QApplication::translate("main", "directory"));
		parser.addOption(inputDirectoryOption);

		QCommandLineOption outputDirectoryOption(QStringList() << "o" << "output",
			QApplication::translate("main", "writes MeshIt output state file <directory>."),
			QApplication::translate("main", "directory"));
		parser.addOption(outputDirectoryOption);

		QCommandLineOption exportDirectoryOption(QStringList() << "export-vtu" << "exportvtu",
			QApplication::translate("main", "export MeshIt mesh to vtu <directory>."),
			QApplication::translate("main", "directory"));
		parser.addOption(exportDirectoryOption);

		/* Process the actual command line arguments given by the user */
		parser.process(app);
		C_CommandLine commandLine(&parser);
	}
	else
	{
		app.setWindowIcon(QIcon(":/images/app_logo.png"));
		app.setOrganizationName("GFZ Potsdam");
		app.setApplicationName("MeshIt");
		MainWindow mainWin;
		mainWin.show();
		return app.exec();
	}
}
