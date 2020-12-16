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

#ifndef _COMMANDLINE_H_
#define _COMMANDLINE_H_

#include <QtWidgets/QtWidgets>

class C_CommandLine
{
public:
public:
	C_CommandLine(QCommandLineParser *);
	~C_CommandLine();
	void runThreadPool(QString, int, int, int, int);
	void preMeshJob();
	void MeshJob();
};

class C_CmdTask : public QRunnable
{
public:
	C_CmdTask(C_CommandLine *commandLine, const QString& Attribute, const int& Object1, const int& Object2, const int& currentStep, const int& totalSteps) :
		commandLine(commandLine), Attribute_(Attribute), Object1_(Object1), Object2_(Object2), currentStep_(currentStep), totalSteps_(totalSteps)
	{};
protected:
	C_CommandLine *commandLine;
	void run()
	{
		commandLine->runThreadPool(Attribute_, Object1_, Object2_, currentStep_, totalSteps_);
	}
private:
	QString Attribute_;
	int Object1_;
	int Object2_;
	int currentStep_;
	int totalSteps_;
};

#endif	// _COMMANDLINE_H_


