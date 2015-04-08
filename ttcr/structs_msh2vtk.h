//
//  structs_msh2vtk.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-10-18.
//  Copyright (c) 2014 Bernard Giroux. All rights reserved.
//

#ifndef ttcr_structs_msh2vtk_h
#define ttcr_structs_msh2vtk_h

#include <string>

struct input_parameters {
	bool rectilinear;
	bool crt;
	bool saveSlowness;
	bool saveReflectors;
	bool saveXYZ;
	int verbose;
	int nSecondary;
	double d;
	std::string mshFile;
	std::string vtkFile;
	std::string velFile;
	std::string sloFile;
	
	input_parameters() : rectilinear(false), crt(false), saveSlowness(false),
	saveReflectors(false), saveXYZ(false), verbose(0), nSecondary(0), d(0.0),
	mshFile(""), vtkFile(""), velFile(""), sloFile("") {}
};

#endif
