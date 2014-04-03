//
//  structs_spmrt.h
//  ttcr_u
//
//  Created by Bernard Giroux on 2012-11-19.
//  Copyright (c) 2012 Bernard Giroux. All rights reserved.
//

#ifndef ttcr_u_structs_spmrt_h
#define ttcr_u_structs_spmrt_h

#include <map>
#include <string>

enum raytracing_method { SHORTEST_PATH, FAST_MARCHING, FAST_SWEEPING };

struct input_parameters {
	uint32_t nn[3];
	int nt;
	int verbose;
	int order;                    // order of l metric
	int nitermax;
    bool inverseDistance;
	bool singlePrecision;
	bool saveRaypaths;
	bool saveModelVTK;
    bool saveGridTT;
	bool time;
	bool processReflectors;
	double epsilon;
    raytracing_method method;
	std::string basename;
	std::string modelfile;
	std::string velfile;
    std::string slofile;
	std::string rcvfile;
    std::vector<std::string> srcfiles;
	
	input_parameters() : nn(), nt(0), verbose(0), order(2), nitermax(20),
    inverseDistance(false),	singlePrecision(false), saveRaypaths(false),
    saveModelVTK(false), saveGridTT(false), time(false), processReflectors(false),
	epsilon(1.e-15), method(SHORTEST_PATH), basename(),
    modelfile(), velfile(), slofile(), rcvfile(), srcfiles() {}

};

#endif
