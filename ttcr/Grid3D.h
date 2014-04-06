//
//  Grid.h
//  ttcr_u
//
//  Created by Bernard Giroux on 2013-01-10.
//  Copyright (c) 2013 Bernard Giroux. All rights reserved.
//

#ifndef ttcr_u_Grid3D_h
#define ttcr_u_Grid3D_h

#include "ttcr_t.h"

template<typename T1, typename T2>
class Grid3D {
public:
    virtual int raytrace(const std::vector<sxyz<T1>>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<sxyz<T1>>& Rx,
                         std::vector<T1>& traveltimes,
                         const size_t threadNo=0) const { return 0; }
	
	virtual int raytrace(const std::vector<sxyz<T1>>&,
						 const std::vector<T1>&,
						 const std::vector<const std::vector<sxyz<T1>>*>&,
						 std::vector<std::vector<T1>*>&,
						 const size_t=0) const { return 0; }
	
    virtual int raytrace(const std::vector<sxyz<T1>>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<sxyz<T1>>& Rx,
                         std::vector<T1>& traveltimes,
                         std::vector<std::vector<sxyz<T1>>>& r_data,
                         const size_t threadNo=0) const { return 0; }
    
    virtual int raytrace(const std::vector<sxyz<T1>>&,
						 const std::vector<T1>&,
						 const std::vector<const std::vector<sxyz<T1>>*>&,
						 std::vector<std::vector<T1>*>&,
						 std::vector<std::vector<std::vector<sxyz<T1>>>*>&,
						 const size_t=0) const { return 0; }
	
    virtual int setSlowness(const std::vector<T1>& s) { return 0; }
    
    virtual size_t getNumberOfNodes() const { return 0; }
	
	virtual void saveTT(const std::string &, const int, const size_t nt=0,
						const bool vtkFormat=0) const {}
	
	virtual T1 getXmin() const { return 0; }
	virtual T1 getXmax() const { return 0; }
	virtual T1 getYmin() const { return 0; }
	virtual T1 getYmax() const { return 0; }
	virtual T1 getZmin() const { return 0; }
	virtual T1 getZmax() const { return 0; }

#ifdef VTK
    virtual void saveModelVTU(const std::string &, const bool saveSlowness=true) const {}
#endif
};

#endif
