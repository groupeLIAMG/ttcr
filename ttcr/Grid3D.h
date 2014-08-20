// Grid3D is the class template of an object containing the 3D grid and the function
// to perform raytracing
//
//  Grid.h
//  ttcr_u
//
//  Created by Bernard Giroux on 2013-01-10.
//  Copyright (c) 2013 Bernard Giroux. All rights reserved.
//

/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef ttcr_u_Grid3D_h
#define ttcr_u_Grid3D_h

#include "ttcr_t.h"

template<typename T1, typename T2>
class Grid3D {
public:

	virtual ~Grid3D() {}

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
	
	virtual void setSourceRadius(const double) {}
    
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
