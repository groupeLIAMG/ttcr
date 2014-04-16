//
//  Grid2Dui.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-04-11.
//  Copyright (c) 2014 Bernard Giroux. All rights reserved.
//

#ifndef ttcr_Grid2Dui_h
#define ttcr_Grid2Dui_h

#include "Grid2D.h"
#include "Interpolator.h"

template<typename T1, typename T2, typename NODE, typename S>
class Grid2Dui : public Grid2D<T1,T2,S> {
public:
	Grid2Dui(const std::vector<S>& no,
             const std::vector<triangleElem<T2>>& tri,
             const size_t nt=1) :
	nThreads(nt),
	nPrimary(static_cast<T2>(no.size())),
	nodes(std::vector<NODE>(no.size(), NODE(nt))),
	neighbors(std::vector<std::vector<T2>>(tri.size())),
	triangles(tri)
	{}
    
    void setSlowness(const T1 s) {
        for ( size_t n=0; n<nodes.size(); ++n ) {
            nodes[n].setSlowness( s );
        }
    }
    
	int setSlowness(const T1 *s, const size_t ns) {
        if ( nodes.size() != ns ) {
            std::cerr << "Error: slowness vectors of incompatible size.";
            return 1;
        }
        for ( size_t n=0; n<nodes.size(); ++n ) {
            nodes[n].setSlowness( s[n] );
        }
        return 0;
    }
    
	virtual int setSlowness(const std::vector<T1>& s) {
        if ( nodes.size() != s.size() ) {
            std::cerr << "Error: slowness vectors of incompatible size.";
            return 1;
        }
        for ( size_t n=0; n<nodes.size(); ++n ) {
            nodes[n].setNodeSlowness( s[n] );
        }
        return 0;
    }
    
    void setTT(const T1 tt, const size_t nn, const size_t nt=0) {
		nodes[nn].setTT(tt, nt);
	}
    
    virtual int raytrace(const std::vector<S>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<S>& Rx,
                         std::vector<T1>& traveltimes,
                         const size_t threadNo=0) const { return 0; }
    
    virtual int raytrace(const std::vector<S>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<const std::vector<S>*>& Rx,
                         std::vector<std::vector<T1>*>& traveltimes,
                         const size_t threadNo=0) const { return 0; }
    
    virtual int raytrace(const std::vector<S>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<S>& Rx,
                         std::vector<T1>& traveltimes,
                         std::vector<std::vector<S>>& r_data,
                         const size_t threadNo=0) const { return 0; }
    
    virtual int raytrace(const std::vector<S>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<const std::vector<S>*>& Rx,
                         std::vector<std::vector<T1>*>& traveltimes,
                         std::vector<std::vector<std::vector<S>>*>& r_data,
                         const size_t threadNo=0) const { return 0; }
	
	size_t getNumberOfNodes() const { return nodes.size(); }
    
	T1 getXmin() const {
		T1 xmin = nodes[0].getX();
		for ( auto it=nodes.begin(); it!=nodes.end(); ++it )
			xmin = xmin<it->getX() ? xmin : it->getX();
		return xmin;
	}
	T1 getXmax() const {
		T1 xmax = nodes[0].getX();
		for ( auto it=nodes.begin(); it!=nodes.end(); ++it )
			xmax = xmax>it->getX() ? xmax : it->getX();
		return xmax;
	}
	T1 getZmin() const {
		T1 zmin = nodes[0].getZ();
		for ( auto it=nodes.begin(); it!=nodes.end(); ++it )
			zmin = zmin<it->getZ() ? zmin : it->getZ();
		return zmin;
	}
	T1 getZmax() const {
		T1 zmax = nodes[0].getZ();
		for ( auto it=nodes.begin(); it!=nodes.end(); ++it )
			zmax = zmax>it->getZ() ? zmax : it->getZ();
		return zmax;
	}
    
	void saveTT(const std::string &, const int, const size_t nt=0,
                const bool vtkFormat=0) const;
	
#ifdef VTK
    void saveModelVTU(const std::string &, const bool saveSlowness=true) const;
    void saveModelVTR(const std::string &, const double*,
					  const bool saveSlowness=true) const;
#endif
	
protected:
	const size_t nThreads;
	T2 nPrimary;
	mutable std::vector<NODE> nodes;
    std::vector<std::vector<T2>> neighbors;  // nodes common to a cell
	std::vector<triangleElem<T2>> triangles;
    
	void buildGridNeighbors() {
		// Index the neighbors nodes of each cell
		for ( T2 n=0; n<nodes.size(); ++n ) {
			for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
				neighbors[ nodes[n].getOwners()[n2] ].push_back(n);
			}
		}
	}
	
    T1 computeDt(const NODE& source, const NODE& node) const {
		return (node.getNodeSlowness()+source.getNodeSlowness())/2 * source.getDistance( node );
	}
    
	T1 computeDt(const NODE& source, const S& node, T1 slo) const {
		return (slo+source.getNodeSlowness())/2 * source.getDistance( node );
	}
    
    T1 computeSlowness(const S& Rx ) const;
	
	T2 getCellNo(const S& pt) const {
		for ( T2 n=0; n<triangles.size(); ++n ) {
			if ( insideTriangle(pt, n) ) {
				return n;
			}
		}
		return -1;
	}
	
	T1 getTraveltime(const S& Rx,
					 const std::vector<NODE>& nodes,
					 const size_t threadNo) const;
	
	T1 getTraveltime(const S& Rx,
                     const std::vector<NODE>& nodes,
					 T2& nodeParentRx,
                     T2& cellParentRx,
					 const size_t threadNo) const;
	
	
	int check_pts(const std::vector<sxz<T1>>&) const;
	int check_pts(const std::vector<sxyz<T1>>&) const;
	
	bool insideTriangle(const sxz<T1>&, const T2) const;
	bool insideTriangle(const sxyz<T1>&, const T2) const;
    
};

template<typename T1, typename T2, typename NODE, typename S>
T1 Grid2Dui<T1,T2,NODE,S>::getTraveltime(const S& Rx,
                                         const std::vector<NODE>& nodes,
                                         const size_t threadNo) const {
    
    for ( size_t nn=0; nn<nodes.size(); ++nn ) {
        if ( nodes[nn] == Rx ) {
            return nodes[nn].getTT(threadNo);
        }
    }
    //If Rx is not on a node:
    T1 slo = computeSlowness( Rx );
    
    T2 cellNo = this->getCellNo( Rx );
    T2 neibNo = neighbors[cellNo][0];
    T1 dt = computeDt(nodes[neibNo], Rx, slo);
    
    T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
    for ( size_t k=1; k< neighbors[cellNo].size(); ++k ) {
        neibNo = neighbors[cellNo][k];
        dt = computeDt(nodes[neibNo], Rx, slo);
        if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
            traveltime =  nodes[neibNo].getTT(threadNo)+dt;
        }
    }
    return traveltime;
}


template<typename T1, typename T2, typename NODE, typename S>
T1 Grid2Dui<T1,T2,NODE,S>::getTraveltime(const S& Rx,
                                         const std::vector<NODE>& nodes,
                                         T2& nodeParentRx, T2& cellParentRx,
                                         const size_t threadNo) const {
    
    for ( size_t nn=0; nn<nodes.size(); ++nn ) {
        if ( nodes[nn] == Rx ) {
            nodeParentRx = nodes[nn].getNodeParent(threadNo);
            cellParentRx = nodes[nn].getCellParent(threadNo);
            return nodes[nn].getTT(threadNo);
        }
    }
    //If Rx is not on a node:
    T1 slo = computeSlowness( Rx );
    
    T2 cellNo = getCellNo( Rx );
    T2 neibNo = neighbors[cellNo][0];
    T1 dt = computeDt(nodes[neibNo], Rx, slo);
    
    T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
    nodeParentRx = neibNo;
    cellParentRx = cellNo;
    for ( size_t k=1; k< neighbors[cellNo].size(); ++k ) {
        neibNo = neighbors[cellNo][k];
        dt = computeDt(nodes[neibNo], Rx, cellNo);
        if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
            traveltime =  nodes[neibNo].getTT(threadNo)+dt;
            nodeParentRx = neibNo;
        }
    }
    return traveltime;
}

template<typename T1, typename T2, typename NODE, typename S>
T1 Grid2Dui<T1,T2,NODE,S>::computeSlowness( const S& Rx ) const {
    
    //Calculate the slowness of any point that is not on a node
    
    T2 cellNo = this->getCellNo( Rx );
    
    //We calculate the Slowness at the point
    std::vector<T2> list;
    
    for (size_t n3=0; n3 < neighbors[ cellNo ].size(); n3++){
        if ( nodes[neighbors[ cellNo ][n3] ].getPrimary() == 5 ){
            list.push_back(neighbors[ cellNo ][n3]);
        }
    }
    
    std::vector<size_t>::iterator it;
    
    std::vector<NODE*> interpNodes;
    
    for ( size_t nn=0; nn<list.size(); ++nn )
        interpNodes.push_back( &(nodes[list[nn] ]) );
    
    return Interpolator<T1>::inverseDistance( Rx, interpNodes );

}

template<typename T1, typename T2, typename NODE, typename S>
int Grid2Dui<T1,T2,NODE,S>::check_pts(const std::vector<sxz<T1>>& pts) const {
	
    for (size_t n=0; n<pts.size(); ++n) {
		bool found = false;
		// check first if point is on a node
		for ( T2 nt=0; nt<nodes.size(); ++nt ) {
			if ( nodes[nt] == pts[n]) {
				found = true;
				break;
			}
		}
		if ( found == false ) {
            for ( T2 nt=0; nt<triangles.size(); ++nt ) {
                if ( insideTriangle(pts[n], nt) ) {
                    found = true;
                }
            }
        }
		if ( found == false ) {
			std::cerr << "Error: point no " << (n+1)
            << " outside the grid.\n";
			return 1;
		}
    }
    return 0;
}

template<typename T1, typename T2, typename NODE, typename S>
int Grid2Dui<T1,T2,NODE,S>::check_pts(const std::vector<sxyz<T1>>& pts) const {
	
    for (size_t n=0; n<pts.size(); ++n) {
		bool found = false;
		// check first if point is on a node
		for ( T2 nt=0; nt<nodes.size(); ++nt ) {
			if ( nodes[nt] == pts[n]) {
				found = true;
				break;
			}
		}
		if ( found == false ) {
            for ( T2 nt=0; nt<triangles.size(); ++nt ) {
                if ( insideTriangle(pts[n], nt) ) {
                    found = true;
                }
            }
        }
		if ( found == false ) {
			std::cerr << "Error: point no " << (n+1)
            << " outside the grid.\n";
			return 1;
		}
    }
    return 0;
}


template<typename T1, typename T2, typename NODE, typename S>
bool Grid2Dui<T1,T2,NODE,S>::insideTriangle(const sxz<T1>& v, const T2 nt) const {
	
	
	// from http://mathworld.wolfram.com/TriangleInterior.html
	
	sxz<T1> v0 = { nodes[ triangles[nt].i[0] ].getX(),
		nodes[ triangles[nt].i[0] ].getZ() };
	sxz<T1> v1 = { nodes[ triangles[nt].i[1] ].getX()-v0.x,
		nodes[ triangles[nt].i[1] ].getZ()-v0.z };
	sxz<T1> v2 = { nodes[ triangles[nt].i[2] ].getX()-v0.x,
		nodes[ triangles[nt].i[2] ].getZ()-v0.z };
	
	T1 invDenom = 1. / det(v1, v2);
	T1 a = (det(v, v2) - det(v0, v2)) * invDenom;
	T1 b = -(det(v, v1) - det(v0, v1)) * invDenom;
	return (a >= 0.) && (b >= 0.) && (a + b < 1.);
}

template<typename T1, typename T2, typename NODE, typename S>
bool Grid2Dui<T1,T2,NODE,S>::insideTriangle(const sxyz<T1>& p, const T2 nt) const {
	
	
    
	sxyz<T1> a = { nodes[ triangles[nt].i[0] ].getX(),
		nodes[ triangles[nt].i[0] ].getY(),
		nodes[ triangles[nt].i[0] ].getZ() };
	sxyz<T1> b = { nodes[ triangles[nt].i[1] ].getX(),
		nodes[ triangles[nt].i[1] ].getY(),
		nodes[ triangles[nt].i[1] ].getZ() };
	sxyz<T1> c = { nodes[ triangles[nt].i[2] ].getX(),
		nodes[ triangles[nt].i[2] ].getY(),
		nodes[ triangles[nt].i[2] ].getZ() };
	
	// Translate point and triangle so that point lies at origin
	a -= p; b -= p; c -= p;
	// Compute normal vectors for triangles pab and pbc
	sxyz<T1> u = cross(b, c);
	sxyz<T1> v = cross(c, a);
	// Make sure they are both pointing in the same direction
	if (dot(u, v) < static_cast<T1>(0.0)) return false;
	// Compute normal vector for triangle pca
	sxyz<T1> w = cross(a, b);
	// Make sure it points in the same direction as the first two
	if (dot(u, w) < static_cast<T1>(0.0)) return false;
	// Otherwise P must be in (or on) the triangle
	return true;
}


template<typename T1, typename T2, typename NODE, typename S>
void Grid2Dui<T1,T2,NODE,S>::saveTT(const std::string &fname, const int all,
                                    const size_t nt, const bool vtkFormat) const {
	
	if (vtkFormat) {
#ifdef VTK
		vtkSmartPointer<vtkUnstructuredGrid> ugrid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
		
		vtkSmartPointer<vtkPoints> newPts =
		vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkDoubleArray> newScalars =
		vtkSmartPointer<vtkDoubleArray>::New();
		
		newScalars->SetName("Travel time");
		
		double xyz[3];
		T2 nMax = nPrimary;  // only primary are saved
		for (size_t n=0; n<nMax; ++n) {
			xyz[0] = nodes[n].getX();
			xyz[1] = nodes[n].getY();
			xyz[2] = nodes[n].getZ();
			newPts->InsertPoint(n, xyz);
			newScalars->InsertValue(n, nodes[n].getTT(nt) );
		}
		
		ugrid->SetPoints(newPts);
		ugrid->GetPointData()->SetScalars(newScalars);
		
		vtkSmartPointer<vtkTriangle> tri =
		vtkSmartPointer<vtkTriangle>::New();
		for (size_t n=0; n<triangles.size(); ++n) {
			tri->GetPointIds()->SetId(0, triangles[n].i[0] );
			tri->GetPointIds()->SetId(1, triangles[n].i[1] );
			tri->GetPointIds()->SetId(2, triangles[n].i[2] );
            
			ugrid->InsertNextCell( tri->GetCellType(), tri->GetPointIds() );
		}
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
		vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		
		writer->SetFileName( fname.c_str() );
		writer->SetInput( ugrid );
		writer->SetDataModeToBinary();
		writer->Update();
#else
		std::cerr << "VTK not included during compilation.\nNothing saved.\n";
#endif
	} else {
        std::ofstream fout(fname.c_str());
        fout.precision(12);
        T2 nMax = nPrimary;
        if ( all == 1 ) {
            nMax = static_cast<T2>(nodes.size());
        }
        for ( T2 n=0; n<nMax; ++n ) {
            fout << nodes[n].getX() << '\t'
            << nodes[n].getZ() << '\t'
            << nodes[n].getTT(nt) << '\n';
        }
        fout.close();
	}
}

template<typename T1, typename T2, typename NODE, typename S>
void Grid2Dui<T1,T2,NODE,S>::saveModelVTU(const std::string &fname,
                                          const bool saveSlowness) const {
    
    vtkSmartPointer<vtkUnstructuredGrid> ugrid =
    vtkSmartPointer<vtkUnstructuredGrid>::New();
    
    vtkSmartPointer<vtkPoints> newPts =
    vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkDoubleArray> newScalars =
    vtkSmartPointer<vtkDoubleArray>::New();
    
    newScalars->SetName("Slowness");

    double xyz[3];
    T2 nMax = nPrimary;  // only primary are saved
    for (size_t n=0; n<nMax; ++n) {
        xyz[0] = nodes[n].getX();
        xyz[1] = nodes[n].getY();
        xyz[2] = nodes[n].getZ();
        newPts->InsertPoint(n, xyz);
        newScalars->InsertValue(n, nodes[n].getNodeSlowness() );
    }
    
    ugrid->SetPoints(newPts);
    ugrid->GetPointData()->SetScalars(newScalars);
        
    vtkSmartPointer<vtkTriangle> tri =
    vtkSmartPointer<vtkTriangle>::New();
    for (size_t n=0; n<triangles.size(); ++n) {
        tri->GetPointIds()->SetId(0, triangles[n].i[0] );
        tri->GetPointIds()->SetId(1, triangles[n].i[1] );
        tri->GetPointIds()->SetId(2, triangles[n].i[2] );
        
        ugrid->InsertNextCell( tri->GetCellType(), tri->GetPointIds() );
    }
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    
    writer->SetFileName( fname.c_str() );
    writer->SetInput( ugrid );
    writer->SetDataModeToBinary();
    writer->Update();
}


#endif
