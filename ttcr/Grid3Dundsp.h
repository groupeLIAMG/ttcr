//
//  Grid3Dundsp.h
//  ttcr
//
//  Created by Bernard Giroux on 18-10-09.
//  Copyright © 2018 Bernard Giroux. All rights reserved.
//

#ifndef Grid3Dundsp_h
#define Grid3Dundsp_h

#include "Grid3Dun.h"

namespace ttcr {
    
    template<typename T1, typename T2>
    class Grid3Dundsp : public Grid3Dun<T1,T2,Node3Dn<T1,T2>> {
    public:
        Grid3Dundsp(const std::vector<sxyz<T1>>& no,
                   const std::vector<tetrahedronElem<T2>>& tet,
                   const int ns, const int nd, const size_t nt=1,
                    const int verbose=0) :
        Grid3Dun<T1,T2,Node3Dn<T1,T2>>(no, tet, nt), nsecondary(ns), ndynamic(nd)
        {
            this->buildGridNodes(no, ns, nt, verbose);
            this->buildGridNeighbors();
        }
        
        ~Grid3Dundsp() {
        }
        
    private:
        T2 ndynamic;
        
    };
}
#endif /* Grid3Dundsp_h */
