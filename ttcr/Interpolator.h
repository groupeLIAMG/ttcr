//
//  Interpolator.h
//  ttcr.v2
//
//  Created by Bernard Giroux on 2012-11-23.
//  Copyright (c) 2012 INRS-ETE. All rights reserved.
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

#ifndef ttcr_v2_Interpolator_h
#define ttcr_v2_Interpolator_h

#include "Node3Disp.h"

template<class T> class Interpolator
{
public:
	inline static T linear(const T x[], const T y[], const T s[]) {
        
        // evaluate s @ x[0]
        // with
        // s[0] @ x[1]
        // s[1] @ x[2]
        
        return (s[0]*(x[2]-x[0]) - s[1]*(x[1]-x[0]))/(x[1]-x[2]);
    }
	
	inline static T bilinear(const T x[], const T y[], const T s[]) {
        
        // evaluate s @ (x[0], y[0])
        // with
        // s[0] @ (x[1], y[1])
        // s[1] @ (x[1], y[2])
        // s[2] @ (x[2], y[1])
        // s[3] @ (x[2], y[2])
        
        return (s[0]*(x[2]-x[0])*(y[2]-y[0]) +
                s[1]*(x[2]-x[0])*(y[0]-y[1]) +
                s[2]*(x[0]-x[1])*(y[2]-y[0]) +
                s[3]*(x[0]-x[1])*(y[0]-y[1]))
        /((x[2]-x[1])*(y[2]-y[1]));
    }
	
	inline static T trilinear(const T x[], const T y[], const T z[], const T s[]) {
        
        // evaluate s @ (x[0], y[0], z[0])
        // with
        // s[0] @ (x[1], y[1], z[1])
        // s[1] @ (x[1], y[1], z[2])
        // s[2] @ (x[1], y[2], z[1])
        // s[3] @ (x[1], y[2], z[2])
        // s[4] @ (x[2], y[1], z[1])
        // s[5] @ (x[2], y[1], z[2])
        // s[6] @ (x[2], y[2], z[1])
        // s[7] @ (x[2], y[2], z[2])
        
        return (s[0]*(x[2]-x[0])*(y[2]-y[0])*(z[2]-z[0]) +
                s[1]*(x[2]-x[0])*(y[2]-y[0])*(z[0]-z[1]) +
                s[2]*(x[2]-x[0])*(y[0]-y[1])*(z[2]-z[0]) +
                s[3]*(x[2]-x[0])*(y[0]-y[1])*(z[0]-z[1]) +
                s[4]*(x[0]-x[1])*(y[2]-y[0])*(z[2]-z[0]) +
                s[5]*(x[0]-x[1])*(y[2]-y[0])*(z[0]-z[1]) +
                s[6]*(x[0]-x[1])*(y[0]-y[1])*(z[2]-z[0]) +
                s[7]*(x[0]-x[1])*(y[0]-y[1])*(z[0]-z[1]))
        /((x[2]-x[1])*(y[2]-y[1])*(z[2]-z[1]));
    }
	
//	inline static T lagrange2d(const T x[], const T y[], const T s[]) {
//        
//        // evaluate s @ (x[0], y[0])
//        // with
//        // s[0] @ (x[1], y[1])
//        // s[1] @ (x[1], y[2])
//        // s[2] @ (x[2], y[1])
//        // s[3] @ (x[2], y[2])
//        
//        T w[4] = {
//            ((x[0]-x[2])*(y[0]-y[2]))/((x[1]-x[2])*(y[1]-y[2])),
//            ((x[0]-x[2])*(y[0]-y[1]))/((x[1]-x[2])*(y[2]-y[1])),
//            ((x[0]-x[1])*(y[0]-y[2]))/((x[2]-x[1])*(y[1]-y[2])),
//            ((x[0]-x[1])*(y[0]-y[1]))/((x[2]-x[1])*(y[2]-y[1]))
//        };
//        
//        return s[0]*w[0] + s[1]*w[1] + s[2]*w[2] + s[3]*w[3];
//    }
//	
//	inline static T lagrange3d(const T x[], const T y[], const T z[], const T s[]) {
//        
//        // evaluate s @ (x[0], y[0], z[0])
//        // with
//        // s[0] @ (x[1], y[1], z[1])
//        // s[1] @ (x[1], y[1], z[2])
//        // s[2] @ (x[1], y[2], z[1])
//        // s[3] @ (x[1], y[2], z[2])
//        // s[4] @ (x[2], y[1], z[1])
//        // s[5] @ (x[2], y[1], z[2])
//        // s[6] @ (x[2], y[2], z[1])
//        // s[7] @ (x[2], y[2], z[2])
//        
//        T w[8] = {
//            ((x[0]-x[2]) * (y[0]-y[2]) * (z[0]-z[2])) / ((x[1]-x[2]) * (y[1]-y[2]) * (z[1]-z[2])),
//            ((x[0]-x[2]) * (y[0]-y[2]) * (z[0]-z[1])) / ((x[1]-x[2]) * (y[1]-y[2]) * (z[2]-z[1])),
//            ((x[0]-x[2]) * (y[0]-y[1]) * (z[0]-z[2])) / ((x[1]-x[2]) * (y[2]-y[1]) * (z[1]-z[2])),
//            ((x[0]-x[2]) * (y[0]-y[1]) * (z[0]-z[1])) / ((x[1]-x[2]) * (y[2]-y[1]) * (z[2]-z[1])),
//            ((x[0]-x[1]) * (y[0]-y[2]) * (z[0]-z[2])) / ((x[2]-x[1]) * (y[1]-y[2]) * (z[1]-z[2])),
//            ((x[0]-x[1]) * (y[0]-y[2]) * (z[0]-z[1])) / ((x[2]-x[1]) * (y[1]-y[2]) * (z[2]-z[1])),
//            ((x[0]-x[1]) * (y[0]-y[1]) * (z[0]-z[2])) / ((x[2]-x[1]) * (y[2]-y[1]) * (z[1]-z[2])),
//            ((x[0]-x[1]) * (y[0]-y[1]) * (z[0]-z[1])) / ((x[2]-x[1]) * (y[2]-y[1]) * (z[2]-z[1]))
//        };
//        
//        return s[0]*w[0] + s[1]*w[1] + s[2]*w[2] + s[3]*w[3] + s[4]*w[4] + s[5]*w[5] + s[6]*w[6] + s[7]*w[7];
//    }
	
	template<typename NODE>
    static T inverseDistance(const NODE &node,
                             const std::vector<NODE*> &inodes) {
        
        T num=0.;
        T den=0.;
        T w;
        
        for ( size_t n=0; n<inodes.size(); ++n ) {
            w = 1./inodes[n]->getDistance( node );
            num += w*inodes[n]->getNodeSlowness();
            den += w;
        }
        
        return num/den;
    }
    
	template<typename NODE, typename S>
    static T inverseDistance(const S& node,
                             const std::vector<NODE*> &inodes) {
        
        T num=0.;
        T den=0.;
        T w;
        
        for ( size_t n=0; n<inodes.size(); ++n ) {
            w = 1./inodes[n]->getDistance( node );
            num += w*inodes[n]->getNodeSlowness();
            den += w;
        }
        
        return num/den;
    }

};

#endif
