//
//  Interpolator.h
//  ttcr.v2
//
//  Created by Bernard Giroux on 2012-11-23.
//  Copyright (c) 2012 INRS-ETE. All rights reserved.
//

#ifndef ttcr_v2_Interpolator_h
#define ttcr_v2_Interpolator_h

#include "Node3Di.h"

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
	
	template<typename T2>
    static T inverseDistance(const Node3Di<T,T2> &node,
                             const std::vector<Node3Di<T,T2>*> &inodes) {
        
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
    
	template<typename T2>
    static T inverseDistance(const sxyz<T>& node,
                             const std::vector<Node3Di<T,T2>*> &inodes) {
        
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
