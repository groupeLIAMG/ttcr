//
//  Interpolator.h
//  ttcr.v2
//
//  Created by Bernard Giroux on 2012-11-23.
//  Copyright (c) 2012 Bernard Giroux. All rights reserved.
//  Copyright (c) 2018 Bernard Giroux, Maher Nasr. All rights reserved.
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

#ifndef ttcr_Interpolator_h
#define ttcr_Interpolator_h

#include <array>
#include <cmath>

namespace ttcr {

    template<class T> class Interpolator
    {
    public:
        inline static T linear(const T x[], const T s[]) {

            // evaluate s @ x[0]
            // with
            // s[0] @ x[1]
            // s[1] @ x[2]

            return (s[0]*(x[2]-x[0]) + s[1]*(x[0]-x[1]))/(x[2]-x[1]);
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

        template<typename NODE>
        static T linear(const sxyz<T>& node,
                        const NODE& node1,
                        const NODE& node2) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> IA = {node1.getX()-node.x, node1.getY()-node.y, node1.getZ()-node.z};
            sxyz<T> IB = {node2.getX()-node.x, node2.getY()-node.y, node2.getZ()-node.z};

            T nAB = norm(AB);
            T w1 = norm(IA)/nAB;
            T w2 = norm(IB)/nAB;

            return w2*node1.getNodeSlowness() + w1*node2.getNodeSlowness();
        }

        template<typename NODE>
        static T linearVel(const sxyz<T>& node,
                           const NODE& node1,
                           const NODE& node2) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> IA = {node1.getX()-node.x, node1.getY()-node.y, node1.getZ()-node.z};
            sxyz<T> IB = {node2.getX()-node.x, node2.getY()-node.y, node2.getZ()-node.z};

            T nAB = norm(AB);
            T w1 = norm(IA)/nAB;
            T w2 = norm(IB)/nAB;

            return 1.0/(w2/node1.getNodeSlowness() + w1/node2.getNodeSlowness());
        }

        template<typename NODE>
        static T barycentricTriangle(const sxz<T>& node,
                                     const NODE& node1,
                                     const NODE& node2,
                                     const NODE& node3) {
            T den = (node2.getZ()-node3.getZ())*(node1.getX()-node3.getX()) +
            (node3.getX()-node2.getX())*(node1.getZ()-node3.getZ());
            T w1 = ((node2.getZ()-node3.getZ())*(node.x-node3.getX()) +
                    (node3.getX()-node2.getX())*(node.z-node3.getZ())) / den;
            T w2 = ((node3.getZ()-node1.getZ())*(node.x-node3.getX()) +
                    (node1.getX()-node3.getX())*(node.z-node3.getZ())) / den;
            T w3 = 1. - w1 - w2;

            return w1*node1.getNodeSlowness() + w2*node2.getNodeSlowness() + w3*node3.getNodeSlowness();
        }

        template<typename NODE>
        static T barycentricTriangle(const sxz<T>& node,
                                     const std::vector<NODE*> &inodes) {
            T den = (inodes[1]->getZ()-inodes[2]->getZ())*(inodes[0]->getX()-inodes[2]->getX()) +
            (inodes[2]->getX()-inodes[1]->getX())*(inodes[0]->getZ()-inodes[2]->getZ());
            T w1 = ((inodes[1]->getZ()-inodes[2]->getZ())*(node.x-inodes[2]->getX()) +
                    (inodes[2]->getX()-inodes[1]->getX())*(node.z-inodes[2]->getZ())) / den;
            T w2 = ((inodes[2]->getZ()-inodes[0]->getZ())*(node.x-inodes[2]->getX()) +
                    (inodes[0]->getX()-inodes[2]->getX())*(node.z-inodes[2]->getZ())) / den;
            T w3 = 1. - w1 - w2;

            return w1*inodes[0]->getNodeSlowness() + w2*inodes[1]->getNodeSlowness() + w3*inodes[2]->getNodeSlowness();
        }

        template<typename NODE>
        static T barycentricTriangleVel(const sxz<T>& node,
                                        const NODE& node1,
                                        const NODE& node2,
                                        const NODE& node3) {
            T den = (node2.getZ()-node3.getZ())*(node1.getX()-node3.getX()) +
            (node3.getX()-node2.getX())*(node1.getZ()-node3.getZ());
            T w1 = ((node2.getZ()-node3.getZ())*(node.x-node3.getX()) +
                    (node3.getX()-node2.getX())*(node.z-node3.getZ())) / den;
            T w2 = ((node3.getZ()-node1.getZ())*(node.x-node3.getX()) +
                    (node1.getX()-node3.getX())*(node.z-node3.getZ())) / den;
            T w3 = 1. - w1 - w2;

            return (1.0/(w1/node1.getNodeSlowness() + w2/node2.getNodeSlowness() + w3/node3.getNodeSlowness()));
        }

        template<typename NODE>
        static T bilinearTriangleVel(const sxyz<T>& node,
                                     const NODE& node1,
                                     const NODE& node2,
                                     const NODE& node3) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};

            T S = norm(cross(AB,AC));

            sxyz<T> IB = {node2.getX()-node.x, node2.getY()-node.y, node2.getZ()-node.z};
            sxyz<T> IC = {node3.getX()-node.x, node3.getY()-node.y, node3.getZ()-node.z};

            T w1 = norm(cross(IB,IC))/S;
            T w2 = norm(cross(AC,IC))/S;
            T w3 = norm(cross(IB,AB))/S;

            return (1.0/(w1/node1.getNodeSlowness() + w2/node2.getNodeSlowness() + w3/node3.getNodeSlowness()));
        }

        template<typename NODE>
        static T bilinearTriangle(const sxyz<T>& node,
                                  const NODE& node1,
                                  const NODE& node2,
                                  const NODE& node3) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};

            T S = norm(cross(AB,AC));

            sxyz<T> IB = {node2.getX()-node.x, node2.getY()-node.y, node2.getZ()-node.z};
            sxyz<T> IC = {node3.getX()-node.x, node3.getY()-node.y, node3.getZ()-node.z};

            T w1 = norm(cross(IB,IC))/S;
            T w2 = norm(cross(AC,IC))/S;
            T w3 = norm(cross(IB,AB))/S;

            return w1*node1.getNodeSlowness() + w2*node2.getNodeSlowness() + w3*node3.getNodeSlowness();
        }

        template<typename NODE>
        static T bilinearTriangleVel(const NODE& node,
                                     const NODE& node1,
                                     const NODE& node2,
                                     const NODE& node3) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};

            T S = norm(cross(AB,AC));

            sxyz<T> IB = {node2.getX()-node.getX(), node2.getY()-node.getY(), node2.getZ()-node.getZ()};
            sxyz<T> IC = {node3.getX()-node.getX(), node3.getY()-node.getY(), node3.getZ()-node.getZ()};

            T w1 = norm(cross(IB,IC))/S;
            T w2 = norm(cross(AC,IC))/S;
            T w3 = norm(cross(IB,AB))/S;

            return (1.0/(w1/node1.getNodeSlowness() + w2/node2.getNodeSlowness() + w3/node3.getNodeSlowness()));
        }

        template<typename NODE1, typename NODE2>
        static T bilinearTriangleVel(const NODE1& node, const std::vector<NODE2*> &inodes) {
            sxyz<T> AB = {inodes[1]->getX()-inodes[0]->getX(), inodes[1]->getY()-inodes[0]->getY(), inodes[1]->getZ()-inodes[0]->getZ()};
            sxyz<T> AC = {inodes[2]->getX()-inodes[0]->getX(), inodes[2]->getY()-inodes[0]->getY(), inodes[2]->getZ()-inodes[0]->getZ()};

            T S = norm(cross(AB,AC));

            sxyz<T> IB = {inodes[1]->getX()-node.getX(), inodes[1]->getY()-node.getY(), inodes[1]->getZ()-node.getZ()};
            sxyz<T> IC = {inodes[2]->getX()-node.getX(), inodes[2]->getY()-node.getY(), inodes[2]->getZ()-node.getZ()};

            T w1 = norm(cross(IB,IC))/S;
            T w2 = norm(cross(AC,IC))/S;
            T w3 = norm(cross(IB,AB))/S;

            return (1.0/(w1/inodes[0]->getNodeSlowness() + w2/inodes[1]->getNodeSlowness() + w3/inodes[2]->getNodeSlowness()));
        }

        template<typename NODE1, typename NODE2>
        static T bilinearTriangle(const NODE1& node, const std::vector<NODE2*> &inodes) {
            sxyz<T> AB = {inodes[1]->getX()-inodes[0]->getX(), inodes[1]->getY()-inodes[0]->getY(), inodes[1]->getZ()-inodes[0]->getZ()};
            sxyz<T> AC = {inodes[2]->getX()-inodes[0]->getX(), inodes[2]->getY()-inodes[0]->getY(), inodes[2]->getZ()-inodes[0]->getZ()};

            T S = norm(cross(AB,AC));

            sxyz<T> IB = {inodes[1]->getX()-node.getX(), inodes[1]->getY()-node.getY(), inodes[1]->getZ()-node.getZ()};
            sxyz<T> IC = {inodes[2]->getX()-node.getX(), inodes[2]->getY()-node.getY(), inodes[2]->getZ()-node.getZ()};

            T w1 = norm(cross(IB,IC))/S;
            T w2 = norm(cross(AC,IC))/S;
            T w3 = norm(cross(IB,AB))/S;

            return w1*inodes[0]->getNodeSlowness() + w2*inodes[1]->getNodeSlowness() + w3*inodes[2]->getNodeSlowness();
        }

        template<typename NODE>
        static void bilinearTriangleWeight(const sxyz<T>& node, const NODE* node1,
                                           const NODE* node2, const NODE* node3,
                                           std::array<T,3>& w) {
            sxyz<T> AB = {node2->getX()-node1->getX(), node2->getY()-node1->getY(), node2->getZ()-node1->getZ()};
            sxyz<T> AC = {node3->getX()-node1->getX(), node3->getY()-node1->getY(), node3->getZ()-node1->getZ()};

            T S = norm(cross(AB,AC));

            sxyz<T> IB = {node2->getX()-node.x, node2->getY()-node.y, node2->getZ()-node.z};
            sxyz<T> IC = {node3->getX()-node.x, node3->getY()-node.y, node3->getZ()-node.z};

            w[0] = norm(cross(IB,IC))/S;
            w[1] = norm(cross(AC,IC))/S;
            w[2] = norm(cross(IB,AB))/S;
        }

        template<typename NODE>
        static T bilinearTime(const sxyz<T> &node, const NODE &node1,
                              const NODE &node2, NODE &node3, const size_t nt) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};

            T S = norm(cross(AB,AC));

            sxyz<T> IB = {node2.getX()-node.x, node2.getY()-node.y, node2.getZ()-node.z};
            sxyz<T> IC = {node3.getX()-node.x, node3.getY()-node.y, node3.getZ()-node.z};

            T w1 = norm(cross(IB,IC))/S;
            T w2 = norm(cross(AC,IC))/S;
            T w3 = norm(cross(IB,AB))/S;

            return w1*node1.getTT(nt) + w2*node2.getTT(nt) + w3*node3.getTT(nt);
        }


        template<typename NODE>
        static void trilinearTriangleWeight(const sxyz<T>& node, const NODE* node1,
                                            const NODE* node2, const NODE* node3,
                                            const NODE* node4, std::array<T,4>& w) {
            sxyz<T> AB = {node2->getX()-node1->getX(), node2->getY()-node1->getY(), node2->getZ()-node1->getZ()};
            sxyz<T> AC = {node3->getX()-node1->getX(), node3->getY()-node1->getY(), node3->getZ()-node1->getZ()};
            sxyz<T> AD = {node4->getX()-node1->getX(), node4->getY()-node1->getY(), node4->getZ()-node1->getZ()};

            T V = std::abs(det(AB, AC, AD));

            sxyz<T> IA = {node1->getX()-node.x, node1->getY()-node.y, node1->getZ()-node.z};
            sxyz<T> IB = {node2->getX()-node.x, node2->getY()-node.y, node2->getZ()-node.z};
            sxyz<T> IC = {node3->getX()-node.x, node3->getY()-node.y, node3->getZ()-node.z};
            sxyz<T> ID = {node4->getX()-node.x, node4->getY()-node.y, node4->getZ()-node.z};

            w[0] = std::abs(det(IB, IC, ID))/V;
            w[1] = std::abs(det(IC, IA, ID))/V;
            w[2] = std::abs(det(IB, IA, ID))/V;
            w[3] = std::abs(det(IB, IA, IC))/V;
        }

        template<typename NODE>
        static T trilinearTriangle(const sxyz<T>& node, const NODE& node1,
                                   const NODE& node2, const NODE& node3,
                                   const NODE& node4) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};
            sxyz<T> AD = {node4.getX()-node1.getX(), node4.getY()-node1.getY(), node4.getZ()-node1.getZ()};

            T V = std::abs(det(AB, AC, AD));

            sxyz<T> IA = {node1.getX()-node.x, node1.getY()-node.y, node1.getZ()-node.z};
            sxyz<T> IB = {node2.getX()-node.x, node2.getY()-node.y, node2.getZ()-node.z};
            sxyz<T> IC = {node3.getX()-node.x, node3.getY()-node.y, node3.getZ()-node.z};
            sxyz<T> ID = {node4.getX()-node.x, node4.getY()-node.y, node4.getZ()-node.z};

            T w1 = std::abs(det(IB, IC, ID))/V;
            T w2 = std::abs(det(IC, IA, ID))/V;
            T w3 = std::abs(det(IB, IA, ID))/V;
            T w4 = std::abs(det(IB, IA, IC))/V;

            return (w1*node1.getNodeSlowness() + w2*node2.getNodeSlowness() + w3*node3.getNodeSlowness() + w4*node4.getNodeSlowness());
        }

        template<typename NODE>
        static T trilinearTriangleVel(const sxyz<T>& node, const NODE& node1,
                                      const NODE& node2, const NODE& node3,
                                      const NODE& node4) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};
            sxyz<T> AD = {node4.getX()-node1.getX(), node4.getY()-node1.getY(), node4.getZ()-node1.getZ()};

            T V = std::abs(det(AB, AC, AD));

            sxyz<T> IA = {node1.getX()-node.x, node1.getY()-node.y, node1.getZ()-node.z};
            sxyz<T> IB = {node2.getX()-node.x, node2.getY()-node.y, node2.getZ()-node.z};
            sxyz<T> IC = {node3.getX()-node.x, node3.getY()-node.y, node3.getZ()-node.z};
            sxyz<T> ID = {node4.getX()-node.x, node4.getY()-node.y, node4.getZ()-node.z};

            T w1 = std::abs(det(IB, IC, ID))/V;
            T w2 = std::abs(det(IC, IA, ID))/V;
            T w3 = std::abs(det(IB, IA, ID))/V;
            T w4 = std::abs(det(IB, IA, IC))/V;

            return 1.0/(w1/node1.getNodeSlowness() + w2/node2.getNodeSlowness() + w3/node3.getNodeSlowness() + w4/node4.getNodeSlowness());
        }

        template<typename NODE>
        static T trilinearTriangle(const NODE& node, const NODE& node1,
                                   const NODE& node2, const NODE& node3,
                                   const NODE& node4) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};
            sxyz<T> AD = {node4.getX()-node1.getX(), node4.getY()-node1.getY(), node4.getZ()-node1.getZ()};

            T V = std::abs(det(AB, AC, AD));

            sxyz<T> IA = {node1.getX()-node.getX(), node1.getY()-node.getY(), node1.getZ()-node.getZ()};
            sxyz<T> IB = {node2.getX()-node.getX(), node2.getY()-node.getY(), node2.getZ()-node.getZ()};
            sxyz<T> IC = {node3.getX()-node.getX(), node3.getY()-node.getY(), node3.getZ()-node.getZ()};
            sxyz<T> ID = {node4.getX()-node.getX(), node4.getY()-node.getY(), node4.getZ()-node.getZ()};

            T w1 = std::abs(det(IB, IC, ID))/V;
            T w2 = std::abs(det(IC, IA, ID))/V;
            T w3 = std::abs(det(IB, IA, ID))/V;
            T w4 = std::abs(det(IB, IA, IC))/V;

            return (w1*node1.getNodeSlowness() + w2*node2.getNodeSlowness() + w3*node3.getNodeSlowness() + w4*node4.getNodeSlowness());
        }

        template<typename NODE>
        static T trilinearTriangle(const sxyz<T>& node,
                                   const std::vector<NODE*> nodes) {

            assert(nodes.size()==4);
            sxyz<T> AB = {nodes[1]->getX()-nodes[0]->getX(), nodes[1]->getY()-nodes[0]->getY(), nodes[1]->getZ()-nodes[0]->getZ()};
            sxyz<T> AC = {nodes[2]->getX()-nodes[0]->getX(), nodes[2]->getY()-nodes[0]->getY(), nodes[2]->getZ()-nodes[0]->getZ()};
            sxyz<T> AD = {nodes[3]->getX()-nodes[0]->getX(), nodes[3]->getY()-nodes[0]->getY(), nodes[3]->getZ()-nodes[0]->getZ()};

            T V = std::abs(det(AB, AC, AD));

            sxyz<T> IA = {nodes[0]->getX()-node.x, nodes[0]->getY()-node.y, nodes[0]->getZ()-node.z};
            sxyz<T> IB = {nodes[1]->getX()-node.x, nodes[1]->getY()-node.y, nodes[1]->getZ()-node.z};
            sxyz<T> IC = {nodes[2]->getX()-node.x, nodes[2]->getY()-node.y, nodes[2]->getZ()-node.z};
            sxyz<T> ID = {nodes[3]->getX()-node.x, nodes[3]->getY()-node.y, nodes[3]->getZ()-node.z};

            T w1 = std::abs(det(IB, IC, ID))/V;
            T w2 = std::abs(det(IC, IA, ID))/V;
            T w3 = std::abs(det(IB, IA, ID))/V;
            T w4 = std::abs(det(IB, IA, IC))/V;

            return (w1*nodes[0]->getNodeSlowness() + w2*nodes[1]->getNodeSlowness() + w3*nodes[2]->getNodeSlowness() + w4*nodes[3]->getNodeSlowness());
        }

        template<typename NODE>
        static T trilinearTriangleVel(const sxyz<T>& node,
                                      const std::vector<NODE*> nodes) {

            assert(nodes.size()==4);
            sxyz<T> AB = {nodes[1]->getX()-nodes[0]->getX(), nodes[1]->getY()-nodes[0]->getY(), nodes[1]->getZ()-nodes[0]->getZ()};
            sxyz<T> AC = {nodes[2]->getX()-nodes[0]->getX(), nodes[2]->getY()-nodes[0]->getY(), nodes[2]->getZ()-nodes[0]->getZ()};
            sxyz<T> AD = {nodes[3]->getX()-nodes[0]->getX(), nodes[3]->getY()-nodes[0]->getY(), nodes[3]->getZ()-nodes[0]->getZ()};

            T V = std::abs(det(AB, AC, AD));

            sxyz<T> IA = {nodes[0]->getX()-node.x, nodes[0]->getY()-node.y, nodes[0]->getZ()-node.z};
            sxyz<T> IB = {nodes[1]->getX()-node.x, nodes[1]->getY()-node.y, nodes[1]->getZ()-node.z};
            sxyz<T> IC = {nodes[2]->getX()-node.x, nodes[2]->getY()-node.y, nodes[2]->getZ()-node.z};
            sxyz<T> ID = {nodes[3]->getX()-node.x, nodes[3]->getY()-node.y, nodes[3]->getZ()-node.z};

            T w1 = std::abs(det(IB, IC, ID))/V;
            T w2 = std::abs(det(IC, IA, ID))/V;
            T w3 = std::abs(det(IB, IA, ID))/V;
            T w4 = std::abs(det(IB, IA, IC))/V;

            return 1.0/(w1/nodes[0]->getNodeSlowness() + w2/nodes[1]->getNodeSlowness() + w3/nodes[2]->getNodeSlowness() + w4/nodes[3]->getNodeSlowness());
        }

        template<typename NODE>
        static T trilinearTriangleVel(const NODE& node, const NODE& node1,
                                      const NODE& node2, const NODE& node3,
                                      const NODE& node4) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};
            sxyz<T> AD = {node4.getX()-node1.getX(), node4.getY()-node1.getY(), node4.getZ()-node1.getZ()};

            T V = std::abs(det(AB, AC, AD));

            sxyz<T> IA = {node1.getX()-node.getX(), node1.getY()-node.getY(), node1.getZ()-node.getZ()};
            sxyz<T> IB = {node2.getX()-node.getX(), node2.getY()-node.getY(), node2.getZ()-node.getZ()};
            sxyz<T> IC = {node3.getX()-node.getX(), node3.getY()-node.getY(), node3.getZ()-node.getZ()};
            sxyz<T> ID = {node4.getX()-node.getX(), node4.getY()-node.getY(), node4.getZ()-node.getZ()};

            T w1 = std::abs(det(IB, IC, ID))/V;
            T w2 = std::abs(det(IC, IA, ID))/V;
            T w3 = std::abs(det(IB, IA, ID))/V;
            T w4 = std::abs(det(IB, IA, IC))/V;

            return 1.0/(w1/node1.getNodeSlowness() + w2/node2.getNodeSlowness() + w3/node3.getNodeSlowness() + w4/node4.getNodeSlowness());
        }

        template<typename NODE>
        static T trilinearTime(const sxyz<T> &node, const NODE &node1,
                               const NODE &node2, const NODE &node3,
                               const NODE &node4, const size_t nt) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};
            sxyz<T> AD = {node4.getX()-node1.getX(), node4.getY()-node1.getY(), node4.getZ()-node1.getZ()};

            T V = std::abs(det(AB,AC,AD));

            sxyz<T> IA = {node1.getX()-node.x, node1.getY()-node.y, node1.getZ()-node.z};
            sxyz<T> IB = {node2.getX()-node.x, node2.getY()-node.y, node2.getZ()-node.z};
            sxyz<T> IC = {node3.getX()-node.x, node3.getY()-node.y, node3.getZ()-node.z};
            sxyz<T> ID = {node4.getX()-node.x, node4.getY()-node.y, node4.getZ()-node.z};

            T w1 = std::abs(det(IB,IC,ID))/V;
            T w2 = std::abs(det(IC,IA,ID))/V;
            T w3 = std::abs(det(IB,IA,ID))/V;
            T w4 = std::abs(det(IB,IA,IC))/V;

            return w1*node1.getTT(nt) + w2*node2.getTT(nt) + w3*node3.getTT(nt) + w4*node4.getTT(nt);
        }
    };

}

#endif
