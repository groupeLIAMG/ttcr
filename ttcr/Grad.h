//
//  Grad.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-03-04.
//  Copyright (c) 2014 Bernard Giroux. All rights reserved.
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

#ifndef ttcr_Grad_h
#define ttcr_Grad_h

#include <cassert>

#include <array>
#include <cmath>
#include <set>
#include <vector>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdocumentation"
#include <Eigen/Dense>
#pragma clang diagnostic pop

#include "ttcr_t.h"
#include "Interpolator.h"
#include "Node.h"

namespace ttcr {

    /**
     * Compute traveltime gradient for a triangle, with first-order least-squares
     *
     *  @tparam T underlying type of node objects
     */
    template <typename T>
    class Grad2D_ls_fo {
    public:
        Grad2D_ls_fo() = default;

        /**
         * Compute gradient given triangle nodes
         *
         * @param n0 first node making the triangle
         * @param n1 second node making the triangle
         * @param n2 third node making the triangle
         * @param nt thread number
         * @returns value of the travetime gradient
         */
        sxz<T> compute(const Node<T> &n0,
                       const Node<T> &n1,
                       const Node<T> &n2,
                       const size_t nt);

    private:
        sxz<T> g;
        Eigen::Matrix<T, 3, 2> A;
        Eigen::Matrix<T, 2, 1> x;
        Eigen::Matrix<T, 3, 1> b;
    };

    template <typename T>
    sxz<T> Grad2D_ls_fo<T>::compute(const Node<T> &n0,
                                    const Node<T> &n1,
                                    const Node<T> &n2,
                                    const size_t nt) {

        // find centroid of triangle
        sxz<T> cent = {(n0.getX()+n1.getX()+n2.getX())/static_cast<T>(3.),
            (n0.getZ()+n1.getZ()+n2.getZ())/static_cast<T>(3.)};

        // time at centroid is inverse distance weeighted
        T w = 1./sqrt( (n0.getX()-cent.x)*(n0.getX()-cent.x) + (n0.getZ()-cent.z)*(n0.getZ()-cent.z) );
        T t = w*n0.getTT(nt);
        T den = w;
        w = 1./sqrt( (n1.getX()-cent.x)*(n1.getX()-cent.x) + (n1.getZ()-cent.z)*(n1.getZ()-cent.z) );
        t += w*n1.getTT(nt);
        den += w;
        w = 1./sqrt( (n2.getX()-cent.x)*(n2.getX()-cent.x) + (n2.getZ()-cent.z)*(n2.getZ()-cent.z) );
        t += w*n2.getTT(nt);
        den += w;

        t /= den;

        A(0,0) = n0.getX() - cent.x;
        A(0,1) = n0.getZ() - cent.z;
        b(0,0) = t - n0.getTT(nt);

        A(1,0) = n1.getX() - cent.x;
        A(1,1) = n1.getZ() - cent.z;
        b(1,0) = t - n1.getTT(nt);

        A(2,0) = n2.getX() - cent.x;
        A(2,1) = n2.getZ() - cent.z;
        b(2,0) = t - n2.getTT(nt);

        // solve Ax = b with least squares
        x = A.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b);

        g.x = x[0];
        g.z = x[1];
        return g;
    }


    /**
     * Compute traveltime gradient for a triangle, with second-order least-squares
     *
     * Nodes that are neighbours to the triangle are used
     *
     *  @tparam T underlying type of node object
     *  @tparam NODE node objects making the mesh
     */
    template <typename T, typename NODE>
    class Grad2D_ls_so {
    public:
        Grad2D_ls_so() = default;

        /**
         * Compute gradient given  nodes
         *
         * @param nodes set of nodes neighbours to the triangle
         * @param nt thread number
         * @returns value of the travetime gradient
         */
        sxz<T> compute(const std::set<NODE*> &nodes,
                       const size_t nt);

    private:
        sxz<T> g;
        Eigen::Matrix<T, Eigen::Dynamic, 5> A;
        Eigen::Matrix<T, 5, 1> x;
        Eigen::Matrix<T, Eigen::Dynamic, 1> b;
    };


    template <typename T, typename NODE>
    sxz<T> Grad2D_ls_so<T,NODE>::compute(const std::set<NODE*> &nodes,
                                         const size_t nt) {

        assert(nodes.size()>=5);

        sxz<T> cent = { 0.0, 0.0 };
        for ( auto n=nodes.cbegin(); n!=nodes.cend(); ++n ) {
            cent.x += (*n)->getX();
            cent.z += (*n)->getZ();
        }
        T den = 1./nodes.size();
        cent.x *= den;
        cent.z *= den;

        T t = 0.0;
        den = 0.0;

        for ( auto n=nodes.cbegin(); n!=nodes.cend(); ++n ) {
            T w = 1./sqrt(((*n)->getX()-cent.x)*((*n)->getX()-cent.x) +
                          ((*n)->getZ()-cent.z)*((*n)->getZ()-cent.z) );
            t += w*(*n)->getTT(nt);
            den += w;
        }
        t /= den;

        A.resize( nodes.size(), 5 );
        b.resize( nodes.size(), 1 );

        size_t i=0;
        for ( auto n=nodes.cbegin(); n!=nodes.cend(); ++n ) {
            T dx = (*n)->getX()-cent.x;
            T dz = (*n)->getZ()-cent.z;

            A(i,0) = dx;
            A(i,1) = dz;
            A(i,2) = dx*dx;
            A(i,3) = dz*dz;
            A(i,4) = dx*dz;

            b(i,0) = t - (*n)->getTT(nt);
            i++;
        }

        // solve Ax = b with least squares
        //        x = (A.transpose() * A).ldlt().solve(A.transpose() * b);
        x = A.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b);

        //	Eigen::Matrix<T, Eigen::Dynamic, 1> e = b-A*x;

        g.x = x[0];
        g.z = x[1];

        return g;
    }

    /**
     * Interface to 3D traveltime gradient computation classes
     *
     * @tparam T underlying type of node object
     * @tparam NODE node objects making the mesh
     */
    template <typename T, typename NODE>
    class Grad3D {
    public:

        virtual ~Grad3D() = default;

        /**
         * Compute gradient at a given point
         *
         * @param pt point where to compute gradient
         * @param t traveltime at point pt
         * @param nodes nodes surrounding point pt
         * @param nt thread number
         * @returns value of the travetime gradient
         */
        virtual sxyz<T> compute(const sxyz<T> &pt,
                                const T t,
                                const std::set<NODE*> &nodes,
                                const size_t nt) = 0;
    };

    /**
     * Compute traveltime gradient for a tetrahedron, with first-order least-squares
     *
     * @tparam T underlying type of node objects
     * @tparam NODE node objects making the mesh
     */
    template <typename T, typename NODE>
    class Grad3D_ls_fo : public Grad3D<T, NODE> {
    public:
        Grad3D_ls_fo() = default;
        ~Grad3D_ls_fo() = default;

        sxyz<T> compute(const sxyz<T> &pt,
                        const T t,
                        const std::set<NODE*> &nodes,
                        const size_t nt) override;

    private:
        sxyz<T> g;

        Eigen::Matrix<T, Eigen::Dynamic, 3> A;
        Eigen::Matrix<T, 3, 1> x;
        Eigen::Matrix<T, Eigen::Dynamic, 1> b;
    };


    template <typename T, typename NODE>
    sxyz<T> Grad3D_ls_fo<T,NODE>::compute(const sxyz<T> &pt,
                                          const T t,
                                          const std::set<NODE*> &nodes,
                                          const size_t nt) {

        assert(nodes.size()>=4);

        A.resize( nodes.size(), 3 );
        b.resize( nodes.size(), 1 );

        size_t i=0;
        for ( auto n=nodes.cbegin(); n!=nodes.cend(); ++n ) {
            A(i,0) = (*n)->getX()-pt.x;
            A(i,1) = (*n)->getY()-pt.y;
            A(i,2) = (*n)->getZ()-pt.z;

            b(i,0) = t - (*n)->getTT(nt);
            i++;
        }

        // solve Ax = b with least squares
        x = A.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b);

        //	Eigen::Matrix<T, Eigen::Dynamic, 1> e = b-A*x;

        //	if ( isnan(x[0]) || isnan(x[1]) || isnan(x[2]) ) {
        //		g.x = g.y = g.z = 0;
        //	}

        g.x = x[0];
        g.y = x[1];
        g.z = x[2];

        return g;
    }


    /**
     * Compute traveltime gradient for a tetrahedron, with second-order least-squares
     *
     * @tparam T underlying type of node objects
     * @tparam NODE node objects making the mesh
     */
    template <typename T, typename NODE>
    class Grad3D_ls_so : public Grad3D<T, NODE> {
    public:
        Grad3D_ls_so() = default;
        ~Grad3D_ls_so() = default;

        sxyz<T> compute(const sxyz<T> &pt,
                        const T t,
                        const std::set<NODE*> &nodes,
                        const size_t nt) override;

    private:
        sxyz<T> g;
        Eigen::Matrix<T, Eigen::Dynamic, 9> A;
        Eigen::Matrix<T, 9, 1> x;
        Eigen::Matrix<T, Eigen::Dynamic, 1> b;
    };


    template <typename T, typename NODE>
    sxyz<T> Grad3D_ls_so<T,NODE>::compute(const sxyz<T> &pt,
                                          const T t,
                                          const std::set<NODE*> &nodes,
                                          const size_t nt) {
        // evaluate gradient are center of gravity
        assert(nodes.size()>=9);

        A.resize( nodes.size(), 9 );
        b.resize( nodes.size(), 1 );

        size_t i=0;
        for ( auto n=nodes.cbegin(); n!=nodes.cend(); ++n ) {
            T dx = (*n)->getX()-pt.x;
            T dy = (*n)->getY()-pt.y;
            T dz = (*n)->getZ()-pt.z;

            A(i,0) = dx;
            A(i,1) = dy;
            A(i,2) = dz;
            A(i,3) = 0.5*dx*dx;
            A(i,4) = 0.5*dy*dy;
            A(i,5) = 0.5*dz*dz;
            A(i,6) = dx*dy;
            A(i,7) = dx*dz;
            A(i,8) = dy*dz;

            b(i,0) = t - (*n)->getTT(nt);
            i++;
        }

        // solve Ax = b with least squares
        x = A.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b);

        //	Eigen::Matrix<T, Eigen::Dynamic, 1> e = b-A*x;

        g.x = x[0];
        g.y = x[1];
        g.z = x[2];

        return g;
    }


    /**
     * Compute traveltime gradient for a tetrahedron, with the Averaging-based method
     *
     * @tparam T underlying type of node objects
     * @tparam NODE node objects making the mesh
     */
    template <typename T, typename NODE>
    class Grad3D_ab : public Grad3D<T, NODE> {
    public:
        Grad3D_ab() = default;
        ~Grad3D_ab() = default;

        /**
         * Compute gradient at a given point
         *
         * @param pt point where to compute gradient
         * @param ref_pt reference points for tetrahedra connected to pt
         * @param opp_pts points opposed to reference points in corresponding tetrahedra
         * @param nt thread number
         * @returns value of the travetime gradient
         */
        sxyz<T> compute(const sxyz<T> &pt,
                        const std::vector<NODE*>& ref_pt,
                        const std::vector<std::vector<std::array<NODE*,3>>>& opp_pts,
                        const size_t nt);

        sxyz<T> compute(const sxyz<T> &pt,
                        const T t,
                        const std::set<NODE*> &nodes,
                        const size_t nt) override {
            return {0.0, 0.0, 0.0};   // should never be called
        }

    private:
        Eigen::Matrix<T, 3, 3> A;
        Eigen::Matrix<T, 3, 1> x;
        Eigen::Matrix<T, 3, 1> b;

        sxyz<T> solve(const NODE* n0,
                      const NODE* n1,
                      const NODE* n2,
                      const NODE* n3,
                      const size_t nt);
    };

    template <typename T, typename NODE>
    sxyz<T> Grad3D_ab<T,NODE>::compute(const sxyz<T> &pt,
                                       const std::vector<NODE*>& ref_pt,
                                       const std::vector<std::vector<std::array<NODE*,3>>>& opp_pts,
                                       const size_t nt) {
        sxyz<T> g = {0.0, 0.0, 0.0};

        if ( ref_pt.size() == 1 ) {
            T sum_wi = 0.0;
            for ( size_t n=0; n<opp_pts[0].size(); ++n ) {
                sxyz<T> centroid = { static_cast<T>(0.25) * (ref_pt[0]->getX() + opp_pts[0][n][0]->getX() + opp_pts[0][n][1]->getX() + opp_pts[0][n][2]->getX()),
                    static_cast<T>(0.25) * (ref_pt[0]->getY() + opp_pts[0][n][0]->getY() + opp_pts[0][n][1]->getY() + opp_pts[0][n][2]->getY()),
                    static_cast<T>(0.25) * (ref_pt[0]->getZ() + opp_pts[0][n][0]->getZ() + opp_pts[0][n][1]->getZ() + opp_pts[0][n][2]->getZ()) };
                T wi = 1.0 / ref_pt[0]->getDistance(centroid);
                sum_wi += wi;
                g += wi * solve(ref_pt[0], opp_pts[0][n][0], opp_pts[0][n][1], opp_pts[0][n][2], nt);
            }
            g /= sum_wi;
        } else if ( ref_pt.size() == 2 ) {

            T AB = ref_pt[0]->getDistance( *(ref_pt[1]) );
            std::array<T,2> weights = { ref_pt[1]->getDistance(pt)/AB,
                ref_pt[0]->getDistance(pt)/AB };

            for ( size_t nr=0; nr<ref_pt.size(); ++nr ) {
                T sum_wi = 0.0;
                sxyz<T> g2 = {0.0, 0.0, 0.0};
                for ( size_t n=0; n<opp_pts[nr].size(); ++n ) {
                    sxyz<T> centroid = { static_cast<T>(0.25) * (ref_pt[nr]->getX() + opp_pts[nr][n][0]->getX() + opp_pts[nr][n][1]->getX() + opp_pts[nr][n][2]->getX()),
                        static_cast<T>(0.25) * (ref_pt[nr]->getY() + opp_pts[nr][n][0]->getY() + opp_pts[nr][n][1]->getY() + opp_pts[nr][n][2]->getY()),
                        static_cast<T>(0.25) * (ref_pt[nr]->getZ() + opp_pts[nr][n][0]->getZ() + opp_pts[nr][n][1]->getZ() + opp_pts[nr][n][2]->getZ()) };
                    T wi = 1.0 / ref_pt[nr]->getDistance(centroid);
                    sum_wi += wi;
                    g2 += wi * solve(ref_pt[nr], opp_pts[nr][n][0], opp_pts[nr][n][1], opp_pts[nr][n][2], nt);
                }
                g2 /= sum_wi;
                g += weights[nr] * g2;
            }

        } else if ( ref_pt.size() == 3 ) {

            std::array<T,3> weights;
            Interpolator<T>::bilinearTriangleWeight(pt, ref_pt[0], ref_pt[1], ref_pt[2], weights);

            for ( size_t nr=0; nr<ref_pt.size(); ++nr ) {
                T sum_wi = 0.0;
                sxyz<T> g2 = {0.0, 0.0, 0.0};
                for ( size_t n=0; n<opp_pts[nr].size(); ++n ) {
                    sxyz<T> centroid = { static_cast<T>(0.25) * (ref_pt[nr]->getX() + opp_pts[nr][n][0]->getX() + opp_pts[nr][n][1]->getX() + opp_pts[nr][n][2]->getX()),
                        static_cast<T>(0.25) * (ref_pt[nr]->getY() + opp_pts[nr][n][0]->getY() + opp_pts[nr][n][1]->getY() + opp_pts[nr][n][2]->getY()),
                        static_cast<T>(0.25) * (ref_pt[nr]->getZ() + opp_pts[nr][n][0]->getZ() + opp_pts[nr][n][1]->getZ() + opp_pts[nr][n][2]->getZ()) };
                    T wi = 1.0 / ref_pt[nr]->getDistance(centroid);
                    sum_wi += wi;
                    g2 += wi * solve(ref_pt[nr], opp_pts[nr][n][0], opp_pts[nr][n][1], opp_pts[nr][n][2], nt);
                }
                g2 /= sum_wi;
                g += weights[nr] * g2;
            }

        } else if ( ref_pt.size() == 4 ) {

            std::array<T,4> weights;
            Interpolator<T>::trilinearTriangleWeight(pt, ref_pt[0], ref_pt[1], ref_pt[2], ref_pt[3], weights);

            for ( size_t nr=0; nr<ref_pt.size(); ++nr ) {
                T sum_wi = 0.0;
                sxyz<T> g2 = {0.0, 0.0, 0.0};
                for ( size_t n=0; n<opp_pts[nr].size(); ++n ) {
                    sxyz<T> centroid = { static_cast<T>(0.25) * (ref_pt[nr]->getX() + opp_pts[nr][n][0]->getX() + opp_pts[nr][n][1]->getX() + opp_pts[nr][n][2]->getX()),
                        static_cast<T>(0.25) * (ref_pt[nr]->getY() + opp_pts[nr][n][0]->getY() + opp_pts[nr][n][1]->getY() + opp_pts[nr][n][2]->getY()),
                        static_cast<T>(0.25) * (ref_pt[nr]->getZ() + opp_pts[nr][n][0]->getZ() + opp_pts[nr][n][1]->getZ() + opp_pts[nr][n][2]->getZ()) };
                    T wi = 1.0 / ref_pt[nr]->getDistance(centroid);
                    sum_wi += wi;
                    g2 += wi * solve(ref_pt[nr], opp_pts[nr][n][0], opp_pts[nr][n][1], opp_pts[nr][n][2], nt);
                }
                g2 /= sum_wi;
                g += weights[nr] * g2;
            }

        }

        return g;
    }

    template <typename T, typename NODE>
    sxyz<T> Grad3D_ab<T,NODE>::solve(const NODE* n0,
                                     const NODE* n1,
                                     const NODE* n2,
                                     const NODE* n3,
                                     const size_t nt) {
        A(0,0) = n1->getX()-n0->getX();
        A(0,1) = n1->getY()-n0->getY();
        A(0,2) = n1->getZ()-n0->getZ();
        b(0,0) = n0->getTT(nt)-n1->getTT(nt);

        A(1,0) = n2->getX()-n0->getX();
        A(1,1) = n2->getY()-n0->getY();
        A(1,2) = n2->getZ()-n0->getZ();
        b(1,0) = n0->getTT(nt)-n2->getTT(nt);

        A(2,0) = n3->getX()-n0->getX();
        A(2,1) = n3->getY()-n0->getY();
        A(2,2) = n3->getZ()-n0->getZ();
        b(2,0) = n0->getTT(nt)-n3->getTT(nt);

        x = A.colPivHouseholderQr().solve(b);

        return {x[0], x[1], x[2]};
    }
}

#endif
