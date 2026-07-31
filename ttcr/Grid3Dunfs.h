//
//  Grid3Dunfs.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-04-21.
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

/**
 * @file Grid3Dunfs.h
 * @brief Fast-sweeping solver on a 3-D tetrahedral mesh with node slowness.
 *
 * Declares ttcr::Grid3Dunfs, which solves the eikonal equation by Gauss-Seidel
 * iteration with alternating sweep orders (Qian et al., 2007). It is the
 * node-slowness counterpart of ttcr::Grid3Ducfs and the 3-D counterpart of
 * ttcr::Grid2Dunfs.
 *
 * @section g3dunfs_ordering Sweep orderings
 * On a rectilinear grid the sweeps run along the axes, and the @f$2^3@f$
 * combinations of directions between them cover every possible ray direction.
 * An unstructured mesh has no axes, so the orderings must be built: the nodes
 * are sorted by their distance from each of a set of reference points, and each
 * such ordering is swept forwards and backwards. The reference points are
 * usually the corners of the model's bounding box.
 *
 * ttcr::Grid3Dunfs::initOrdering builds them, using a ttcr::Metric to measure
 * distance — @f$L_1@f$ or @f$L_2@f$, selected by
 * ttcr::input_parameters::order. Since it may be called after construction, the
 * class offers a constructor that takes the reference points and one that does
 * not.
 *
 * @section g3dunfs_converge Convergence
 * Sweeping repeats until the summed change in traveltime falls below a
 * tolerance or @c nitermax iterations have run. The tolerance given to the
 * constructor is **per node**; it is multiplied by the node count so that the
 * test can be applied to the @f$L_1@f$ sum directly, which keeps the criterion
 * independent of mesh size.
 *
 * @sa Grid3Dun.h, Grid3Ducfs.h, Grid2Dunfs.h, Grid3Dunfm.h, Metric.h
 */

#ifndef ttcr_Grid3Dunfs_h
#define ttcr_Grid3Dunfs_h

#include <cmath>
#include <fstream>
#include <queue>
#include <vector>

#include "Grid3Dun.h"
#include "Node3Dn.h"
#include "Metric.h"

namespace ttcr {

    /**
     * @brief Fast-sweeping traveltimes on a tetrahedral mesh with node slowness.
     *
     * @tparam T1 floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2 unsigned integer type used for node and cell indices.
     *
     * Uses ttcr::Node3Dn: no graph adjacency is needed, only the mesh
     * connectivity, so the lighter node type suffices. There are no secondary
     * nodes either — accuracy comes from iterating rather than from refining the
     * discretisation.
     *
     * @sa Grid3Dun.h, Grid3Dunfm.h, Grid3Ducfs.h
     */
    template<typename T1, typename T2>
    class Grid3Dunfs : public Grid3Dun<T1,T2,Node3Dn<T1,T2>> {
    public:
        /**
         * @brief Build the mesh, leaving the sweep orderings unset.
         * @param no    primary node coordinates.
         * @param tet   tetrahedra, as quadruples of indices into @p no.
         * @param eps   convergence tolerance **per node**; see
         *              @ref g3dunfs_converge.
         * @param maxit maximum number of sweep iterations.
         * @param rp    raypath gradient estimator; see @ref g3dun_raypath.
         * @param iv    interpolate velocity rather than slowness.
         * @param rptt  compute traveltimes by integrating along the raypath.
         * @param md    minimum distance for snapping a raypath point onto a node.
         * @param nt    number of threads.
         * @param _translateOrigin shift the mesh to the origin.
         *
         * @pre ttcr::Grid3Dunfs::initOrdering must be called before raytracing;
         *      with no orderings there is nothing to sweep.
         */
        Grid3Dunfs(const std::vector<sxyz<T1>>& no,
                   const std::vector<tetrahedronElem<T2>>& tet,
                   const T1 eps,
                   const int maxit,
                   const int rp,
                   const bool iv,
                   const bool rptt,
                   const T1 md,
                   const size_t nt=1,
                   const bool _translateOrigin=false) :
        Grid3Dun<T1,T2,Node3Dn<T1,T2>>(no, tet, rp, iv, rptt, md, nt, _translateOrigin),
        epsilon(eps), nitermax(maxit), S(), niter_final(0)
        {
            this->buildGridNodes(no, nt);
            this->template buildGridNeighbors<Node3Dn<T1,T2>>(this->nodes);
            epsilon *= static_cast<T1>(this->nodes.size());  // per-node tol -> L1-sum threshold (nodes built)
        }
        /**
         * @brief Build the mesh and its sweep orderings.
         * @param no     primary node coordinates.
         * @param tet    tetrahedra.
         * @param eps    convergence tolerance per node.
         * @param maxit  maximum number of sweep iterations.
         * @param refPts reference points the orderings are built from; see
         *               @ref g3dunfs_ordering.
         * @param order  1 for the @f$L_1@f$ metric, otherwise @f$L_2@f$.
         * @param rp     raypath gradient estimator.
         * @param iv     interpolate velocity rather than slowness.
         * @param rptt   compute traveltimes by integrating along the raypath.
         * @param md     minimum distance for snapping a raypath point.
         * @param nt     number of threads.
         * @param _translateOrigin shift the mesh to the origin.
         */
        Grid3Dunfs(const std::vector<sxyz<T1>>& no,
                   const std::vector<tetrahedronElem<T2>>& tet,
                   const T1 eps,
                   const int maxit,
                   const std::vector<sxyz<T1>>& refPts,
                   const int order,
                   const int rp,
                   const bool iv,
                   const bool rptt,
                   const T1 md,
                   const size_t nt=1,
                   const bool _translateOrigin=false) :
        Grid3Dun<T1,T2,Node3Dn<T1,T2>>(no, tet, rp, iv, rptt, md, nt, _translateOrigin),
        epsilon(eps), nitermax(maxit), S(), niter_final(0)
        {
            this->buildGridNodes(no, nt);
            this->buildGridNeighbors(this->nodes);
            initOrdering(refPts, order);
            epsilon *= static_cast<T1>(this->nodes.size());  // per-node tol -> L1-sum threshold (nodes built)
        }

        ~Grid3Dunfs() {
        }

        /**
         * @brief Build the sweep orderings from a set of reference points.
         * @param refPts the reference points, typically the corners of the
         *               bounding box.
         * @param order  1 for the @f$L_1@f$ metric, otherwise @f$L_2@f$.
         *
         * Sorts the nodes by distance from each reference point, giving one
         * ordering per point. Each is swept forwards and backwards, so
         * @p refPts.size() points yield twice that many sweep directions. See
         * @ref g3dunfs_ordering.
         */
        void initOrdering(const std::vector<sxyz<T1>>& refPts, const int order);

        /**
         * @brief Number of sweep iterations the last solve took.
         * @return the count; equal to @c nitermax if it did not converge.
         */
        const int get_niter() const { return niter_final; }

    private:
        T1 epsilon;    ///< Convergence threshold on the @f$L_1@f$ change, already scaled by the node count.
        int nitermax;  ///< Iteration cap.
        std::vector<std::vector<Node3Dn<T1,T2>*>> S; ///< One node ordering per reference point.
        mutable int niter_final; ///< Iterations taken; @c mutable so the @c const solver can record it.

        /**
         * @brief Set the traveltimes at the source nodes and freeze them.
         * @param Tx           source coordinates.
         * @param t0           source excitation times.
         * @param[out] frozen  flags marking nodes that must not be updated.
         * @param threadNo     thread number.
         *
         * Also seeds the nodes around each source, out to @c source_radius; see
         * ttcr::Grid3Dun::setSourceRadius.
         */
        void initTx(const std::vector<sxyz<T1>>& Tx, const std::vector<T1>& t0,
                    std::vector<bool>& frozen, const size_t threadNo) const;

        /**
         * @brief Solve the traveltime field by sweeping.
         * @param Tx       source coordinates, already origin-translated.
         * @param t0       source excitation times.
         * @param Rx       receiver coordinates; checked to lie in the mesh, but
         *                 they do not bound the work — the whole field is
         *                 solved.
         * @param threadNo thread number.
         *
         * Overrides the base class's solver hook, which the public @c raytrace
         * overloads call; see ttcr::Grid3D. Sweeps each ordering forwards then
         * backwards, repeating until the change falls below @c epsilon or
         * @c nitermax is reached.
         */
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<sxyz<T1>>& Rx,
                      const size_t threadNo=0) const;

        /**
         * @brief Solver hook, grouped-receiver form.
         * @param Tx       source coordinates.
         * @param t0       source excitation times.
         * @param Rx       receiver coordinates, one group per source.
         * @param threadNo thread number.
         */
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      const size_t threadNo=0) const;

    };

    template<typename T1, typename T2>
    void Grid3Dunfs<T1,T2>::initOrdering(const std::vector<sxyz<T1>>& refPts,
                                         const int order) {
        S.resize( refPts.size() );

        Metric<T1> *m;
        if ( order == 1 )
            m = new Metric1<T1>();
        else
            m = new Metric2<T1>();

        std::priority_queue<siv<T1>,std::vector<siv<T1>>,CompareSiv_vr<T1>> queue;

        for ( size_t np=0; np<refPts.size(); ++np ) {

            for ( size_t n=0; n<this->nodes.size(); ++n ) {
                queue.push( {n, m->l(this->nodes[n], refPts[np])} );
            }

            while ( !queue.empty() ) {
                siv<T1> s = queue.top();
                queue.pop();
                S[np].push_back( &(this->nodes[s.i]) );
            }
        }

        delete m;
    }

    template<typename T1, typename T2>
    void Grid3Dunfs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<sxyz<T1>>& Rx,
                                     const size_t threadNo) const {

        this->checkPts(Tx, true);
        this->checkPts(Rx, true);

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        std::vector<bool> frozen( this->nodes.size(), false );
        initTx(Tx, t0, frozen, threadNo);

        std::vector<T1> times( this->nodes.size() );
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            times[n] = this->nodes[n].getTT( threadNo );
        }

        int niter = 0;
        T1 change = std::numeric_limits<T1>::max();
        while ( change >= epsilon && niter<nitermax ) {

            for ( size_t i=0; i<S.size(); ++i ) {

                // ascending
                for ( auto vertexC=S[i].begin(); vertexC!=S[i].end(); ++vertexC ) {
                    if ( !frozen[(*vertexC)->getGridIndex()] )
                        this->localUpdate3D(*vertexC, threadNo);
                }

                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                if ( change < epsilon ) {
                    break;
                }

                // descending
                for ( auto vertexC=S[i].rbegin(); vertexC!=S[i].rend(); ++vertexC ) {
                    if ( !frozen[(*vertexC)->getGridIndex()] )
                        this->localUpdate3D(*vertexC, threadNo);
                }

                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                if ( change < epsilon ) {
                    break;
                }
            }
            niter++;
        }
        niter_final = niter;
    }

    template<typename T1, typename T2>
    void Grid3Dunfs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<std::vector<sxyz<T1>>>& Rx,
                                     const size_t threadNo) const {

        this->checkPts(Tx, true);
        for ( size_t n=0; n<Rx.size(); ++n ) {
            this->checkPts(Rx[n], true);
        }

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        std::vector<bool> frozen( this->nodes.size(), false );
        initTx(Tx, t0, frozen, threadNo);

        std::vector<T1> times( this->nodes.size() );
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            times[n] = this->nodes[n].getTT( threadNo );
        }

        int niter = 0;
        T1 change = std::numeric_limits<T1>::max();
        while ( change >= epsilon && niter<nitermax ) {

            for ( size_t i=0; i<S.size(); ++i ) {

                // ascending
                for ( auto vertexC=S[i].begin(); vertexC!=S[i].end(); ++vertexC ) {
                    if ( !frozen[(*vertexC)->getGridIndex()] )
                        this->localUpdate3D(*vertexC, threadNo);
                }

                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                if ( change < epsilon ) {
                    break;
                }

                // descending
                for ( auto vertexC=S[i].rbegin(); vertexC!=S[i].rend(); ++vertexC ) {
                    if ( !frozen[(*vertexC)->getGridIndex()] )
                        this->localUpdate3D(*vertexC, threadNo);
                }

                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                if ( change < epsilon ) {
                    break;
                }

            }
            niter++;
        }
        niter_final = niter;
    }


    template<typename T1, typename T2>
    void Grid3Dunfs<T1,T2>::initTx(const std::vector<sxyz<T1>>& Tx,
                                   const std::vector<T1>& t0,
                                   std::vector<bool>& frozen,
                                   const size_t threadNo) const {

        for (size_t n=0; n<Tx.size(); ++n) {
            bool found = false;
            for ( size_t nn=0; nn<this->nodes.size(); ++nn ) {
                if ( this->nodes[nn] == Tx[n] ) {
                    found = true;
                    this->nodes[nn].setTT( t0[n], threadNo );
                    frozen[nn] = true;

                    if ( Grid3Dun<T1,T2,Node3Dn<T1,T2>>::source_radius == 0.0 ) {
                        // populate around Tx
                        for ( size_t no=0; no<this->nodes[nn].getOwners().size(); ++no ) {

                            T2 cellNo = this->nodes[nn].getOwners()[no];
                            for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                                T2 neibNo = this->neighbors[cellNo][k];
                                if ( neibNo == nn ) continue;
                                T1 dt = this->computeDt(this->nodes[nn], this->nodes[neibNo]);

                                if ( t0[n]+dt < this->nodes[neibNo].getTT(threadNo) ) {
                                    this->nodes[neibNo].setTT( t0[n]+dt, threadNo );
                                }
                            }
                        }
                    } else {
                        // find nodes within source radius
                        size_t nodes_added = 0;
                        for ( size_t no=0; no<this->nodes.size(); ++no ) {

                            if ( no == nn ) continue;

                            T1 d = this->nodes[nn].getDistance( this->nodes[no] );
                            if ( d <= Grid3Dun<T1,T2,Node3Dn<T1,T2>>::source_radius ) {

                                T1 dt = this->computeDt(this->nodes[nn], this->nodes[no] );

                                if ( t0[n]+dt < this->nodes[no].getTT(threadNo) ) {
                                    if ( this->nodes[no].getTT(threadNo) == std::numeric_limits<T1>::max() ) nodes_added++;
                                    this->nodes[no].setTT( t0[n]+dt, threadNo );
                                }
                            }
                        }
                        if ( nodes_added == 0 ) {
                            std::cerr << "Error: no nodes found within source radius, aborting" << std::endl;
                            abort();
                        } else {
                            std::cout << "(found " << nodes_added << " nodes around Tx point)\n";
                        }
                    }

                    break;
                }
            }
            if ( found==false ) {

                T2 sTx = this->computeSlowness( Tx[n], true );

                T2 cellNo = this->getCellNo(Tx[n]);
                if ( Grid3Dun<T1,T2,Node3Dn<T1,T2>>::source_radius == 0.0 ) {
                    for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                        T2 neibNo = this->neighbors[cellNo][k];
                        // compute dt
                        T1 dt = this->computeDt(this->nodes[neibNo], Tx[n], sTx);

                        this->nodes[neibNo].setTT( t0[n]+dt, threadNo );
                        frozen[neibNo] = true;
                    }
                } else if ( Tx.size()==1 ) { // look into source radius only for point sources
                    // find nodes within source radius
                    size_t nodes_added = 0;
                    for ( size_t no=0; no<this->nodes.size(); ++no ) {

                        T1 d = this->nodes[no].getDistance( Tx[n] );
                        if ( d <= Grid3Dun<T1,T2,Node3Dn<T1,T2>>::source_radius ) {

                            T1 dt = this->computeDt(this->nodes[no], Tx[n], sTx);

                            if ( t0[n]+dt < this->nodes[no].getTT(threadNo) ) {
                                if ( this->nodes[no].getTT(threadNo) == std::numeric_limits<T1>::max() ) nodes_added++;
                                this->nodes[no].setTT( t0[n]+dt, threadNo );
                            }
                        }
                    }
                    if ( nodes_added == 0 ) {
                        std::cerr << "Error: no nodes found within source radius, aborting" << std::endl;
                        abort();
                    } else {
                        std::cout << "(found " << nodes_added << " nodes around Tx point)\n";
                    }
                }
            }
        }
    }

}

#endif
