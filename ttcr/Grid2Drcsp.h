/*
 *  Grid2Drcsp.h
 *  ttcr
 *
 *  Created by Bernard Giroux on 08-04-24.
 *  Copyright 2008 Bernard Giroux.
 */

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

/*
 * Reference paper
 *
 * @article{gruber:1062,
 *  author = {Thomas Gruber and Stewart A. Greenhalgh},
 *  collaboration = {},
 *  title = {Precision analysis of first-break times in grid models},
 *  publisher = {SEG},
 *  year = {1998},
 *  journal = {Geophysics},
 *  volume = {63},
 *  number = {3},
 *  pages = {1062-1065},
 *  url = {http://link.aip.org/link/?GPY/63/1062/1},
 *  doi = {10.1190/1.1444384}
 * }
 *
 *
 */

/**
 * @file Grid2Drcsp.h
 * @brief Shortest-path solver on a 2-D rectilinear grid with cell-based slowness.
 *
 * Declares ttcr::Grid2Drcsp, which solves the eikonal equation by treating the
 * grid as a graph and running a Dijkstra-style relaxation over it. It derives
 * from ttcr::Grid2Drc and adds the two things the shortest-path method needs:
 * **secondary nodes** along the cell edges, and a priority-queue propagation.
 *
 * @section g2drcsp_secondary Secondary nodes
 * With nodes only at the cell corners, a graph method can only propagate along
 * the few directions those corners allow, and the traveltime error is dominated
 * by that angular discretisation. The remedy is to add @c nsnx nodes along each
 * horizontal cell edge and @c nsnz along each vertical one, which multiplies the
 * number of available ray directions and drives the error down — at a cost in
 * memory and time that grows quickly with the counts. Accuracy as a function of
 * these counts is the subject of the reference paper above.
 *
 * The secondary nodes are appended after the primary ones in the node vector,
 * as described in @ref g2drc_numbering, so
 * @ref ttcr::Grid2Drcsp::getTT still returns a plain
 * @f$(n_{cx}+1)\times(n_{cz}+1)@f$ field.
 *
 * @sa Grid2Drc.h, Node2Dcsp.h, Grid2Drcdsp.h
 */

#ifndef ttcr_Grid2Drcsp_h
#define ttcr_Grid2Drcsp_h

#include "Grid2Drc.h"
#include "Node2Dcsp.h"

namespace ttcr {

    /**
     * @brief Whether a cell class reports its sensitivity into a given container
     *
     * The number of values a cell reports is the number of medium parameters it
     * describes, so each cell class provides computeDistance() for one container
     * only.  A grid, however, has to override every raytrace() of Grid2D,
     * whichever cells it holds, and the overrides are emitted with the vtable
     * whether or not they are ever called.  This tells them apart, so that the
     * ones the cells cannot serve fail when called rather than when compiled.
     */
    template<typename CELL, typename SIV, typename NODE, typename S,
             typename = void>
    struct cell_reports_into : std::false_type {};

    template<typename CELL, typename SIV, typename NODE, typename S>
    struct cell_reports_into<CELL, SIV, NODE, S,
        std::void_t<decltype(std::declval<const CELL&>().computeDistance(
            std::declval<const NODE&>(), std::declval<const S&>(),
            std::declval<SIV&>()))>> : std::true_type {};


    /**
     * @brief Shortest-path eikonal solver on a rectilinear cell-slowness grid.
     *
     * @tparam T1   floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2   integer type of node and cell indices.
     * @tparam S    point type, @ref sxz or @ref sxyz.
     * @tparam CELL cell policy from Cell.h; this class stays templated on it, so
     *              it supports the anisotropic models as well as the isotropic
     *              one.
     *
     * Uses ttcr::Node2Dcsp, whose parent-node and parent-cell arrays let a
     * raypath be recovered by walking back from a receiver to the source.
     */
    template<typename T1, typename T2, typename S, typename CELL>
    class Grid2Drcsp : public Grid2Drc<T1,T2,S,Node2Dcsp<T1,T2>,CELL> {
    public:
        /**
         * @brief Build the grid and its secondary nodes.
         *
         * @param nx  number of cells along x.
         * @param nz  number of cells along z.
         * @param ddx cell size along x.
         * @param ddz cell size along z.
         * @param minx x coordinate of the grid origin.
         * @param minz z coordinate of the grid origin.
         * @param nnx number of secondary nodes per cell edge along x.
         * @param nnz number of secondary nodes per cell edge along z.
         * @param ttrp recompute receiver traveltimes along the raypath.
         * @param nt  number of threads.
         *
         * @note Raising @p nnx and @p nnz improves accuracy but costs memory and
         *       time; see @ref g2drcsp_secondary. They come from
         *       ttcr::input_parameters::nn.
         */
        Grid2Drcsp(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
                   const T1 minx, const T1 minz, const T2 nnx, const T2 nnz,
                   const bool ttrp, const size_t nt=1) :
        Grid2Drc<T1,T2,S,Node2Dcsp<T1,T2>,CELL>(nx,nz,ddx,ddz,minx,minz,ttrp,nt),
        nsnx(nnx), nsnz(nnz), nsgx(0), nsgz(0)
        {
            buildGridNodes();
            this->template buildGridNeighbors<Node2Dcsp<T1,T2>>(this->nodes);
        }

        /// Destructor.
        virtual ~Grid2Drcsp() {
        }

        /**
         * @name Raytracing
         *
         * Seven overloads of the same computation, differing only in how much
         * they report. All of them propagate the traveltime field from @p Tx
         * over the whole grid, then evaluate it at the receivers; the extra
         * output arguments are recovered afterwards by walking the parent
         * pointers of ttcr::Node2Dcsp back from each receiver.
         *
         * Two axes of variation:
         * - **Receiver grouping** — a plain @c Rx vector, or a vector of
         *   pointers to receiver vectors, which lets one propagation serve
         *   several receiver sets (used for reflected arrivals, where each
         *   reflector is its own set).
         * - **What is returned** — traveltimes alone; plus @c r_data, the
         *   raypath geometry; plus @c l_data, the length travelled in each cell,
         *   which forms one row of the sensitivity matrix.
         *
         * @c l_data comes in two flavours: @ref siv holds one value per cell and
         * @ref siv2 holds two, the second being used by the anisotropic cell
         * policies that need separate along- and across-axis path lengths.
         *
         * @param[in]  Tx          source positions.
         * @param[in]  t0          origin time of each source, parallel to @p Tx.
         * @param[in]  Rx          receiver positions.
         * @param[out] traveltimes traveltime at each receiver.
         * @param[in]  threadNo    thread to compute on; selects which per-node
         *                         traveltime slot is written, so concurrent
         *                         calls must pass distinct values.
         *
         * @pre @p Tx and @p t0 have the same length, and every point lies inside
         *      the grid — @ref checkPts is called and throws otherwise.
         * @note @c const, but they do mutate the node traveltimes for
         *       @p threadNo; that is what @c mutable on @ref nodes is for.
         * @{
         */

        /// Traveltimes only.
        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<S>& Rx,
                      std::vector<T1>& traveltimes,
                      const size_t threadNo=0) const;

        /// Traveltimes, for several receiver sets from one propagation.
        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<const std::vector<S>*>& Rx,
                      std::vector<std::vector<T1>*>& traveltimes,
                      const size_t threadNo=0) const;

        /// Traveltimes and raypaths.
        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<S>& Rx,
                      std::vector<T1>& traveltimes,
                      std::vector<std::vector<S>>& r_data,
                      const size_t threadNo=0) const;

        /// Traveltimes and raypaths, for several receiver sets.
        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<const std::vector<S>*>& Rx,
                      std::vector<std::vector<T1>*>& traveltimes,
                      std::vector<std::vector<std::vector<S>>*>& r_data,
                      const size_t threadNo=0) const;

        /// @name Traveltimes and per-cell sensitivities
        ///
        /// The cells report the sensitivity of the traveltime to their medium
        /// parameters, one value per parameter, so the container widens with
        /// the number of parameters the medium has: @ref siv for an isotropic
        /// cell, @ref siv2, @ref siv4 or @ref siv5 for the anisotropic ones.
        /// All of these walk the same shortest-path tree, through
        /// raytraceSensitivity().
        /// @{
        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<S>& Rx,
                      std::vector<T1>& traveltimes,
                      std::vector<std::vector<S>>& r_data,
                      std::vector<std::vector<siv<T1>>>& l_data,
                      const size_t threadNo) const {
            raytraceSensitivity(Tx, t0, Rx, traveltimes, &r_data, l_data, threadNo);
        }

        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<S>& Rx,
                      std::vector<T1>& traveltimes,
                      std::vector<std::vector<S>>& r_data,
                      std::vector<std::vector<siv2<T1>>>& l_data,
                      const size_t threadNo=0) const {
            raytraceSensitivity(Tx, t0, Rx, traveltimes, &r_data, l_data, threadNo);
        }

        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<S>& Rx,
                      std::vector<T1>& traveltimes,
                      std::vector<std::vector<S>>& r_data,
                      std::vector<std::vector<siv4<T1>>>& l_data,
                      const size_t threadNo=0) const {
            raytraceSensitivity(Tx, t0, Rx, traveltimes, &r_data, l_data, threadNo);
        }

        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<S>& Rx,
                      std::vector<T1>& traveltimes,
                      std::vector<std::vector<S>>& r_data,
                      std::vector<std::vector<siv5<T1>>>& l_data,
                      const size_t threadNo=0) const {
            raytraceSensitivity(Tx, t0, Rx, traveltimes, &r_data, l_data, threadNo);
        }

        /// The same, without the raypaths.
        ///
        /// Without these overrides the calls would resolve to
        /// Grid2D::raytrace(), which traces the raypath by following the
        /// gradient of the traveltime field instead of walking the
        /// shortest-path tree, yielding a longer path and a traveltime
        /// inconsistent with the other raytrace() overloads of this class.
        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<S>& Rx,
                      std::vector<T1>& traveltimes,
                      std::vector<std::vector<siv<T1>>>& l_data,
                      const size_t threadNo=0) const {
            raytraceSensitivity(Tx, t0, Rx, traveltimes, nullptr, l_data, threadNo);
        }

        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<S>& Rx,
                      std::vector<T1>& traveltimes,
                      std::vector<std::vector<siv2<T1>>>& l_data,
                      const size_t threadNo=0) const {
            raytraceSensitivity(Tx, t0, Rx, traveltimes, nullptr, l_data, threadNo);
        }

        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<S>& Rx,
                      std::vector<T1>& traveltimes,
                      std::vector<std::vector<siv4<T1>>>& l_data,
                      const size_t threadNo=0) const {
            raytraceSensitivity(Tx, t0, Rx, traveltimes, nullptr, l_data, threadNo);
        }

        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<S>& Rx,
                      std::vector<T1>& traveltimes,
                      std::vector<std::vector<siv5<T1>>>& l_data,
                      const size_t threadNo=0) const {
            raytraceSensitivity(Tx, t0, Rx, traveltimes, nullptr, l_data, threadNo);
        }
        /// @}

        /// @return Number of secondary nodes per cell edge along x.
        const T2 getNsnx() const { return nsnx; }
        /// @return Number of secondary nodes per cell edge along z.
        const T2 getNsnz() const { return nsnz; }

        /**
         * @brief Collect the traveltimes computed at the primary nodes.
         * @param[out] tt       resized to the primary-node count and filled.
         * @param[in]  threadNo thread whose solution to read.
         * @note Overrides ttcr::Grid2Drc::getTT with an identical implementation;
         *       the secondary nodes this class adds are skipped either way.
         */
        void getTT(std::vector<T1>& tt, const size_t threadNo=0) const final {
            size_t nPrimary = static_cast<size_t>(this->ncx+1) * (this->ncz+1);
            tt.resize(nPrimary);
            for ( size_t n=0, n2=0; n<this->nodes.size(); ++n ) {
                if (this->nodes[n].isPrimary()) {
                    tt[n2++] = this->nodes[n].getTT(threadNo);
                }
            }
        }

    private:
        T2 nsnx;    ///< number of secondary nodes in x
        T2 nsnz;    ///< number of secondary nodes in z
        T2 nsgx;    ///< number of subgrid cells in x
        T2 nsgz;    ///< number of subgrid cells in z

        /// @name Non-copyable
        /// Copy assignment is deleted, so any attempt to assign one grid to
        /// another is a compile error. The default and copy constructors still
        /// use the older private-and-defined idiom; they are unreachable from
        /// outside the class, but unlike the assignment operator they would
        /// produce an object with uninitialised members if ever called from
        /// within it.
        /// @{
        Grid2Drcsp() {}
        Grid2Drcsp(const Grid2Drcsp<T1,T2,S,CELL>& g) {}
        Grid2Drcsp<T1,T2,S,CELL>& operator=(const Grid2Drcsp<T1,T2,S,CELL>& g) = delete;
        /// @}

        /**
         * @brief Create the primary and secondary nodes and set their positions.
         * @post @ref nodes holds the cell corners followed by @ref nsnx and
         *       @ref nsnz secondary nodes per edge; every node knows which cells
         *       own it. @sa @ref g2drcsp_secondary
         */
        void buildGridNodes();

        void propagate(std::priority_queue<Node2Dcsp<T1,T2>*,
                       std::vector<Node2Dcsp<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;

        void propagate_lw(std::priority_queue<Node2Dcsp<T1,T2>*,
                          std::vector<Node2Dcsp<T1,T2>*>,
                          CompareNodePtr<T1>>& queue,
                          std::vector<bool>& inQueue,
                          std::vector<bool>& frozen,
                          const size_t threadNo) const;

        void initQueue(const std::vector<S>& Tx,
                       const std::vector<T1>& t0,
                       std::priority_queue<Node2Dcsp<T1,T2>*,
                       std::vector<Node2Dcsp<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<Node2Dcsp<T1,T2>>& txNodes,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;

        void initBand(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      std::priority_queue<Node2Dcsp<T1,T2>*,
                      std::vector<Node2Dcsp<T1,T2>*>,
                      CompareNodePtr<T1>>& queue,
                      std::vector<Node2Dcsp<T1,T2>>& txNodes,
                      std::vector<bool>& inQueue,
                      std::vector<bool>& frozen,
                      const size_t threadNo) const;

        /**
         * @brief Traveltimes and per-cell sensitivities, whatever the container
         *
         * Walks the shortest-path tree from each receiver back to the source,
         * asking the cells for the sensitivity of the traveltime of every
         * segment and accumulating it.  The four public overloads differ only
         * in the container the cells report into, and in whether the raypath is
         * wanted, so they share this one body.
         *
         * @param Tx          source coordinates
         * @param t0          source excitation times
         * @param Rx          receiver coordinates
         * @param traveltimes traveltimes at the receivers
         * @param r_data      raypaths, or nullptr if they are not wanted
         * @param l_data      per-cell sensitivities
         * @param threadNo    thread on which to run
         */
        template<typename SIV>
        void raytraceSensitivity(const std::vector<S>& Tx,
                                 const std::vector<T1>& t0,
                                 const std::vector<S>& Rx,
                                 std::vector<T1>& traveltimes,
                                 std::vector<std::vector<S>>* r_data,
                                 std::vector<std::vector<SIV>>& l_data,
                                 const size_t threadNo) const;

        /// @brief Add the contribution of a segment to the cell it crosses
        ///
        /// A ray may cross the same cell more than once, passing through the
        /// secondary nodes along an edge, so the contributions are summed.
        template<typename SIV>
        static void accumulate(std::vector<SIV>& row, const SIV& cell) {
            for ( size_t nc=0; nc<row.size(); ++nc ) {
                if ( row[nc].i == cell.i ) {
                    row[nc] += cell;
                    return;
                }
            }
            row.push_back( cell );
        }

        T1 get_tt_corr(const siv2<T1>& cell,
                       const Grid2Drcsp<T1,T2,S,CELL> *grid,
                       const size_t i) const {
            return cell.v*(this->slowness[cell.i] - grid->slowness[i]);
        }

        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<S>& Rx,
                      const size_t threadNo) const;

        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<const std::vector<S>*>& Rx,
                      const size_t threadNo=0) const;

        using Grid2Drc<T1,T2,S,Node2Dcsp<T1,T2>,CELL>::getTraveltime;
        
        T1 getTraveltime(const S& Rx, const std::vector<Node2Dcsp<T1,T2>>& nodes,
                         T2& nodeParentRx, T2& cellParentRx,
                         const size_t threadNo) const;

    };

    template<typename T1, typename T2, typename S, typename CELL>
    void Grid2Drcsp<T1,T2,S,CELL>::buildGridNodes() {

        this->nodes.resize( // noeuds secondaires
                           this->ncx*nsnx*(this->ncz+1) +
                           this->ncz*nsnz*(this->ncx+1) +
                           // noeuds primaires
                           (this->ncx+1) * (this->ncz+1),
                           Node2Dcsp<T1,T2>(this->nThreads));

        T1 dxs = this->dx/(nsnx+1);
        T1 dzs = this->dz/(nsnz+1);

        T2 cell_upLeft = std::numeric_limits<T2>::max();
        T2 cell_upRight = std::numeric_limits<T2>::max();
        T2 cell_downLeft = 0;
        T2 cell_downRight = 0;

        T2 n = 0;
        for ( T2 nc=0; nc<=this->ncx; ++nc ) {

            T1 x = this->xmin + nc*this->dx;

            for ( T2 nr=0; nr<=this->ncz; ++nr ) {

                T1 z = this->zmin + nr*this->dz;

                if ( nr < this->ncz && nc < this->ncx ) {
                    cell_downRight = nc*this->ncz + nr;
                }
                else {
                    cell_downRight = std::numeric_limits<T2>::max();
                }

                if ( nr > 0 && nc < this->ncx ) {
                    cell_upRight = nc*this->ncz + nr - 1;
                }
                else {
                    cell_upRight = std::numeric_limits<T2>::max();
                }

                if ( nr < this->ncz && nc > 0 ) {
                    cell_downLeft = (nc-1)*this->ncz + nr;
                }
                else {
                    cell_downLeft = std::numeric_limits<T2>::max();
                }

                if ( nr > 0 && nc > 0 ) {
                    cell_upLeft = (nc-1)*this->ncz + nr - 1;
                }
                else {
                    cell_upLeft = std::numeric_limits<T2>::max();
                }

                if ( cell_upLeft != std::numeric_limits<T2>::max() ) {
                    this->nodes[n].pushOwner( cell_upLeft );
                }
                if ( cell_downLeft != std::numeric_limits<T2>::max() ) {
                    this->nodes[n].pushOwner( cell_downLeft );
                }
                if ( cell_upRight != std::numeric_limits<T2>::max() ) {
                    this->nodes[n].pushOwner( cell_upRight );
                }
                if ( cell_downRight != std::numeric_limits<T2>::max() ) {
                    this->nodes[n].pushOwner( cell_downRight );
                }

                this->nodes[n].setX( x );
                this->nodes[n].setZ( z );
                this->nodes[n].setGridIndex( n );
                this->nodes[n].setPrimary(true);

                ++n;
            }
        }

        // continue with secondary nodes
        for ( T2 nc=0; nc<=this->ncx; ++nc ) {

            T1 x = this->xmin + nc*this->dx;

            for ( T2 nr=0; nr<=this->ncz; ++nr ) {

                T1 z = this->zmin + nr*this->dz;

                if ( nr < this->ncz && nc < this->ncx ) {
                    cell_downRight = nc*this->ncz + nr;
                }
                else {
                    cell_downRight = std::numeric_limits<T2>::max();
                }

                if ( nr > 0 && nc < this->ncx ) {
                    cell_upRight = nc*this->ncz + nr - 1;
                }
                else {
                    cell_upRight = std::numeric_limits<T2>::max();
                }

                if ( nr < this->ncz && nc > 0 ) {
                    cell_downLeft = (nc-1)*this->ncz + nr;
                }
                else {
                    cell_downLeft = std::numeric_limits<T2>::max();
                }

                if ( nr > 0 && nc > 0 ) {
                    cell_upLeft = (nc-1)*this->ncz + nr - 1;
                }
                else {
                    cell_upLeft = std::numeric_limits<T2>::max();
                }

                // secondary nodes on the vertical
                if ( nr < this->ncz ) {
                    for (T2 ns=0; ns<nsnz; ++ns, ++n ) {

                        T1 zsv = this->zmin + nr*this->dz + (ns+1)*dzs;

                        if ( cell_downLeft != std::numeric_limits<T2>::max() ) {
                            this->nodes[n].pushOwner( cell_downLeft );
                        }
                        if ( cell_downRight != std::numeric_limits<T2>::max() ) {
                            this->nodes[n].pushOwner( cell_downRight );
                        }

                        this->nodes[n].setX( x );
                        this->nodes[n].setZ( zsv );
                        this->nodes[n].setGridIndex( n );
                        this->nodes[n].setPrimary(false);
                    }
                }

                // secondary nodes on the horizontal
                if ( nc < this->ncx ) {
                    for ( T2 ns=0; ns<nsnx; ++ns, ++n ) {

                        T1 xsh = this->xmin + nc*this->dx + (ns+1)*dxs;

                        if ( cell_upRight != std::numeric_limits<T2>::max() ) {
                            this->nodes[n].pushOwner( cell_upRight );
                        }
                        if ( cell_downRight != std::numeric_limits<T2>::max() ) {
                            this->nodes[n].pushOwner( cell_downRight );
                        }

                        this->nodes[n].setX( xsh );
                        this->nodes[n].setZ( z );
                        this->nodes[n].setGridIndex( n );
                        this->nodes[n].setPrimary(false);
                    }
                }
            }
        }
    }


    template<typename T1, typename T2, typename S, typename CELL>
    void Grid2Drcsp<T1,T2,S,CELL>::initQueue(const std::vector<S>& Tx,
                                             const std::vector<T1>& t0,
                                             std::priority_queue<Node2Dcsp<T1,T2>*,
                                             std::vector<Node2Dcsp<T1,T2>*>,
                                             CompareNodePtr<T1>>& queue,
                                             std::vector<Node2Dcsp<T1,T2>>& txNodes,
                                             std::vector<bool>& inQueue,
                                             std::vector<bool>& frozen,
                                             const size_t threadNo) const {

        for (size_t n=0; n<Tx.size(); ++n) {
            bool found = false;
            for ( size_t nn=0; nn<this->nodes.size(); ++nn ) {
                if ( this->nodes[nn] == Tx[n] ) {
                    found = true;
                    this->nodes[nn].setTT( t0[n], threadNo );
                    queue.push( &(this->nodes[nn]) );
                    inQueue[nn] = true;
                    frozen[nn] = true;
                    break;
                }
            }
            if ( found==false ) {
                txNodes.push_back( Node2Dcsp<T1,T2>(t0[n], Tx[n].x, Tx[n].z,
                                                    this->nThreads, threadNo) );
                // we belong to cell index no
                txNodes.back().pushOwner( this->getCellNo(Tx[n]) );
                txNodes.back().setGridIndex( static_cast<T2>(this->nodes.size()+
                                                             txNodes.size()-1) );

                queue.push( &(txNodes.back()) );
                inQueue.push_back( true );
                frozen.push_back( true );
            }
        }
    }


    template<typename T1, typename T2, typename S, typename CELL>
    void Grid2Drcsp<T1,T2,S,CELL>::initBand(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            std::priority_queue<Node2Dcsp<T1,T2>*,
                                            std::vector<Node2Dcsp<T1,T2>*>,
                                            CompareNodePtr<T1>>& narrow_band,
                                            std::vector<Node2Dcsp<T1,T2>>& txNodes,
                                            std::vector<bool>& inBand,
                                            std::vector<bool>& frozen,
                                            const size_t threadNo) const {

        for (size_t n=0; n<Tx.size(); ++n) {
            bool found = false;
            for ( size_t nn=0; nn<this->nodes.size(); ++nn ) {
                if ( this->nodes[nn] == Tx[n] ) {
                    found = true;
                    this->nodes[nn].setTT( t0[n], threadNo );
                    narrow_band.push( &(this->nodes[nn]) );
                    inBand[nn] = true;
                    frozen[nn] = true;

                    if ( Tx.size()==1 ) {
                        // populate around Tx
                        for ( size_t no=0; no<this->nodes[nn].getOwners().size(); ++no ) {

                            T2 cellNo = this->nodes[nn].getOwners()[no];
                            for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                                T2 neibNo = this->neighbors[cellNo][k];
                                if ( neibNo == nn ) continue;
                                T1 dt = this->cells.computeDt(this->nodes[nn], this->nodes[neibNo], cellNo);

                                if ( t0[n]+dt < this->nodes[neibNo].getTT(threadNo) ) {
                                    this->nodes[neibNo].setTT( t0[n]+dt, threadNo );
                                    this->nodes[neibNo].setnodeParent(this->nodes[nn].getGridIndex(),threadNo);
                                    this->nodes[neibNo].setCellParent(cellNo, threadNo );

                                    if ( !inBand[neibNo] ) {
                                        narrow_band.push( &(this->nodes[neibNo]) );
                                        inBand[neibNo] = true;
                                        frozen[neibNo] = true;
                                    }
                                }
                            }
                        }
                    }
                    break;
                }
            }
            if ( found==false ) {

                T2 cellNo = getCellNo(Tx[n]);
                for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                    T2 neibNo = this->neighbors[cellNo][k];

                    // compute dt
                    T1 dt = this->cells.computeDt(this->nodes[neibNo], Tx[n], cellNo);

                    this->nodes[neibNo].setTT( t0[n]+dt, threadNo );
                    narrow_band.push( &(this->nodes[neibNo]) );
                    inBand[neibNo] = true;
                    frozen[neibNo] = true;

                }
            }
        }
    }


    template<typename T1, typename T2, typename S, typename CELL>
    void Grid2Drcsp<T1,T2,S,CELL>::raytrace(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const std::vector<S>& Rx,
                                            const size_t threadNo) const {

        this->checkPts(Tx);
        this->checkPts(Rx);

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node2Dcsp<T1,T2>*, std::vector<Node2Dcsp<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );

        std::vector<Node2Dcsp<T1,T2>> txNodes;
        std::vector<bool> inQueue( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );

        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);

        propagate_lw(queue, inQueue, frozen, threadNo);
    }

    template<typename T1, typename T2, typename S, typename CELL>
    void Grid2Drcsp<T1,T2,S,CELL>::raytrace(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const std::vector<const std::vector<S>*>& Rx,
                                            const size_t threadNo) const {
        this->checkPts(Tx);
        for ( size_t n=0; n<Rx.size(); ++n ) {
            this->checkPts(*Rx[n]);
        }

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node2Dcsp<T1,T2>*, std::vector<Node2Dcsp<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );

        std::vector<Node2Dcsp<T1,T2>> txNodes;
        std::vector<bool> inQueue( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );

        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);

        propagate_lw(queue, inQueue, frozen, threadNo);

    }

    template<typename T1, typename T2, typename S, typename CELL>
    void Grid2Drcsp<T1,T2,S,CELL>::raytrace(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const std::vector<S>& Rx,
                                            std::vector<T1>& traveltimes,
                                            const size_t threadNo) const {
        // this to make sure that we are no use tt_from _rp
        raytrace(Tx, t0, Rx, threadNo);

        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        for (size_t n=0; n<Rx.size(); ++n) {
            traveltimes[n] = Grid2Drc<T1,T2,S,Node2Dcsp<T1,T2>,CELL>::getTraveltime(Rx[n], threadNo);
        }
    }

    template<typename T1, typename T2, typename S, typename CELL>
    void Grid2Drcsp<T1,T2,S,CELL>::raytrace(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const std::vector<const std::vector<S>*>& Rx,
                                            std::vector<std::vector<T1>*>& traveltimes,
                                            const size_t threadNo) const {
        // this to make sure that we are no use tt_from _rp
        raytrace(Tx, t0, Rx, threadNo);

        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        for (size_t nr=0; nr<Rx.size(); ++nr) {
            traveltimes[nr]->resize( Rx[nr]->size() );
            for (size_t n=0; n<Rx[nr]->size(); ++n) {
                (*traveltimes[nr])[n] = Grid2Drc<T1,T2,S,Node2Dcsp<T1,T2>,CELL>::getTraveltime((*Rx[nr])[n], threadNo);
            }
        }
    }

    template<typename T1, typename T2, typename S, typename CELL>
    void Grid2Drcsp<T1,T2,S,CELL>::raytrace(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const std::vector<S>& Rx,
                                            std::vector<T1>& traveltimes,
                                            std::vector<std::vector<S>>& r_data,
                                            const size_t threadNo) const {

        this->checkPts(Tx);
        this->checkPts(Rx);

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node2Dcsp<T1,T2>*, std::vector<Node2Dcsp<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );
        std::vector<Node2Dcsp<T1,T2>> txNodes;
        std::vector<bool> inQueue( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );

        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);

        propagate(queue, inQueue, frozen, threadNo);

        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        if ( r_data.size() != Rx.size() ) {
            r_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<r_data.size(); ++ni ) {
            r_data[ni].resize( 0 );
        }
        T2 nodeParentRx;
        T2 cellParentRx;

        for (size_t n=0; n<Rx.size(); ++n) {

            // check if Rx is on one of Tx nodes
            bool foundTx = false;
            for (size_t nt=0; nt<Tx.size(); ++nt) {
                if (Tx[nt] == Rx[n]) {
                    traveltimes[n] = t0[nt];
                    r_data[n].push_back(Tx[nt]);
                    foundTx = true;
                    break;
                }
            }
            if (foundTx)
                continue;

            traveltimes[n] = getTraveltime(Rx[n], this->nodes, nodeParentRx, cellParentRx,
                                           threadNo);

            // Rx are in nodes (not txNodes)
            std::vector<Node2Dcsp<T1,T2>> *node_p;
            node_p = &(this->nodes);

            std::vector<sxz<double>> r_tmp;
            T2 iChild, iParent = nodeParentRx;
            sxz<double> child;

            // store the son's coord
            child.x = Rx[n].x;
            child.z = Rx[n].z;
            while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {

                r_tmp.push_back( child );

                // we now go up in time - parent becomes the child of grand'pa
                iChild = iParent;
                child.x = (*node_p)[iChild].getX();
                child.z = (*node_p)[iChild].getZ();

                // grand'pa is now papa
                iParent = (*node_p)[iChild].getNodeParent(threadNo);
                if ( iParent >= this->nodes.size() ) {
                    node_p = &txNodes;
                    iParent -= this->nodes.size();
                }
                else {
                    node_p = &(this->nodes);
                }
            }

            // parent is now at Tx
            r_tmp.push_back( child );

            // finally, store Tx position
            child.x = (*node_p)[iParent].getX();
            child.z = (*node_p)[iParent].getZ();
            r_tmp.push_back( child );

            // the order should be from Tx to Rx, so we reorder...
            iParent = static_cast<T2>(r_tmp.size());
            r_data[n].resize( r_tmp.size() );
            for ( size_t nn=0; nn<r_data[n].size(); ++nn ) {
                r_data[n][nn].x = r_tmp[ iParent-1-nn ].x;
                r_data[n][nn].z = r_tmp[ iParent-1-nn ].z;
            }
        }
    }

    template<typename T1, typename T2, typename S, typename CELL>
    void Grid2Drcsp<T1,T2,S,CELL>::raytrace(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const std::vector<const std::vector<S>*>& Rx,
                                            std::vector<std::vector<T1>*>& traveltimes,
                                            std::vector<std::vector<std::vector<S>>*>& r_data,
                                            const size_t threadNo) const {

        this->checkPts(Tx);
        for ( size_t n=0; n<Rx.size(); ++n ) {
            this->checkPts(*Rx[n]);
        }

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node2Dcsp<T1,T2>*, std::vector<Node2Dcsp<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );

        std::vector<Node2Dcsp<T1,T2>> txNodes;
        std::vector<bool> inQueue( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );

        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);

        propagate(queue, inQueue, frozen, threadNo);

        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        if ( r_data.size() != Rx.size() ) {
            r_data.resize( Rx.size() );
        }

        for (size_t nr=0; nr<Rx.size(); ++nr) {

            traveltimes[nr]->resize( Rx[nr]->size() );
            r_data[nr]->resize( Rx[nr]->size() );
            for ( size_t ni=0; ni<r_data[nr]->size(); ++ni ) {
                (*r_data[nr])[ni].resize( 0 );
            }

            T2 nodeParentRx;
            T2 cellParentRx;

            for (size_t n=0; n<Rx[nr]->size(); ++n) {

                // check if Rx is on one of Tx nodes
                bool foundTx = false;
                for (size_t nt=0; nt<Tx.size(); ++nt) {
                    if (Tx[nt] == (*Rx[nr])[n]) {
                        (*traveltimes[nr])[n] = t0[nt];
                        (*r_data[nr])[n].push_back(Tx[nt]);
                        foundTx = true;
                        break;
                    }
                }
                if (foundTx)
                    continue;

                (*traveltimes[nr])[n] = getTraveltime((*Rx[nr])[n], this->nodes,
                                                      nodeParentRx, cellParentRx,
                                                      threadNo);

                bool flag=false;
                for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                    if ( (*Rx[nr])[n] == Tx[ns] ) {

                        (*r_data[nr])[n].resize( 1 );
                        (*r_data[nr])[n][0].x = (*Rx[nr])[n].x;
                        (*r_data[nr])[n][0].z = (*Rx[nr])[n].z;

                        flag = true;
                        break;
                    }
                }
                if ( flag ) continue;

                // Rx are in nodes (not txNodes)
                std::vector<Node2Dcsp<T1,T2>> *node_p;
                node_p = &(this->nodes);

                std::vector<S> r_tmp;
                T2 iChild, iParent = nodeParentRx;
                S child;

                // store the son's coord
                child.x = (*Rx[nr])[n].x;
                child.z = (*Rx[nr])[n].z;
                while ( (*node_p)[iParent].getNodeParent(threadNo) !=
                       std::numeric_limits<T2>::max() ) {

                    r_tmp.push_back( child );

                    // we now go up in time - parent becomes the child of grand'pa
                    iChild = iParent;
                    child.x = (*node_p)[iChild].getX();
                    child.z = (*node_p)[iChild].getZ();

                    // grand'pa is now papa
                    iParent = (*node_p)[iChild].getNodeParent(threadNo);
                    if ( iParent >= this->nodes.size() ) {
                        node_p = &txNodes;
                        iParent -= this->nodes.size();
                    }
                    else {
                        node_p = &(this->nodes);
                    }
                }

                // parent is now at Tx
                r_tmp.push_back( child );

                // finally, store Tx position
                child.x = (*node_p)[iParent].getX();
                child.z = (*node_p)[iParent].getZ();
                r_tmp.push_back( child );

                // the order should be from Tx to Rx, so we reorder...
                iParent = static_cast<T2>(r_tmp.size());
                (*r_data[nr])[n].resize( r_tmp.size() );
                for ( size_t nn=0; nn<(*r_data[nr])[n].size(); ++nn ) {
                    (*r_data[nr])[n][nn].x = r_tmp[ iParent-1-nn ].x;
                    (*r_data[nr])[n][nn].z = r_tmp[ iParent-1-nn ].z;
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename CELL>
    template<typename SIV>
    void Grid2Drcsp<T1,T2,S,CELL>::raytraceSensitivity(const std::vector<S>& Tx,
                                                       const std::vector<T1>& t0,
                                                       const std::vector<S>& Rx,
                                                       std::vector<T1>& traveltimes,
                                                       std::vector<std::vector<S>>* r_data,
                                                       std::vector<std::vector<SIV>>& l_data,
                                                       const size_t threadNo) const {

        if constexpr ( !cell_reports_into<CELL, SIV, Node2Dcsp<T1,T2>, S>::value ) {
            throw std::logic_error("Error: these cells do not report their "
                                   "sensitivity in the container requested.");
        } else {

        this->checkPts(Tx);
        this->checkPts(Rx);

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node2Dcsp<T1,T2>*, std::vector<Node2Dcsp<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );
        std::vector<Node2Dcsp<T1,T2>> txNodes;
        std::vector<bool> inQueue( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );

        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);

        propagate(queue, inQueue, frozen, threadNo);

        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        if ( l_data.size() != Rx.size() ) {
            l_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<l_data.size(); ++ni ) {
            l_data[ni].resize( 0 );
        }
        if ( r_data != nullptr ) {
            if ( r_data->size() != Rx.size() ) {
                r_data->resize( Rx.size() );
            }
            for ( size_t ni=0; ni<r_data->size(); ++ni ) {
                (*r_data)[ni].resize( 0 );
            }
        }
        T2 nodeParentRx;
        T2 cellParentRx;

        for (size_t n=0; n<Rx.size(); ++n) {

            // check if Rx is on one of Tx nodes
            bool foundTx = false;
            for (size_t nt=0; nt<Tx.size(); ++nt) {
                if (Tx[nt] == Rx[n]) {
                    traveltimes[n] = t0[nt];
                    if ( r_data != nullptr ) {
                        (*r_data)[n].push_back(Tx[nt]);
                    }
                    foundTx = true;
                    break;
                }
            }
            if (foundTx)
                continue;

            traveltimes[n] = getTraveltime(Rx[n], this->nodes, nodeParentRx,
                                           cellParentRx, threadNo);

            // Rx are in nodes (not txNodes)
            std::vector<Node2Dcsp<T1,T2>> *node_p;
            node_p = &(this->nodes);

            std::vector<S> r_tmp;
            T2 iChild, iParent = nodeParentRx;
            S child;
            SIV cell;

            // store the son's coord
            child.x = Rx[n].x;
            child.z = Rx[n].z;
            cell.i = cellParentRx;
            while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {

                if ( r_data != nullptr ) r_tmp.push_back( child );

                this->cells.computeDistance( (*node_p)[iParent], child, cell);
                accumulate(l_data[n], cell);

                // we now go up in time - parent becomes the child of grand'pa
                iChild = iParent;
                child.x = (*node_p)[iChild].getX();
                child.z = (*node_p)[iChild].getZ();
                cell.i = (*node_p)[iChild].getCellParent(threadNo);

                // grand'pa is now papa
                iParent = (*node_p)[iChild].getNodeParent(threadNo);
                if ( iParent >= this->nodes.size() ) {
                    node_p = &txNodes;
                    iParent -= this->nodes.size();
                }
                else {
                    node_p = &(this->nodes);
                }
            }

            // parent is now at Tx
            if ( r_data != nullptr ) r_tmp.push_back( child );

            this->cells.computeDistance( (*node_p)[iParent], child, cell);
            accumulate(l_data[n], cell);

            // finally, store Tx position
            child.x = (*node_p)[iParent].getX();
            child.z = (*node_p)[iParent].getZ();
            if ( r_data != nullptr ) r_tmp.push_back( child );

            //  must be sorted to build matrix L
            std::sort(l_data[n].begin(), l_data[n].end(),
                      [](const SIV& a, const SIV& b) { return a.i < b.i; });

            if ( r_data != nullptr ) {
                // the order should be from Tx to Rx, so we reorder...
                iParent = static_cast<T2>(r_tmp.size());
                (*r_data)[n].resize( r_tmp.size() );
                for ( size_t nn=0; nn<(*r_data)[n].size(); ++nn ) {
                    (*r_data)[n][nn].x = r_tmp[ iParent-1-nn ].x;
                    (*r_data)[n][nn].z = r_tmp[ iParent-1-nn ].z;
                }
            }
        }
        }
    }


    template<typename T1, typename T2, typename S, typename CELL>
    void Grid2Drcsp<T1,T2,S,CELL>::propagate(std::priority_queue<Node2Dcsp<T1,T2>*,
                                             std::vector<Node2Dcsp<T1,T2>*>,
                                             CompareNodePtr<T1>>& queue,
                                             std::vector<bool>& inQueue,
                                             std::vector<bool>& frozen,
                                             const size_t threadNo) const {

        while ( !queue.empty() ) {
            const Node2Dcsp<T1,T2>* source = queue.top();
            queue.pop();
            inQueue[ source->getGridIndex() ] = false;
            frozen[ source->getGridIndex() ] = true;

            for ( size_t no=0; no<source->getOwners().size(); ++no ) {

                T2 cellNo = source->getOwners()[no];

                for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                    T2 neibNo = this->neighbors[cellNo][k];
                    if ( neibNo == source->getGridIndex() || frozen[neibNo] ) {
                        continue;
                    }

                    // compute dt
                    T1 dt = this->cells.computeDt(*source, this->nodes[neibNo], cellNo);

                    if ( source->getTT(threadNo)+dt < this->nodes[neibNo].getTT(threadNo) ) {
                        this->nodes[neibNo].setTT( source->getTT(threadNo)+dt, threadNo );
                        this->nodes[neibNo].setnodeParent( source->getGridIndex(), threadNo );
                        this->nodes[neibNo].setCellParent( cellNo, threadNo );

                        if ( !inQueue[neibNo] ) {
                            queue.push( &(this->nodes[neibNo]) );
                            inQueue[neibNo] = true;
                        }
                    }
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename CELL>
    void Grid2Drcsp<T1,T2,S,CELL>::propagate_lw(std::priority_queue<Node2Dcsp<T1,T2>*,
                                                std::vector<Node2Dcsp<T1,T2>*>,
                                                CompareNodePtr<T1>>& queue,
                                                std::vector<bool>& inQueue,
                                                std::vector<bool>& frozen,
                                                const size_t threadNo) const {
        // lightweight method where cell/node parent are not stored
        while ( !queue.empty() ) {
            const Node2Dcsp<T1,T2>* source = queue.top();
            queue.pop();
            inQueue[ source->getGridIndex() ] = false;
            frozen[ source->getGridIndex() ] = true;

            for ( size_t no=0; no<source->getOwners().size(); ++no ) {
                T2 cellNo = source->getOwners()[no];
                for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                    T2 neibNo = this->neighbors[cellNo][k];
                    if ( neibNo == source->getGridIndex() || frozen[neibNo] ) {
                        continue;
                    }

                    // compute dt
                    T1 dt = this->cells.computeDt(*source, this->nodes[neibNo], cellNo);

                    if ( source->getTT(threadNo)+dt < this->nodes[neibNo].getTT(threadNo) ) {
                        this->nodes[neibNo].setTT( source->getTT(threadNo)+dt, threadNo );
                        if ( !inQueue[neibNo] ) {
                            queue.push( &(this->nodes[neibNo]) );
                            inQueue[neibNo] = true;
                        }
                    }
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename CELL>
    T1 Grid2Drcsp<T1,T2,S,CELL>::getTraveltime(const S& Rx,
                                               const std::vector<Node2Dcsp<T1,T2>>& nodes,
                                               T2& nodeParentRx, T2& cellParentRx,
                                               const size_t threadNo) const {

        for ( size_t nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == Rx ) {
                nodeParentRx = nodes[nn].getNodeParent(threadNo);
                cellParentRx = nodes[nn].getCellParent(threadNo);
                return nodes[nn].getTT(threadNo);
            }
        }

        T2 cellNo = this->getCellNo( Rx );
        T2 neibNo = this->neighbors[cellNo][0];
        T1 dt = this->cells.computeDt(nodes[neibNo], Rx, cellNo);

        T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
        nodeParentRx = neibNo;
        cellParentRx = cellNo;
        for ( size_t k=1; k< this->neighbors[cellNo].size(); ++k ) {
            neibNo = this->neighbors[cellNo][k];
            dt = this->cells.computeDt(nodes[neibNo], Rx, cellNo);
            if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
                traveltime =  nodes[neibNo].getTT(threadNo)+dt;
                nodeParentRx = neibNo;
            }
        }
        return traveltime;
    }

}

#endif
