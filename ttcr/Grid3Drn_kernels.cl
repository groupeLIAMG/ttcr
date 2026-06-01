/*
 * OpenCL Kernels for Grid3Drn Fast Marching Sweeps
 * 
 * This file contains optimized OpenCL kernels for:
 * 1. Basic sweep update (first-order)
 * 2. WENO3 sweep update (third-order)
 *
 * Usage:
 *   - Launch one kernel per sweep direction (8 total)
 *   - Use 3D work groups for optimal memory access
 *   - Synchronize between sweep directions
 */

// =============================================================================
// TYPE DEFINITIONS
// =============================================================================

// Support both float and double precision.
// The host (OpenCLSweepSolver<T1>::buildKernels) passes -DUSE_DOUBLE=1 when the
// grid's T1 is double, so real_t always matches the host scalar type T1.  In
// the float build the host instead passes -cl-single-precision-constant so the
// unsuffixed double literals below are computed in single precision (no
// float<->double promotion); in the double build they stay double.  Therefore
// floating-point literals here must be written WITHOUT an 'f' suffix -- a hard
// 'f' would pin them to float even in the double build.
//
// WENO_EPS is the WENO smoothness regularizer; it is scaled to the working
// precision (matching the CPU which uses numeric_limits<T1>::epsilon()).
#ifdef USE_DOUBLE
    #if USE_DOUBLE
        #pragma OPENCL EXTENSION cl_khr_fp64 : enable
        typedef double real_t;
        #define REAL_MAX DBL_MAX
        #define WENO_EPS 1.0e-15
    #else
        typedef float real_t;
        #define REAL_MAX FLT_MAX
        #define WENO_EPS 1.2e-7
    #endif
#else
    // Default to float for backward compatibility
    typedef float real_t;
    #define REAL_MAX FLT_MAX
    #define WENO_EPS 1.2e-7
#endif

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

/**
 * Inline device function to compute 1D index from 3D coordinates
 */
inline size_t get_index_3d(size_t i, size_t j, size_t k, 
                           size_t nx, size_t ny) {
    return (k * (ny + 1) + j) * (nx + 1) + i;
}

/**
 * Inline device function for min of two values
 */
inline real_t fmin2(real_t a, real_t b) {
    return (a < b) ? a : b;
}

/**
 * Inline device function to swap two values
 */
inline void swap(real_t *a, real_t *b) {
    real_t temp = *a;
    *a = *b;
    *b = temp;
}

/**
 * Sort three values in ascending order
 */
inline void sort3(real_t *a1, real_t *a2, real_t *a3) {
    if (*a1 > *a2) swap(a1, a2);
    if (*a1 > *a3) swap(a1, a3);
    if (*a2 > *a3) swap(a2, a3);
}

// =============================================================================
// KERNEL 1: BASIC SWEEP (First-Order)
// =============================================================================

/**
 * Basic sweep kernel - corresponds to Grid3Drn::sweep()
 * 
 * This kernel updates travel times using a first-order finite difference
 * scheme. It must be called 8 times (once per sweep direction).
 * 
 * @param tt_in      Input travel time array (read-only)
 * @param tt_out     Output travel time array (write)
 * @param slowness   Slowness values at grid nodes
 * @param frozen     Boolean array indicating frozen nodes
 * @param dx         Grid spacing in x
 * @param dy         Grid spacing in y  
 * @param dz         Grid spacing in z
 * @param ncx        Number of cells in x (nodes = ncx+1)
 * @param ncy        Number of cells in y (nodes = ncy+1)
 * @param ncz        Number of cells in z (nodes = ncz+1)
 * @param i_start    Starting index in i (for sweep direction)
 * @param j_start    Starting index in j
 * @param k_start    Starting index in k
 * @param i_dir      Direction in i (-1 or +1)
 * @param j_dir      Direction in j (-1 or +1)
 * @param k_dir      Direction in k (-1 or +1)
 */
__kernel void sweep_update_basic(
    __global const real_t *tt_in,
    __global real_t *tt_out,
    __global const real_t *slowness,
    __global const uchar *frozen,
    const real_t dx,
    const real_t dy,
    const real_t dz,
    const uint ncx,
    const uint ncy,
    const uint ncz,
    // OLD CODE (3-D full-box dispatch + per-plane gating):
    //
    //     const int i_start,
    //     const int j_start,
    //     const int k_start,
    //     const int i_dir,
    //     const int j_dir,
    //     const int k_dir,
    //     const int level)
    //
    // Each plane launch dispatched the entire (ncx+1)x(ncy+1)x(ncz+1) box and
    // every off-plane work-item returned at the gate below -- O(n) launches
    // each spawning O(n^3) threads for only O(n^2) useful updates.  The 6
    // direction arguments existed solely to drive that gate; the eikonal solver
    // never used them.  They are replaced by an explicit per-(direction,level)
    // node list: the host precomputes, for each sweep direction, the node
    // indices grouped by level, and each launch dispatches a 1-D range over
    // exactly that plane's nodes (see buildPlaneNodeLists / executeSweep in the
    // host code).
    __global const uint *plane_nodes,
    const uint plane_offset,
    const uint plane_count)
{
    // One work-item per node on the current plane.
    const size_t gid = get_global_id(0);
    if (gid >= plane_count) return;

    // The node's linear index is read straight from the precomputed plane list;
    // decode it back into (i, j, k) for the stencil below.
    const size_t idx = plane_nodes[plane_offset + gid];
    const size_t i = idx % (ncx + 1);
    const size_t j = (idx / (ncx + 1)) % (ncy + 1);
    const size_t k = idx / ((ncx + 1) * (ncy + 1));

    // OLD CODE (3-D thread index + bounds check + Gauss-Seidel plane gating):
    //     const size_t i = get_global_id(0);
    //     const size_t j = get_global_id(1);
    //     const size_t k = get_global_id(2);
    //     if (i > ncx || j > ncy || k > ncz) return;
    //     {
    //         const int ip = (i_dir > 0) ? (int)i : (int)(ncx - i);
    //         const int jp = (j_dir > 0) ? (int)j : (int)(ncy - j);
    //         const int kp = (k_dir > 0) ? (int)k : (int)(ncz - k);
    //         if (ip + jp + kp != level) return;
    //     }
    //     const size_t idx = get_index_3d(i, j, k, ncx, ncy);

    // Skip if frozen
    if (frozen[idx]) {
        tt_out[idx] = tt_in[idx];
        return;
    }

    // Get neighbor values in each direction
    real_t a1, a2, a3, t;
    
    // Z-direction neighbors
    if (k == 0) {
        a1 = tt_in[get_index_3d(i, j, k+1, ncx, ncy)];
    } else if (k == ncz) {
        a1 = tt_in[get_index_3d(i, j, k-1, ncx, ncy)];
    } else {
        a1 = tt_in[get_index_3d(i, j, k-1, ncx, ncy)];
        t  = tt_in[get_index_3d(i, j, k+1, ncx, ncy)];
        a1 = fmin2(a1, t);
    }
    
    // Y-direction neighbors
    if (j == 0) {
        a2 = tt_in[get_index_3d(i, j+1, k, ncx, ncy)];
    } else if (j == ncy) {
        a2 = tt_in[get_index_3d(i, j-1, k, ncx, ncy)];
    } else {
        a2 = tt_in[get_index_3d(i, j-1, k, ncx, ncy)];
        t  = tt_in[get_index_3d(i, j+1, k, ncx, ncy)];
        a2 = fmin2(a2, t);
    }
    
    // X-direction neighbors
    if (i == 0) {
        a3 = tt_in[get_index_3d(i+1, j, k, ncx, ncy)];
    } else if (i == ncx) {
        a3 = tt_in[get_index_3d(i-1, j, k, ncx, ncy)];
    } else {
        a3 = tt_in[get_index_3d(i-1, j, k, ncx, ncy)];
        t  = tt_in[get_index_3d(i+1, j, k, ncx, ncy)];
        a3 = fmin2(a3, t);
    }
    
    // Sort values
    sort3(&a1, &a2, &a3);
    
    // Compute slowness * grid_spacing
    const real_t fh = slowness[idx] * dx;  // Assuming dx == dy == dz
    
    // Solve eikonal equation
    t = a1 + fh;
    
    if (t > a2) {
        // 2D update
        t = 0.5 * (a1 + a2 + sqrt(2.0 * fh * fh - (a1 - a2) * (a1 - a2)));
        
        if (t > a3) {
            // 3D update
            const real_t sum_a = a1 + a2 + a3;
            const real_t term = -2.0 * a1 * a1 + 2.0 * a1 * a2 - 2.0 * a2 * a2 +
                               2.0 * a1 * a3 + 2.0 * a2 * a3 - 2.0 * a3 * a3 + 
                               3.0 * fh * fh;
            
            t = (sum_a + sqrt(term)) / 3.0;
        }
    }
    
    // Update only if smaller
    tt_out[idx] = fmin2(t, tt_in[idx]);
}

// =============================================================================
// KERNEL 2: WENO3 SWEEP (Third-Order)
// =============================================================================

inline real_t weno3_upwind(real_t v0, real_t v1, real_t v2, real_t v3, real_t v4, real_t dx, bool forward)
{
    const real_t eps = WENO_EPS;

    if (forward) {
        // Forward differencing: ap = d/dx approximation
        const real_t num = (v4 - 2.0 * v3 + v2);
        const real_t den = (v3 - 2.0 * v2 + v1);
        const real_t r = (eps + num * num) / (eps + den * den);
        const real_t w = 1.0 / (1.0 + 2.0 * r * r);

        const real_t ap = (1.0 - w) * (v3 - v1) / (2.0 * dx) +
                         w * (-v4 + 4.0 * v3 - 3.0 * v2) / (2.0 * dx);

        return v2 + dx * ap;
    } else {
        // Backward differencing: am = -d/dx approximation
        const real_t num = (v2 - 2.0 * v1 + v0);
        const real_t den = (v3 - 2.0 * v2 + v1);
        const real_t r = (eps + num * num) / (eps + den * den);
        const real_t w = 1.0 / (1.0 + 2.0 * r * r);

        const real_t am = (1.0 - w) * (v3 - v1) / (2.0 * dx) +
                         w * (3.0 * v2 - 4.0 * v1 + v0) / (2.0 * dx);

        return v2 - dx * am;
    }
}

/**
 * WENO3 sweep kernel - corresponds to Grid3Drn::sweep_weno3()
 *
 * This kernel uses third-order WENO interpolation for more accurate
 * travel time updates. Requires more computation but better accuracy.
 *
 * Parameters same as sweep_update_basic
 */
__kernel void sweep_update_weno3(
    __global const real_t *tt_in,
    __global real_t *tt_out,
    __global const real_t *slowness,
    __global const uchar *frozen,
    const real_t dx,
    const real_t dy,
    const real_t dz,
    const uint ncx,
    const uint ncy,
    const uint ncz,
    // OLD CODE (3-D full-box dispatch + per-plane gating): see the matching
    // comment in sweep_update_basic.  The 6 direction arguments + level
    //
    //     const int i_start,
    //     const int j_start,
    //     const int k_start,
    //     const int i_dir,
    //     const int j_dir,
    //     const int k_dir,
    //     const int level)
    //
    // drove a per-plane gate; they are replaced by an explicit
    // per-(direction,level) node list dispatched as a 1-D range.
    __global const uint *plane_nodes,
    const uint plane_offset,
    const uint plane_count)
{
    // One work-item per node on the current plane.
    const size_t gid = get_global_id(0);
    if (gid >= plane_count) return;

    // The node's linear index is read straight from the precomputed plane list;
    // decode it back into (i, j, k) for the WENO stencil below.
    const size_t idx = plane_nodes[plane_offset + gid];
    const size_t i = idx % (ncx + 1);
    const size_t j = (idx / (ncx + 1)) % (ncy + 1);
    const size_t k = idx / ((ncx + 1) * (ncy + 1));

    // OLD CODE (3-D thread index + bounds check + Gauss-Seidel plane gating):
    //     const size_t i = get_global_id(0);
    //     const size_t j = get_global_id(1);
    //     const size_t k = get_global_id(2);
    //     if (i > ncx || j > ncy || k > ncz) return;
    //     {
    //         const int ip = (i_dir > 0) ? (int)i : (int)(ncx - i);
    //         const int jp = (j_dir > 0) ? (int)j : (int)(ncy - j);
    //         const int kp = (k_dir > 0) ? (int)k : (int)(ncz - k);
    //         if (ip + jp + kp != level) return;
    //     }
    //     const size_t idx = get_index_3d(i, j, k, ncx, ncy);

    // Skip if frozen
    if (frozen[idx]) {
        tt_out[idx] = tt_in[idx];
        return;
    }

//    const real_t eps = 1.2e-7;
    real_t a1, a2, a3, t;
    
    // =========================================================================
    // Z-DIRECTION (k) - WENO3 approximation WITH BOUNDS CHECKING
    // =========================================================================
    
    if (k == 0) {
        // Boundary: first order
        a1 = tt_in[get_index_3d(i, j, k+1, ncx, ncy)];
    } else if (k == ncz) {
        a1 = tt_in[get_index_3d(i, j, k-1, ncx, ncy)];
    } else if (k == 1) {
        // Near boundary: first-order + neighbor
        const real_t v0 = 0.0;
        const real_t v1 = tt_in[get_index_3d(i, j, k-1, ncx, ncy)];
        const real_t v2 = tt_in[get_index_3d(i, j, k,   ncx, ncy)];
        const real_t v3 = tt_in[get_index_3d(i, j, k+1, ncx, ncy)];
        const real_t v4 = tt_in[get_index_3d(i, j, k+2, ncx, ncy)];
        
        a1 = weno3_upwind(v0, v1, v2, v3, v4, dx, true);
        a1 = fmin2(a1, v1);
//        real_t num = tt_in[get_index_3d(i, j, k+2, ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i, j, k+1, ncx, ncy)] +
//                     tt_in[get_index_3d(i, j, k  , ncx, ncy)];
//        num *= num;
//        real_t den = tt_in[get_index_3d(i, j, k+1, ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i, j, k  , ncx, ncy)] +
//                     tt_in[get_index_3d(i, j, k-1, ncx, ncy)];
//        den *= den;
//        real_t r = (eps+num)/(eps+den);
//        real_t w = 1./(1.+2.*r*r);
//
//        real_t ap = (1.-w)*(tt_in[get_index_3d(i, j, k+1, ncx, ncy)]-
//                            tt_in[get_index_3d(i, j, k-1, ncx, ncy)])/(2.*dx) +
//        w*(  -tt_in[get_index_3d(i, j, k+2, ncx, ncy)] +
//           4.*tt_in[get_index_3d(i, j, k+1, ncx, ncy)] -
//           3.*tt_in[get_index_3d(i, j, k  , ncx, ncy)])/(2.*dx);
//
//        a1 = tt_in[get_index_3d(i, j, k, ncx, ncy)] + dx*ap;
//        
//        t = tt_in[get_index_3d(i, j, k-1, ncx, ncy)]; // first order for left
//        a1 = a1<t ? a1 : t;
    } else if (k == ncz - 1) {
        // Near boundary: first-order + neighbor
        const real_t v0 = tt_in[get_index_3d(i, j, k-2, ncx, ncy)];
        const real_t v1 = tt_in[get_index_3d(i, j, k-1, ncx, ncy)];
        const real_t v2 = tt_in[get_index_3d(i, j, k  , ncx, ncy)];
        const real_t v3 = tt_in[get_index_3d(i, j, k+1, ncx, ncy)];
        const real_t v4 = 0.0;
        
        a1 = weno3_upwind(v0, v1, v2, v3, v4, dx, false);
        a1 = fmin2(a1, v3);
//        real_t num = tt_in[get_index_3d(i, j, k  , ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i, j, k-1, ncx, ncy)] +
//                     tt_in[get_index_3d(i, j, k-2, ncx, ncy)];
//        num *= num;
//        real_t den = tt_in[get_index_3d(i, j, k+1, ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i, j, k  , ncx, ncy)] +
//                     tt_in[get_index_3d(i, j, k-1, ncx, ncy)];
//        den *= den;
//        real_t r = (eps+num)/(eps+den);
//        real_t w = 1./(1.+2.*r*r);
//
//        real_t am = (1.-w)*(tt_in[get_index_3d(i, j, k+1, ncx, ncy)]-
//                            tt_in[get_index_3d(i, j, k-1, ncx, ncy)])/(2.*dx) +
//        w*(3.*tt_in[get_index_3d(i, j, k  , ncx, ncy)] -
//           4.*tt_in[get_index_3d(i, j, k-1, ncx, ncy)] +
//              tt_in[get_index_3d(i, j, k-2, ncx, ncy)])/(2.*dx);
//
//        a1 = tt_in[get_index_3d(i, j, k  , ncx, ncy)] - dx*am;
//
//        t = tt_in[get_index_3d(i, j, k+1, ncx, ncy)]; // first order for right
//        a1 = a1<t ? a1 : t;
    } else {
                // Interior: WENO3 from both directions, take minimum
                const real_t v0 = tt_in[get_index_3d(i, j, k-2, ncx, ncy)];
                const real_t v1 = tt_in[get_index_3d(i, j, k-1, ncx, ncy)];
                const real_t v2 = tt_in[get_index_3d(i, j, k  , ncx, ncy)];
                const real_t v3 = tt_in[get_index_3d(i, j, k+1, ncx, ncy)];
                const real_t v4 = tt_in[get_index_3d(i, j, k+2, ncx, ncy)];
        
                a1 = weno3_upwind(v0, v1, v2, v3, v4, dx, true);
                t = weno3_upwind(v0, v1, v2, v3, v4, dx, false);
                a1 = fmin2(a1, t);
//        real_t num = tt_in[get_index_3d(i, j, k+2, ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i, j, k+1, ncx, ncy)] +
//                     tt_in[get_index_3d(i, j, k  , ncx, ncy)];
//        num *= num;
//        real_t den = tt_in[get_index_3d(i, j, k+1, ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i, j, k  , ncx, ncy)] +
//                     tt_in[get_index_3d(i, j, k-1, ncx, ncy)];
//        den *= den;
//        real_t r = (eps+num)/(eps+den);
//        real_t w = 1./(1.+2.*r*r);
//        
//        real_t ap = (1.-w)*(tt_in[get_index_3d(i, j, k+1, ncx, ncy)] -
//                            tt_in[get_index_3d(i, j, k-1, ncx, ncy)])/(2.*dx) +
//        w*(  -tt_in[get_index_3d(i, j, k+2, ncx, ncy)] +
//           4.*tt_in[get_index_3d(i, j, k+1, ncx, ncy)] -
//           3.*tt_in[get_index_3d(i, j, k  , ncx, ncy)])/(2.*dx);
//        
//        a1 = tt_in[get_index_3d(i, j, k  , ncx, ncy)] + dx*ap;
//        
//        num = tt_in[get_index_3d(i, j, k  , ncx, ncy)] -
//           2.*tt_in[get_index_3d(i, j, k-1, ncx, ncy)] +
//              tt_in[get_index_3d(i, j, k-2, ncx, ncy)];
//        num *= num;
//        r = (eps+num)/(eps+den);
//        w = 1./(1.+2.*r*r);
//        
//        real_t am = (1.-w)*(tt_in[get_index_3d(i, j, k+1, ncx, ncy)] -
//                            tt_in[get_index_3d(i, j, k-1, ncx, ncy)])/(2.*dx) +
//        w*(3.*tt_in[get_index_3d(i, j, k  , ncx, ncy)] -
//           4.*tt_in[get_index_3d(i, j, k-1, ncx, ncy)] +
//              tt_in[get_index_3d(i, j, k-2, ncx, ncy)])/(2.*dx);
//        
//        t = tt_in[get_index_3d(i, j, k  , ncx, ncy)] - dx*am;
//        
//        a1 = a1<t ? a1 : t;
    }
    
    // =========================================================================
    // Y-DIRECTION (j) - WENO3 approximation WITH BOUNDS CHECKING
    // =========================================================================
    
    if (j == 0) {
        a2 = tt_in[get_index_3d(i, j+1, k, ncx, ncy)];
    } else if (j == ncy) {
        a2 = tt_in[get_index_3d(i, j-1, k, ncx, ncy)];
    } else if (j == 1) {
        const real_t v0 = 0.0;
        const real_t v1 = tt_in[get_index_3d(i, j-1, k, ncx, ncy)];
        const real_t v2 = tt_in[get_index_3d(i, j  , k, ncx, ncy)];
        const real_t v3 = tt_in[get_index_3d(i, j+1, k, ncx, ncy)];
        const real_t v4 = tt_in[get_index_3d(i, j+2, k, ncx, ncy)];
        
        a2 = weno3_upwind(v0, v1, v2, v3, v4, dx, true);
        a2 = fmin2(a2, v1);
//        real_t num = tt_in[get_index_3d(i, j+2, k, ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i, j+1, k, ncx, ncy)] +
//                     tt_in[get_index_3d(i, j  , k, ncx, ncy)];
//        num *= num;
//        real_t den = tt_in[get_index_3d(i, j+1, k, ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i, j  , k, ncx, ncy)] +
//                     tt_in[get_index_3d(i, j-1, k, ncx, ncy)];
//        den *= den;
//        real_t r = (eps+num)/(eps+den);
//        real_t w = 1./(1.+2.*r*r);
//
//        real_t ap = (1.-w)*(tt_in[get_index_3d(i, j+1, k, ncx, ncy)]-
//                            tt_in[get_index_3d(i, j-1, k, ncx, ncy)])/(2.*dx) +
//        w*(  -tt_in[get_index_3d(i, j+2, k, ncx, ncy)] +
//           4.*tt_in[get_index_3d(i, j+1, k, ncx, ncy)] -
//           3.*tt_in[get_index_3d(i, j  , k, ncx, ncy)])/(2.*dx);
//
//        a2 = tt_in[get_index_3d(i, j  , k, ncx, ncy)] + dx*ap;
//
//        t = tt_in[get_index_3d(i, j-1, k, ncx, ncy)]; // first order for left
//        a2 = a2<t ? a2 : t;
    } else if (j == ncy - 1) {
        const real_t v0 = tt_in[get_index_3d(i, j-2, k, ncx, ncy)];
        const real_t v1 = tt_in[get_index_3d(i, j-1, k, ncx, ncy)];
        const real_t v2 = tt_in[get_index_3d(i, j  , k, ncx, ncy)];
        const real_t v3 = tt_in[get_index_3d(i, j+1, k, ncx, ncy)];
        const real_t v4 = 0.0;
        
        a2 = weno3_upwind(v0, v1, v2, v3, v4, dx, false);
        a2 = fmin2(a2, v3);
//        real_t num = tt_in[get_index_3d(i, j  , k, ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i, j-1, k, ncx, ncy)] +
//                     tt_in[get_index_3d(i, j-2, k, ncx, ncy)];
//        num *= num;
//        real_t den = tt_in[get_index_3d(i, j+1, k, ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i, j  , k, ncx, ncy)] +
//                     tt_in[get_index_3d(i, j-1, k, ncx, ncy)];
//        den *= den;
//        real_t r = (eps+num)/(eps+den);
//        real_t w = 1./(1.+2.*r*r);
//
//        real_t am = (1.-w)*(tt_in[get_index_3d(i, j+1, k, ncx, ncy)]-
//                            tt_in[get_index_3d(i, j-1, k, ncx, ncy)])/(2.*dx) +
//        w*(3.*tt_in[get_index_3d(i, j  , k, ncx, ncy)] -
//           4.*tt_in[get_index_3d(i, j-1, k, ncx, ncy)] +
//              tt_in[get_index_3d(i, j-2, k, ncx, ncy)])/(2.*dx);
//
//        a2 = tt_in[get_index_3d(i, j  , k, ncx, ncy)] - dx*am;
//
//        t = tt_in[get_index_3d(i, j+1, k, ncx, ncy)]; // first order for right
//        a2 = a2<t ? a2 : t;
    } else {
        const real_t v0 = tt_in[get_index_3d(i, j-2, k, ncx, ncy)];
        const real_t v1 = tt_in[get_index_3d(i, j-1, k, ncx, ncy)];
        const real_t v2 = tt_in[get_index_3d(i, j  , k, ncx, ncy)];
        const real_t v3 = tt_in[get_index_3d(i, j+1, k, ncx, ncy)];
        const real_t v4 = tt_in[get_index_3d(i, j+2, k, ncx, ncy)];
        
        a2 = weno3_upwind(v0, v1, v2, v3, v4, dx, true);
        t = weno3_upwind(v0, v1, v2, v3, v4, dx, false);
        a2 = fmin2(a2, t);
//        real_t num = tt_in[get_index_3d(i, j+2, k, ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i, j+1, k, ncx, ncy)] +
//                     tt_in[get_index_3d(i, j  , k, ncx, ncy)];
//        num *= num;
//        real_t den = tt_in[get_index_3d(i, j+1, k, ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i, j  , k, ncx, ncy)] +
//                     tt_in[get_index_3d(i, j-1, k, ncx, ncy)];
//        den *= den;
//        real_t r = (eps+num)/(eps+den);
//        real_t w = 1./(1.+2.*r*r);
//
//        real_t ap = (1.-w)*(tt_in[get_index_3d(i, j+1, k, ncx, ncy)]-
//                            tt_in[get_index_3d(i, j-1, k, ncx, ncy)])/(2.*dx) +
//        w*(  -tt_in[get_index_3d(i, j+2, k, ncx, ncy)] +
//           4.*tt_in[get_index_3d(i, j+1, k, ncx, ncy)] -
//           3.*tt_in[get_index_3d(i, j  , k, ncx, ncy)])/(2.*dx);
//
//        a2 = tt_in[get_index_3d(i, j  , k, ncx, ncy)] + dx*ap;
//
//        num = tt_in[get_index_3d(i, j  , k, ncx, ncy)] -
//           2.*tt_in[get_index_3d(i, j-1, k, ncx, ncy)] +
//              tt_in[get_index_3d(i, j-2, k, ncx, ncy)];
//        num *= num;
//        r = (eps+num)/(eps+den);
//        w = 1./(1.+2.*r*r);
//
//        real_t am = (1.-w)*(tt_in[get_index_3d(i, j+1, k, ncx, ncy)]-
//                            tt_in[get_index_3d(i, j-1, k, ncx, ncy)])/(2.*dx) +
//        w*(3.*tt_in[get_index_3d(i, j  , k, ncx, ncy)] -
//           4.*tt_in[get_index_3d(i, j-1, k, ncx, ncy)] +
//              tt_in[get_index_3d(i, j-2, k, ncx, ncy)])/(2.*dx);
//
//        t = tt_in[get_index_3d(i, j  , k, ncx, ncy)] - dx*am;
//
//        a2 = a2<t ? a2 : t;
    }
    
    // =========================================================================
    // X-DIRECTION (i) - WENO3 approximation WITH BOUNDS CHECKING
    // =========================================================================
    
    if (i == 0) {
        a3 = tt_in[get_index_3d(i+1, j, k, ncx, ncy)];
    } else if (i == ncx) {
        a3 = tt_in[get_index_3d(i-1, j, k, ncx, ncy)];
    } else if (i == 1) {
        const real_t v0 = 0.0;
        const real_t v1 = tt_in[get_index_3d(i-1, j, k, ncx, ncy)];
        const real_t v2 = tt_in[get_index_3d(i  , j, k, ncx, ncy)];
        const real_t v3 = tt_in[get_index_3d(i+1, j, k, ncx, ncy)];
        const real_t v4 = tt_in[get_index_3d(i+2, j, k, ncx, ncy)];

        a3 = weno3_upwind(v0, v1, v2, v3, v4, dx, true);
        a3 = fmin2(a3, v1);
//        real_t num = tt_in[get_index_3d(i+2, j, k, ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i+1, j, k, ncx, ncy)] +
//                     tt_in[get_index_3d(i  , j, k, ncx, ncy)];
//        num *= num;
//        real_t den = tt_in[get_index_3d(i+1, j, k, ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i  , j, k, ncx, ncy)] +
//                     tt_in[get_index_3d(i-1, j, k, ncx, ncy)];
//        den *= den;
//        real_t r = (eps+num)/(eps+den);
//        real_t w = 1./(1.+2.*r*r);
//
//        real_t ap = (1.-w)*(tt_in[get_index_3d(i+1, j, k, ncx, ncy)]-
//                            tt_in[get_index_3d(i-1, j, k, ncx, ncy)])/(2.*dx) +
//        w*(  -tt_in[get_index_3d(i+2, j, k, ncx, ncy)] +
//           4.*tt_in[get_index_3d(i+1, j, k, ncx, ncy)] -
//           3.*tt_in[get_index_3d(i  , j, k, ncx, ncy)])/(2.*dx);
//
//        a3 = tt_in[get_index_3d(i  , j, k, ncx, ncy)] + dx*ap;
//
//        t = tt_in[get_index_3d(i-1, j, k, ncx, ncy)]; // first order for left
//        a3 = a3<t ? a3 : t;
    } else if (i == ncx - 1) {
        const real_t v0 = tt_in[get_index_3d(i-2, j, k, ncx, ncy)];
        const real_t v1 = tt_in[get_index_3d(i-1, j, k, ncx, ncy)];
        const real_t v2 = tt_in[get_index_3d(i  , j, k, ncx, ncy)];
        const real_t v3 = tt_in[get_index_3d(i+1, j, k, ncx, ncy)];
        const real_t v4 = 0.0;
        
        a3 = weno3_upwind(v0, v1, v2, v3, v4, dx, false);
        a3 = fmin2(a3, v3);
//        real_t num = tt_in[get_index_3d(i  , j, k, ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i-1, j, k, ncx, ncy)] +
//                     tt_in[get_index_3d(i-2, j, k, ncx, ncy)];
//        num *= num;
//        real_t den = tt_in[get_index_3d(i+1, j, k, ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i  , j, k, ncx, ncy)] +
//                     tt_in[get_index_3d(i-1, j, k, ncx, ncy)];
//        den *= den;
//        real_t r = (eps+num)/(eps+den);
//        real_t w = 1./(1.+2.*r*r);
//
//        real_t am = (1.-w)*(tt_in[get_index_3d(i+1, j, k, ncx, ncy)]-
//                            tt_in[get_index_3d(i-1, j, k, ncx, ncy)])/(2.*dx) +
//        w*(3.*tt_in[get_index_3d(i  , j, k, ncx, ncy)] -
//           4.*tt_in[get_index_3d(i-1, j, k, ncx, ncy)] +
//              tt_in[get_index_3d(i-2, j, k, ncx, ncy)])/(2.*dx);
//
//        a3 = tt_in[get_index_3d(i  , j, k, ncx, ncy)] - dx*am;
//
//        t = tt_in[get_index_3d(i+1, j, k, ncx, ncy)]; // first order for right
//        a3 = a3<t ? a3 : t;
    } else {
        const real_t v0 = tt_in[get_index_3d(i-2, j, k, ncx, ncy)];
        const real_t v1 = tt_in[get_index_3d(i-1, j, k, ncx, ncy)];
        const real_t v2 = tt_in[get_index_3d(i  , j, k, ncx, ncy)];
        const real_t v3 = tt_in[get_index_3d(i+1, j, k, ncx, ncy)];
        const real_t v4 = tt_in[get_index_3d(i+2, j, k, ncx, ncy)];

        a3 = weno3_upwind(v0, v1, v2, v3, v4, dx, true);
        t = weno3_upwind(v0, v1, v2, v3, v4, dx, false);
        a3 = fmin2(a3, t);
//        real_t num = tt_in[get_index_3d(i+2, j, k, ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i+1, j, k, ncx, ncy)] +
//                     tt_in[get_index_3d(i  , j, k, ncx, ncy)];
//        num *= num;
//        real_t den = tt_in[get_index_3d(i+1, j, k, ncx, ncy)] -
//                  2.*tt_in[get_index_3d(i  , j, k, ncx, ncy)] +
//                     tt_in[get_index_3d(i-1, j, k, ncx, ncy)];
//        den *= den;
//        real_t r = (eps+num)/(eps+den);
//        real_t w = 1./(1.+2.*r*r);
//
//        real_t ap = (1.-w)*(tt_in[get_index_3d(i+1, j, k, ncx, ncy)]-
//                            tt_in[get_index_3d(i-1, j, k, ncx, ncy)])/(2.*dx) +
//        w*(  -tt_in[get_index_3d(i+2, j, k, ncx, ncy)] +
//           4.*tt_in[get_index_3d(i+1, j, k, ncx, ncy)] -
//           3.*tt_in[get_index_3d(i  , j, k, ncx, ncy)])/(2.*dx);
//
//        a3 = tt_in[get_index_3d(i  , j, k, ncx, ncy)] + dx*ap;
//
//        num = tt_in[get_index_3d(i  , j, k, ncx, ncy)] -
//           2.*tt_in[get_index_3d(i-1, j, k, ncx, ncy)] +
//              tt_in[get_index_3d(i-2, j, k, ncx, ncy)];
//        num *= num;
//        r = (eps+num)/(eps+den);
//        w = 1./(1.+2.*r*r);
//
//        real_t am = (1.-w)*(tt_in[get_index_3d(i+1, j, k, ncx, ncy)]-
//                            tt_in[get_index_3d(i-1, j, k, ncx, ncy)])/(2.*dx) +
//        w*(3.*tt_in[get_index_3d(i  , j, k, ncx, ncy)] -
//           4.*tt_in[get_index_3d(i-1, j, k, ncx, ncy)] +
//              tt_in[get_index_3d(i-2, j, k, ncx, ncy)])/(2.*dx);
//
//        t = tt_in[get_index_3d(i  , j, k, ncx, ncy)] - dx*am;
//
//        a3 = a3<t ? a3 : t;
    }
    
    // =========================================================================
    // EIKONAL SOLVER (same as basic version)
    // =========================================================================

    sort3(&a1, &a2, &a3);
    
    const real_t fh = slowness[idx] * dx;
    
    t = a1 + fh;
    
    if (t > a2) {
        t = 0.5 * (a1 + a2 + sqrt(2.0 * fh * fh - (a1 - a2) * (a1 - a2)));

        if (t > a3) {
            const real_t sum_a = a1 + a2 + a3;
            const real_t term = -2.0 * a1 * a1 + 2.0 * a1 * a2 - 2.0 * a2 * a2 +
                               2.0 * a1 * a3 + 2.0 * a2 * a3 - 2.0 * a3 * a3 +
                               3.0 * fh * fh;

            t = (sum_a + sqrt(term)) / 3.0;
        }
    }
    
    // Update only if smaller
    tt_out[idx] = fmin2(t, tt_in[idx]);
}
