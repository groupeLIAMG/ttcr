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

/*
 * Godunov upwind eikonal solve at one node, with a spacing per axis.
 *
 * Solves sum_S ((T - a_i)/h_i)^2 = s^2 over the active set S of axes, those
 * whose upwind neighbour the front has already passed.  With w_i = 1/h_i^2 the
 * root is T = (B + sqrt(B^2 - A*C))/A for A = sum w_i, B = sum a_i*w_i,
 * C = sum a_i^2*w_i - s^2, and the active set is found by ordering the axes and
 * widening from one to three while the result exceeds the next value.
 *
 * Two exact-arithmetic rearrangements keep it accurate: the solve is done in
 * the increment u = T - a1, and the discriminant is formed by the Lagrange
 * identity  B^2 - A*C = A*s^2 - sum_{i<j} w_i*w_j*(b_i - b_j)^2  rather than
 * literally, B^2 and A*C being large and nearly equal.  Mirrors
 * Grid3Drn::solve_godunov, whose comment carries the detail.
 *
 * The axes may be given in any order; each spacing must travel with its own
 * value, which is why this sorts internally rather than taking sorted input.
 */
inline real_t solve_godunov(real_t a1, real_t h1, real_t a2, real_t h2,
                            real_t a3, real_t h3, real_t s)
{
    if (a1 > a2) { swap(&a1, &a2); swap(&h1, &h2); }
    if (a1 > a3) { swap(&a1, &a3); swap(&h1, &h3); }
    if (a2 > a3) { swap(&a2, &a3); swap(&h2, &h3); }

    const real_t b2 = a2 - a1;
    const real_t b3 = a3 - a1;

    // one axis active
    const real_t u1 = s * h1;
    if (u1 <= b2) {
        return a1 + u1;
    }

    // two axes active
    const real_t w1 = 1.0 / (h1 * h1);
    const real_t w2 = 1.0 / (h2 * h2);
    real_t A = w1 + w2;
    real_t B = b2 * w2;
real_t d = A * s * s - w1 * w2 * b2 * b2;
if (d < 0.0) d = 0.0;
const real_t u2 = (B + sqrt(d)) / A;
if (u2 <= b3) {
    return a1 + u2;
}

    // three axes active
    const real_t w3 = 1.0 / (h3 * h3);
    const real_t b23 = b2 - b3;
    A += w3;
    B += b3 * w3;
    d = A * s * s - (w1 * w2 * b2 * b2 +
                     w1 * w3 * b3 * b3 +
                     w2 * w3 * b23 * b23);
    if (d < 0.0) {
        return a1 + u2;
    }
    return a1 + (B + sqrt(d)) / A;
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
    
    // Cubic cells take the cheaper closed forms below; anything else goes to
    // the general per-axis solve, which needs the values UNSORTED so each one
    // keeps the spacing of its own axis.  The test is exact rather than
    // toleranced: being wrong towards solve_godunov costs a little arithmetic,
    // being wrong towards the cubic form returns wrong traveltimes.  dx, dy and
    // dz are kernel arguments uniform across the work-group, so this branch
    // does not diverge.  Mirrors Grid3Drn::sweep_auto.
    const real_t s = slowness[idx];
    if (!(dx == dy && dy == dz)) {
        tt_out[idx] = fmin2(solve_godunov(a1, dz, a2, dy, a3, dx, s), tt_in[idx]);
        return;
    }

    // Sort values
    sort3(&a1, &a2, &a3);

    // Compute slowness * grid_spacing
    const real_t fh = s * dx;

    // Solve the eikonal equation in the increment u = T - a1, building each
    // discriminant from increment-sized terms rather than from the a_i
    // themselves.  Both are exact-arithmetic identities, and both matter
    // because the a_i are absolute traveltimes that can dwarf the increment
    // they differ by; written with the a_i directly, as this was, the update
    // loses digits in proportion to the traveltime.  Mirrors
    // Grid3Drn::update_node, whose comment carries the detail.
    const real_t b2 = a2 - a1;
    const real_t b3 = a3 - a1;

    real_t u = fh;

    if (u > b2) {
        // 2D update
        const real_t d2 = 2.0 * fh * fh - b2 * b2;
        if (d2 >= 0.0) {
            u = 0.5 * (b2 + sqrt(d2));

            if (u > b3) {
                // 3D update
                const real_t b23 = b2 - b3;
                const real_t d3 = 3.0 * fh * fh - b2 * b2 - b3 * b3 - b23 * b23;
                if (d3 >= 0.0) {
                    u = (b2 + b3 + sqrt(d3)) / 3.0;
                }
            }
        }
    }
    t = a1 + u;
    
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
        
        a1 = weno3_upwind(v0, v1, v2, v3, v4, dz, true);
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
        
        a1 = weno3_upwind(v0, v1, v2, v3, v4, dz, false);
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
        
                a1 = weno3_upwind(v0, v1, v2, v3, v4, dz, true);
                t = weno3_upwind(v0, v1, v2, v3, v4, dz, false);
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
        
        a2 = weno3_upwind(v0, v1, v2, v3, v4, dy, true);
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
        
        a2 = weno3_upwind(v0, v1, v2, v3, v4, dy, false);
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
        
        a2 = weno3_upwind(v0, v1, v2, v3, v4, dy, true);
        t = weno3_upwind(v0, v1, v2, v3, v4, dy, false);
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

    // Cubic cells take the cheaper closed forms below; anything else goes to
    // the general per-axis solve, which needs the values UNSORTED so each one
    // keeps the spacing of its own axis.  The test is exact rather than
    // toleranced: being wrong towards solve_godunov costs a little arithmetic,
    // being wrong towards the cubic form returns wrong traveltimes.  dx, dy and
    // dz are kernel arguments uniform across the work-group, so this branch
    // does not diverge.  Mirrors Grid3Drn::sweep_auto.
    const real_t s = slowness[idx];
    if (!(dx == dy && dy == dz)) {
        tt_out[idx] = fmin2(solve_godunov(a1, dz, a2, dy, a3, dx, s), tt_in[idx]);
        return;
    }

    // Sort values
    sort3(&a1, &a2, &a3);

    // Compute slowness * grid_spacing
    const real_t fh = s * dx;

    // Solve the eikonal equation in the increment u = T - a1, building each
    // discriminant from increment-sized terms rather than from the a_i
    // themselves.  Both are exact-arithmetic identities, and both matter
    // because the a_i are absolute traveltimes that can dwarf the increment
    // they differ by; written with the a_i directly, as this was, the update
    // loses digits in proportion to the traveltime.  Mirrors
    // Grid3Drn::update_node, whose comment carries the detail.
    const real_t b2 = a2 - a1;
    const real_t b3 = a3 - a1;

    real_t u = fh;

    if (u > b2) {
        // 2D update
        const real_t d2 = 2.0 * fh * fh - b2 * b2;
        if (d2 >= 0.0) {
            u = 0.5 * (b2 + sqrt(d2));

            if (u > b3) {
                // 3D update
                const real_t b23 = b2 - b3;
                const real_t d3 = 3.0 * fh * fh - b2 * b2 - b3 * b3 - b23 * b23;
                if (d3 >= 0.0) {
                    u = (b2 + b3 + sqrt(d3)) / 3.0;
                }
            }
        }
    }
    t = a1 + u;
    
    // Update only if smaller
    tt_out[idx] = fmin2(t, tt_in[idx]);
}
