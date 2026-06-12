/*
 * OpenCL Kernels for Grid2Drn Fast Sweeping Methods
 *
 * Four kernels covering the two sweep variants (square / non-square cells)
 * and their WENO3 counterparts:
 *   sweep_update_basic_2d    — 1st-order, dx == dz
 *   sweep_update_xz_2d       — 1st-order, dx != dz
 *   sweep_update_weno3_2d    — WENO3,     dx == dz
 *   sweep_update_weno3_xz_2d — WENO3,     dx != dz
 *
 * Node index convention: idx = i*(ncz+1)+j  (x outer, z inner).
 * Four Gauss-Seidel sweep directions: (±i, ±j).
 * Diagonal level: ip+jp ∈ [0, ncx+ncz].
 *
 * Precision is controlled by the build flag -DUSE_DOUBLE=1 (same as the 3D
 * kernels).  Floating-point literals must NOT carry an 'f' suffix so that the
 * double build keeps them double.
 */

// =============================================================================
// TYPE DEFINITIONS
// =============================================================================

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
    typedef float real_t;
    #define REAL_MAX FLT_MAX
    #define WENO_EPS 1.2e-7
#endif

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

inline size_t idx2d(size_t i, size_t j, size_t ncz1) {
    return i * ncz1 + j;
}

inline real_t fmin2(real_t a, real_t b) {
    return (a < b) ? a : b;
}

/*
 * WENO3 upwind approximation of the traveltime at the upwind side of node v2.
 * v0..v4 are five consecutive node values with spacing h.
 * forward=1: returns v2 + h*ap  (uses the right-biased stencil; needs v3,v4)
 * forward=0: returns v2 - h*am  (uses the left-biased stencil;  needs v0,v1)
 * Pass 0.0 for the unused endpoint (v0 when forward=1 at a near-boundary,
 * v4 when forward=0); those slots are never read in the respective branch.
 */
inline real_t weno3_upwind(real_t v0, real_t v1, real_t v2,
                            real_t v3, real_t v4, real_t h, int forward)
{
    const real_t eps = WENO_EPS;
    if (forward) {
        const real_t num = v4 - 2.0 * v3 + v2;
        const real_t den = v3 - 2.0 * v2 + v1;
        const real_t r   = (eps + num * num) / (eps + den * den);
        const real_t w   = 1.0 / (1.0 + 2.0 * r * r);
        const real_t ap  = (1.0 - w) * (v3 - v1) / (2.0 * h)
                         + w * (-v4 + 4.0 * v3 - 3.0 * v2) / (2.0 * h);
        return v2 + h * ap;
    } else {
        const real_t num = v2 - 2.0 * v1 + v0;
        const real_t den = v3 - 2.0 * v2 + v1;
        const real_t r   = (eps + num * num) / (eps + den * den);
        const real_t w   = 1.0 / (1.0 + 2.0 * r * r);
        const real_t am  = (1.0 - w) * (v3 - v1) / (2.0 * h)
                         + w * (3.0 * v2 - 4.0 * v1 + v0) / (2.0 * h);
        return v2 - h * am;
    }
}

/*
 * 2D eikonal solver for square cells (dx == dz).
 * Given upwind x-neighbor a, upwind z-neighbor b, and fh = slowness * dx,
 * returns the updated traveltime.
 */
inline real_t eikonal2d_sq(real_t a, real_t b, real_t fh)
{
    real_t diff = a - b;
    if (diff < 0.0) diff = -diff;
    if (diff >= fh)
        return fmin2(a, b) + fh;
    return 0.5 * (a + b + sqrt(2.0 * fh * fh - (a - b) * (a - b)));
}

/*
 * 2D eikonal solver for non-square cells.
 * dx, dz may differ; uses the full quadratic discriminant.
 */
inline real_t eikonal2d_xz(real_t a, real_t b, real_t s,
                             real_t dx, real_t dz)
{
    if (a < b && (b - a) / dx > s)
        return a + s * dx;
    if (a > b && (a - b) / dz > s)
        return b + s * dz;

    const real_t dx2   = dx * dx;
    const real_t dz2   = dz * dz;
    const real_t s2    = s * s;
    const real_t denom = dx2 + dz2;
    return (b * dx2 + a * dz2) / denom
         + sqrt((2.0 * a * b * dx2 * dz2
                 - a * a * dx2 * dz2
                 - b * b * dx2 * dz2
                 + dx2 * dx2 * dz2 * s2
                 + dx2 * dz2 * dz2 * s2) / (denom * denom));
}

// =============================================================================
// KERNEL 1: BASIC sweep, square cells (dx == dz)
// Corresponds to Grid2Drn::sweep() / update_node()
// =============================================================================

__kernel void sweep_update_basic_2d(
    __global const real_t *tt_in,
    __global       real_t *tt_out,
    __global const real_t *slowness,
    __global const uchar  *frozen,
    const real_t dx,
    const real_t dz,
    const uint   ncx,
    const uint   ncz,
    __global const uint *plane_nodes,
    const uint plane_offset,
    const uint plane_count)
{
    const size_t gid = get_global_id(0);
    if (gid >= plane_count) return;

    const size_t flat = plane_nodes[plane_offset + gid];
    const size_t ncz1 = ncz + 1;
    const size_t i    = flat / ncz1;
    const size_t j    = flat % ncz1;

    if (frozen[flat]) return;

    real_t a, b, t;

    if      (i == 0)   a = tt_in[idx2d(i+1, j, ncz1)];
    else if (i == ncx) a = tt_in[idx2d(i-1, j, ncz1)];
    else {
        a = tt_in[idx2d(i-1, j, ncz1)];
        t = tt_in[idx2d(i+1, j, ncz1)];
        a = fmin2(a, t);
    }

    if      (j == 0)   b = tt_in[idx2d(i, j+1, ncz1)];
    else if (j == ncz) b = tt_in[idx2d(i, j-1, ncz1)];
    else {
        b = tt_in[idx2d(i, j-1, ncz1)];
        t = tt_in[idx2d(i, j+1, ncz1)];
        b = fmin2(b, t);
    }

    t = eikonal2d_sq(a, b, slowness[flat] * dx);
    tt_out[flat] = fmin2(t, tt_in[flat]);
}

// =============================================================================
// KERNEL 2: BASIC sweep, non-square cells (dx != dz)
// Corresponds to Grid2Drn::sweep_xz() / update_node_xz()
// =============================================================================

__kernel void sweep_update_xz_2d(
    __global const real_t *tt_in,
    __global       real_t *tt_out,
    __global const real_t *slowness,
    __global const uchar  *frozen,
    const real_t dx,
    const real_t dz,
    const uint   ncx,
    const uint   ncz,
    __global const uint *plane_nodes,
    const uint plane_offset,
    const uint plane_count)
{
    const size_t gid = get_global_id(0);
    if (gid >= plane_count) return;

    const size_t flat = plane_nodes[plane_offset + gid];
    const size_t ncz1 = ncz + 1;
    const size_t i    = flat / ncz1;
    const size_t j    = flat % ncz1;

    if (frozen[flat]) return;

    real_t a, b, t;

    if      (i == 0)   a = tt_in[idx2d(i+1, j, ncz1)];
    else if (i == ncx) a = tt_in[idx2d(i-1, j, ncz1)];
    else {
        a = tt_in[idx2d(i-1, j, ncz1)];
        t = tt_in[idx2d(i+1, j, ncz1)];
        a = fmin2(a, t);
    }

    if      (j == 0)   b = tt_in[idx2d(i, j+1, ncz1)];
    else if (j == ncz) b = tt_in[idx2d(i, j-1, ncz1)];
    else {
        b = tt_in[idx2d(i, j-1, ncz1)];
        t = tt_in[idx2d(i, j+1, ncz1)];
        b = fmin2(b, t);
    }

    t = eikonal2d_xz(a, b, slowness[flat], dx, dz);
    tt_out[flat] = fmin2(t, tt_in[flat]);
}

// =============================================================================
// KERNEL 3: WENO3 sweep, square cells (dx == dz)
// Corresponds to Grid2Drn::sweep_weno3() / update_node_weno3()
// Only valid when dx == dz (both directions use the same step size).
// =============================================================================

__kernel void sweep_update_weno3_2d(
    __global const real_t *tt_in,
    __global       real_t *tt_out,
    __global const real_t *slowness,
    __global const uchar  *frozen,
    const real_t dx,
    const real_t dz,
    const uint   ncx,
    const uint   ncz,
    __global const uint *plane_nodes,
    const uint plane_offset,
    const uint plane_count)
{
    const size_t gid = get_global_id(0);
    if (gid >= plane_count) return;

    const size_t flat = plane_nodes[plane_offset + gid];
    const size_t ncz1 = ncz + 1;
    const size_t i    = flat / ncz1;
    const size_t j    = flat % ncz1;

    if (frozen[flat]) return;

    real_t a, b, t;

    // ---- X direction (spacing dx) ----
    if (i == 0) {
        a = tt_in[idx2d(i+1, j, ncz1)];
    } else if (i == ncx) {
        a = tt_in[idx2d(i-1, j, ncz1)];
    } else if (i == 1) {
        // Near-left boundary: forward WENO3 + first-order backward
        const real_t v1 = tt_in[idx2d(i-1, j, ncz1)];
        const real_t v2 = tt_in[idx2d(i,   j, ncz1)];
        const real_t v3 = tt_in[idx2d(i+1, j, ncz1)];
        const real_t v4 = tt_in[idx2d(i+2, j, ncz1)];
        a = weno3_upwind(0.0, v1, v2, v3, v4, dx, 1);
        a = fmin2(a, v1);
    } else if (i == ncx - 1) {
        // Near-right boundary: backward WENO3 + first-order forward
        const real_t v0 = tt_in[idx2d(i-2, j, ncz1)];
        const real_t v1 = tt_in[idx2d(i-1, j, ncz1)];
        const real_t v2 = tt_in[idx2d(i,   j, ncz1)];
        const real_t v3 = tt_in[idx2d(i+1, j, ncz1)];
        a = weno3_upwind(v0, v1, v2, v3, 0.0, dx, 0);
        a = fmin2(a, v3);
    } else {
        // Interior: min of forward and backward WENO3
        const real_t v0 = tt_in[idx2d(i-2, j, ncz1)];
        const real_t v1 = tt_in[idx2d(i-1, j, ncz1)];
        const real_t v2 = tt_in[idx2d(i,   j, ncz1)];
        const real_t v3 = tt_in[idx2d(i+1, j, ncz1)];
        const real_t v4 = tt_in[idx2d(i+2, j, ncz1)];
        a = weno3_upwind(v0, v1, v2, v3, v4, dx, 1);
        t = weno3_upwind(v0, v1, v2, v3, v4, dx, 0);
        a = fmin2(a, t);
    }

    // ---- Z direction (spacing dx, same as x because dx==dz) ----
    if (j == 0) {
        b = tt_in[idx2d(i, j+1, ncz1)];
    } else if (j == ncz) {
        b = tt_in[idx2d(i, j-1, ncz1)];
    } else if (j == 1) {
        const real_t v1 = tt_in[idx2d(i, j-1, ncz1)];
        const real_t v2 = tt_in[idx2d(i, j,   ncz1)];
        const real_t v3 = tt_in[idx2d(i, j+1, ncz1)];
        const real_t v4 = tt_in[idx2d(i, j+2, ncz1)];
        b = weno3_upwind(0.0, v1, v2, v3, v4, dx, 1);
        b = fmin2(b, v1);
    } else if (j == ncz - 1) {
        const real_t v0 = tt_in[idx2d(i, j-2, ncz1)];
        const real_t v1 = tt_in[idx2d(i, j-1, ncz1)];
        const real_t v2 = tt_in[idx2d(i, j,   ncz1)];
        const real_t v3 = tt_in[idx2d(i, j+1, ncz1)];
        b = weno3_upwind(v0, v1, v2, v3, 0.0, dx, 0);
        b = fmin2(b, v3);
    } else {
        const real_t v0 = tt_in[idx2d(i, j-2, ncz1)];
        const real_t v1 = tt_in[idx2d(i, j-1, ncz1)];
        const real_t v2 = tt_in[idx2d(i, j,   ncz1)];
        const real_t v3 = tt_in[idx2d(i, j+1, ncz1)];
        const real_t v4 = tt_in[idx2d(i, j+2, ncz1)];
        b = weno3_upwind(v0, v1, v2, v3, v4, dx, 1);
        t = weno3_upwind(v0, v1, v2, v3, v4, dx, 0);
        b = fmin2(b, t);
    }

    t = eikonal2d_sq(a, b, slowness[flat] * dx);
    tt_out[flat] = fmin2(t, tt_in[flat]);
}

// =============================================================================
// KERNEL 4: WENO3 sweep, non-square cells (dx != dz)
// Corresponds to Grid2Drn::sweep_weno3_xz() / update_node_weno3_xz()
// Uses dx for the x WENO stencil and dz for the z WENO stencil.
// =============================================================================

__kernel void sweep_update_weno3_xz_2d(
    __global const real_t *tt_in,
    __global       real_t *tt_out,
    __global const real_t *slowness,
    __global const uchar  *frozen,
    const real_t dx,
    const real_t dz,
    const uint   ncx,
    const uint   ncz,
    __global const uint *plane_nodes,
    const uint plane_offset,
    const uint plane_count)
{
    const size_t gid = get_global_id(0);
    if (gid >= plane_count) return;

    const size_t flat = plane_nodes[plane_offset + gid];
    const size_t ncz1 = ncz + 1;
    const size_t i    = flat / ncz1;
    const size_t j    = flat % ncz1;

    if (frozen[flat]) return;

    real_t a, b, t;

    // ---- X direction (spacing dx) ----
    if (i == 0) {
        a = tt_in[idx2d(i+1, j, ncz1)];
    } else if (i == ncx) {
        a = tt_in[idx2d(i-1, j, ncz1)];
    } else if (i == 1) {
        const real_t v1 = tt_in[idx2d(i-1, j, ncz1)];
        const real_t v2 = tt_in[idx2d(i,   j, ncz1)];
        const real_t v3 = tt_in[idx2d(i+1, j, ncz1)];
        const real_t v4 = tt_in[idx2d(i+2, j, ncz1)];
        a = weno3_upwind(0.0, v1, v2, v3, v4, dx, 1);
        a = fmin2(a, v1);
    } else if (i == ncx - 1) {
        const real_t v0 = tt_in[idx2d(i-2, j, ncz1)];
        const real_t v1 = tt_in[idx2d(i-1, j, ncz1)];
        const real_t v2 = tt_in[idx2d(i,   j, ncz1)];
        const real_t v3 = tt_in[idx2d(i+1, j, ncz1)];
        a = weno3_upwind(v0, v1, v2, v3, 0.0, dx, 0);
        a = fmin2(a, v3);
    } else {
        const real_t v0 = tt_in[idx2d(i-2, j, ncz1)];
        const real_t v1 = tt_in[idx2d(i-1, j, ncz1)];
        const real_t v2 = tt_in[idx2d(i,   j, ncz1)];
        const real_t v3 = tt_in[idx2d(i+1, j, ncz1)];
        const real_t v4 = tt_in[idx2d(i+2, j, ncz1)];
        a = weno3_upwind(v0, v1, v2, v3, v4, dx, 1);
        t = weno3_upwind(v0, v1, v2, v3, v4, dx, 0);
        a = fmin2(a, t);
    }

    // ---- Z direction (spacing dz — differs from dx) ----
    if (j == 0) {
        b = tt_in[idx2d(i, j+1, ncz1)];
    } else if (j == ncz) {
        b = tt_in[idx2d(i, j-1, ncz1)];
    } else if (j == 1) {
        const real_t v1 = tt_in[idx2d(i, j-1, ncz1)];
        const real_t v2 = tt_in[idx2d(i, j,   ncz1)];
        const real_t v3 = tt_in[idx2d(i, j+1, ncz1)];
        const real_t v4 = tt_in[idx2d(i, j+2, ncz1)];
        b = weno3_upwind(0.0, v1, v2, v3, v4, dz, 1);
        b = fmin2(b, v1);
    } else if (j == ncz - 1) {
        const real_t v0 = tt_in[idx2d(i, j-2, ncz1)];
        const real_t v1 = tt_in[idx2d(i, j-1, ncz1)];
        const real_t v2 = tt_in[idx2d(i, j,   ncz1)];
        const real_t v3 = tt_in[idx2d(i, j+1, ncz1)];
        b = weno3_upwind(v0, v1, v2, v3, 0.0, dz, 0);
        b = fmin2(b, v3);
    } else {
        const real_t v0 = tt_in[idx2d(i, j-2, ncz1)];
        const real_t v1 = tt_in[idx2d(i, j-1, ncz1)];
        const real_t v2 = tt_in[idx2d(i, j,   ncz1)];
        const real_t v3 = tt_in[idx2d(i, j+1, ncz1)];
        const real_t v4 = tt_in[idx2d(i, j+2, ncz1)];
        b = weno3_upwind(v0, v1, v2, v3, v4, dz, 1);
        t = weno3_upwind(v0, v1, v2, v3, v4, dz, 0);
        b = fmin2(b, t);
    }

    t = eikonal2d_xz(a, b, slowness[flat], dx, dz);
    tt_out[flat] = fmin2(t, tt_in[flat]);
}
