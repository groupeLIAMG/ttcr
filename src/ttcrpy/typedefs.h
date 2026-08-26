#ifndef ttcr_typedefs_h
#define ttcr_typedefs_h

#include "ttcr_t.h"
#include "Node2Dcsp.h"
#include "Node3Dcsp.h"
#include "Cell.h"

namespace ttcr {

typedef Node2Dcsp<double, uint32_t> node2d;
typedef Cell<double, node2d, sxz<double>> cell2d;
typedef CellElliptical<double,node2d,sxz<double>> cell2d_e;
typedef CellTiltedElliptical<double,node2d,sxz<double>> cell2d_te;
typedef CellVTI_PSV<double,node2d,sxz<double>> cell2d_p;
typedef CellVTI_SH<double,node2d,sxz<double>> cell2d_h;
typedef CellTTI_PSV<double,node2d,sxz<double>> cell2d_tp;
typedef CellTTI_SH<double,node2d,sxz<double>> cell2d_th;
typedef CellWeaklyAnelliptical<double,node2d,sxz<double>> cell2d_wa;

// The 3D grids come in both precisions, so each anisotropic cell class has a
// double and a float typedef; the float ones carry an f on the dimension.
typedef Node3Dcsp<double, uint32_t> node3d;
typedef CellElliptical3D<double,node3d,sxyz<double>> cell3d_e;
typedef CellVTI_PSV3D<double,node3d,sxyz<double>> cell3d_p;
typedef CellVTI_SH3D<double,node3d,sxyz<double>> cell3d_h;
typedef CellWeaklyAnelliptical3D<double,node3d,sxyz<double>> cell3d_wa;

typedef Node3Dcsp<float, uint32_t> node3df;
typedef CellElliptical3D<float,node3df,sxyz<float>> cell3df_e;
typedef CellVTI_PSV3D<float,node3df,sxyz<float>> cell3df_p;
typedef CellVTI_SH3D<float,node3df,sxyz<float>> cell3df_h;
typedef CellWeaklyAnelliptical3D<float,node3df,sxyz<float>> cell3df_wa;
}

#endif
