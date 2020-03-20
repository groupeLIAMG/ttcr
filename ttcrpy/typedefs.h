

#ifndef ttcrtypedefs_h
#define ttcrtypedefs_h

#include "ttcr_t.h"
#include "Cell.h"

namespace ttcr {

typedef Node2Dcsp<double, uint32_t> node2d;
typedef Cell<double, node2d, sxz<double>> cell2d;
typedef CellElliptical<double,node2d,sxz<double>> cell2d_e;
typedef CellTiltedElliptical<double,node2d,sxz<double>> cell2d_te;
typedef CellVTI_PSV<double,node2d,sxz<double>> cell2d_p;
typedef CellVTI_SH<double,node2d,sxz<double>> cell2d_h;
}

#endif
