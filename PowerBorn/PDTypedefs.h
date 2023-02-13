/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef PDTYPEDEFS_H_
#define PDTYPEDEFS_H_

#include "Eigen/StdVector"
#include "Eigen/StdDeque"
#include "Eigen/Core"
#include "Eigen/Geometry"

#ifndef NO_POWERDIAGRAM
#include "power_diagram.h"
#endif

namespace powerborn {


#ifndef NO_POWERDIAGRAM
typedef POWER_DIAGRAM::PowerDiagram<float, Eigen::Vector3f, 3> PowerDiagram;
typedef PowerDiagram::cell PowerDiagramCell;
typedef PowerDiagram::vertex PowerDiagramVertex;
#endif


} // end of namespace

#endif /* PDTYPEDEFS_H_ */
