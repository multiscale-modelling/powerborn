/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef ATOMS_H_
#define ATOMS_H_

#if defined(__SSE__) && !defined(POWERBORN_NO_SSE)
#include "Atoms_sse.h"
#else
#include "Atoms_default.h"
#endif

#endif /* ATOMS_H_ */
