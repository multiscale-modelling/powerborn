/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */
#ifndef PGPU_DEFS_H_
#define PGPU_DEFS_H_

#define __CL_ENABLE_EXCEPTIONS

#include <vector>
#include <iostream>
#include <stdint.h>
#include <CL/opencl.hpp>

namespace powerbornGpu {

template<class T>
struct T4
{
    T x,y,z,w;
    T& operator[](unsigned int s)
    {
        switch(s){
        case 0: return x;
        case 1: return y;
        case 2: return z;
        case 3: return w;
        default:
        {
            std::cerr << "Acessing an element >= 4 in T4 is forbidden!" << std::endl;
            throw std::exception();
        }
        }
    }
    const T& operator[](unsigned int s) const
    {
        switch(s){
        case 0: return x;
        case 1: return y;
        case 2: return z;
        case 3: return w;
        default:
        {
            std::cerr << "Acessing an element >= 4 in T4 is forbidden!" << std::endl;
            throw std::exception();
        }
        }
    }
};

typedef T4<cl_float> pGPUPoint4d;
typedef T4<cl_uint> pGPU4u;

// size of this structure has to be a multiple of 16
struct pGPUConstants
{
    pGPUPoint4d root;
    cl_uint atomCnt;
    cl_uint atomPitch;
    cl_uint waterCnt;
    cl_uint atomsOffset;
    cl_uint waterOffset;
    cl_uint stackOffset;
    cl_uint stackSize;
    //cl_uint waterSeedsCnt;
    cl_uint octTreeMaxSize;
    cl_float integration_factor;
    cl_uint maxLevel;
    cl_uint maincelllevel;
    cl_uint maincellcount; // count of non-zero main cells. this is the ammount of octtree nodes to start integrating with
};

struct pGPUEnergyConstants
{
    cl_uint dihCnt;
    cl_float coulomb_param;
    cl_float gb_param;
    cl_float npse_param;
    cl_float excluded_param_co;
    cl_float excluded_param_lj;
    cl_float r1, r2;
};

struct pPGUOctNode64
{
    pGPUPoint4d p;
    cl_float r1,r2,r3,r4; // reserved to make the size == 96b
    //        int allocsize;
    //        int allocstart;
    //        int localnodeids;
    //        int nextnodeids;
    //        int localnodecnt;
    //        int treenodeid;

    cl_int celltype[72];
};
} // end of namespace

#endif
