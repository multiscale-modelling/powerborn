/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#include <PowerBorn/gpu/ForcefieldParams.h>
#include <cmath>

#ifndef M_PI
#define M_PI 3.141592653589793
#endif


namespace powerbornGpu
{

void ForcefieldParams::reweightGrid(std::vector<pGPUPoint4d>& grid, float factor)
{
    // now normalize the weights so that the actual sum really is 4 pi on the unit sphere
    float weight_sum = 0.0f;
    for(unsigned int i=0; i<grid.size(); ++i)
    {
        weight_sum += grid[i].w;
    }
    float reweight = 4.0 * M_PI / (weight_sum) * factor;
    std::cout << "reweight grid with " << grid.size() << " points and reweight factor of  " << reweight << std::endl;
    for(unsigned int i=0; i<grid.size(); ++i)
    {
        grid[i].w *= reweight;
    }
}

void ForcefieldParams::createGridSymmetric(std::vector<pGPUPoint4d>& grid, bool large)
{
    /*
     * 3 way symmetric grid, either 32 or 64 points per 8th of a unit sphere
     * so in total either 256 or 512 points
     */
    grid.clear();
    pGPUPoint4d p;
    p.x = 0;
    p.y = 0;
    p.z = 0;
    p.w = 0;
    const float maxrad = 1.0f;
    float cutoff;

    // empiric values to get the right grid sizes
    if(large)
        cutoff = 0.15f;
    else
        cutoff = 0.22f;

    float theta_0 = acos(1.0f - cutoff * cutoff/(2.0f * maxrad * maxrad));
    unsigned int lines = 0.5f * M_PI / theta_0;
    float delta_theta = 0.5f * M_PI / float(lines);
    for(unsigned int a=0; a < lines; ++a)
    {
        float theta = (0.5f + float(a)) * delta_theta;
        float phi_0=acos(1.0f-cutoff*cutoff/(2.0f*sin(theta)*sin(theta)*maxrad*maxrad));
        unsigned int pointnum= 0.5f * M_PI / phi_0;
        pointnum = pointnum < 2 ? 2 : pointnum;
        float delta_phi = 0.5f * M_PI / float(pointnum);
        for(unsigned int b=0; b < pointnum; ++b)
        {
            float phi=delta_phi * (float(b) + 0.5f);
            p.x = sin(theta)*cos(phi);
            p.y = sin(theta)*sin(phi);
            p.z = cos(theta);
            p.w = fabs(sin(theta) * delta_phi * delta_theta);
            grid.push_back(p);
        }
    }
    p.x = 0.0f;
    p.y = 0.0f;
    p.z = 1.0f;
    p.w = 0.0f;
    grid.push_back(p); // symmetric points are the same, but since point has no weight, the surface area is not affect. Born radii are also not affected by duplicate water spheres
    if(large)
        grid.push_back(p); // dummy point with no weight to get grid size 64 in large case
    this->reweightGrid(grid, 0.125f);
}

void ForcefieldParams::initWaterGrid()
{
    std::vector<pGPUPoint4d> grid_part;
    this->createGridSymmetric(grid_part, false);
    watergrid.resize(8*(64 + 32));
    for(int i = 0; i < 32; ++ i)
    {
        pGPUPoint4d p = grid_part[i];
        watergrid[i] = p;
        p.x = -p.x;
        watergrid[i + 32] = p;
        p.x = -p.x;
        p.y = -p.y;
        watergrid[i + 2 * 32] = p;
        p.y = -p.y;
        p.z = -p.z;
        watergrid[i + 3 * 32] = p;
        p.z = -p.z;
        p.x = -p.x;
        p.y = -p.y;
        watergrid[i + 4 * 32] = p;
        p.y = -p.y;
        p.z = -p.z;
        watergrid[i + 5 * 32] = p;
        p.x = -p.x;
        p.y = -p.y;
        watergrid[i + 6 * 32] = p;
        p.x = -p.x;
        watergrid[i + 7 * 32] = p;
    }
    this->createGridSymmetric(grid_part, true);
    for(int i = 0; i < 64; ++ i)
    {
        pGPUPoint4d p = grid_part[i];
        watergrid[i + 8 * 32] = p;
        p.x = -p.x;
        watergrid[i + 64 + 8 * 32] = p;
        p.x = -p.x;
        p.y = -p.y;
        watergrid[i + 2 * 64 + 8 * 32] = p;
        p.y = -p.y;
        p.z = -p.z;
        watergrid[i + 3 * 64 + 8 * 32] = p;
        p.z = -p.z;
        p.x = -p.x;
        p.y = -p.y;
        watergrid[i + 4 * 64 + 8 * 32] = p;
        p.y = -p.y;
        p.z = -p.z;
        watergrid[i + 5 * 64 + 8 * 32] = p;
        p.x = -p.x;
        p.y = -p.y;
        watergrid[i + 6 * 64 + 8 * 32] = p;
        p.x = -p.x;
        watergrid[i + 7 * 64 + 8 * 32] = p;
    }
}

ForcefieldParams::ForcefieldParams()
{
    //set all OpenCL memory flags
    cl_mem_flags f = CL_MEM_READ_ONLY;
    energyconstants.setFlags(f);
    watergrid.setFlags(f);
    dihparam_v.setFlags(f);
    dihparam_t.setFlags(f);
    dihparam_m.setFlags(f);
    dihdefs.setFlags(f);
    eps.setFlags(f);
    sig.setFlags(f);
    q.setFlags(f);
    excludedmask.setFlags(f);

    // create water grids
    this->initWaterGrid();
}

void ForcefieldParams::copyToDevice(opencl_setup& ocl)
{
    energyconstants.writeToDevice(ocl);
    watergrid.writeToDevice(ocl);
    dihparam_v.writeToDevice(ocl);
    dihparam_t.writeToDevice(ocl);
    dihparam_m.writeToDevice(ocl);
    dihdefs.writeToDevice(ocl);
    eps.writeToDevice(ocl);
    sig.writeToDevice(ocl);
    q.writeToDevice(ocl);
    excludedmask.writeToDevice(ocl);
    // explicit synchronisation
    ocl.getqueue().finish();
}

void ForcefieldParams::initConstants()
{
    pGPUEnergyConstants& c = energyconstants.getHost();
    c.coulomb_param = 0.5f * 1389.3545660404325f / 4.184f ;
    c.gb_param = -c.coulomb_param * (1.0f - 1.0f / 80.0f);
    c.npse_param = 0.00542f;
    c.excluded_param_co = 0.8333f;
    c.excluded_param_lj = 0.5f;
    c.r1 = c.r2 = 1.0f;
}

void ForcefieldParams::defaultParams(unsigned int s)
{
    pGPUEnergyConstants& c = energyconstants.getHost();
    c.coulomb_param = 0.5 * 1389.3545660404325f / 4.184f;
    c.gb_param = - (1.0f - 1.0f/80.0f) * c.coulomb_param;
    c.npse_param = 0.00542f;
    c.excluded_param_co = 0.8333f;
    c.excluded_param_lj = 0.5f;
    c.r1 = c.r2 = 1.0f;

    unsigned int ps = (s + 31) & (~31);
    eps.resize(ps, 0.1f);
    sig.resize(ps, 1.0f);
    q.resize(ps, 0.1f);
    for(unsigned int i=s; i<ps; ++i)
    {
        eps[i] = 0.0f;
        sig[i] = 0.0f;
        q[i] = 0.0f;
    }

    // make sure s is rounded up to multiple of 32
    //std::cout << "padded size " << s << " " << ps << std::endl;
    excludedmask.resize(ps * ps / 16, 0); // 2 bits for each atom
    for(unsigned int i=0; i<s; ++i)
    {
        for(unsigned int j=0; j<s; ++j)
        {
            unsigned int id = i * ps + j;
            unsigned int flag = i == j;
            unsigned int shift = (id % 16) * 2;
            excludedmask[id/16] |= (flag << shift);
        }
        for(unsigned int j=s; j<ps; ++j) // exclude padded
        {
            unsigned int id = i * ps + j;
            unsigned int flag = 1;
            unsigned int shift = (id % 16) * 2;
            excludedmask[id/16] |= (flag << shift);
        }
    }
    // test excluded
    /*for(unsigned int i=0; i<s; ++i)
    {
        for(unsigned int j=0; j<s; ++j)
        {
            unsigned int id = i * ps + j;
            unsigned int mask = excludedmask[id/16];
            mask >>= ((id % 16) * 2);
            if(mask & 3)
            {
                std::cout << i << " " << j << std::endl;
            }
        }
    }*/

    // dihedrals
    c.dihCnt = s - 4;
    dihdefs.resize(c.dihCnt * 4);
    dihparam_v.resize(c.dihCnt * 4);
    dihparam_t.resize(c.dihCnt * 4);
    dihparam_m.resize(c.dihCnt * 4);
    for(unsigned int i=0; i<c.dihCnt; i++)
    {
        dihdefs[i].x = i;
        dihdefs[i].y = i+1;
        dihdefs[i].z = i+2;
        dihdefs[i].w = i+3;

        dihparam_v[i].x = 1.0f;
        dihparam_v[i].y = 1.0f;
        dihparam_v[i].z = 1.0f;
        dihparam_v[i].w = 1.0f;

        dihparam_t[i].x = 0.0f;
        dihparam_t[i].y = 0.0f;
        dihparam_t[i].z = 0.0f;
        dihparam_t[i].w = 0.0f;

        dihparam_m[i].x = 1.0f;
        dihparam_m[i].y = 2.0f;
        dihparam_m[i].z = 3.0f;
        dihparam_m[i].w = 4.0f;
    }
}

void ForcefieldParams::setCharges(float *c, unsigned int s)
{
    this->defaultParams(s);
    for(unsigned int i=0; i< s; ++i)
    {
        q[i] = c[i];
    }
}

} /* namespace powerbornGpu */
