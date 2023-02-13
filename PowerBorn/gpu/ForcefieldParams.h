/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef FORCEFIELDPARAMS_H_
#define FORCEFIELDPARAMS_H_

#include <PowerBorn/gpu/pGPU_defs.h>
#include <PowerBorn/gpu/Buffer.h>

namespace powerbornGpu
{

/*
 * Data that never changes during evaluations(HOST and DEVICE)
 */
class ForcefieldParams
{
    void createGridSymmetric(std::vector<pGPUPoint4d>& grid, bool large);
    void reweightGrid(std::vector<pGPUPoint4d>& grid, float factor);
    void initWaterGrid();
public:
    //this stuff is public for easier access. I am too lazy to write accessors for all of these
    Buffer<pGPUEnergyConstants> energyconstants;
    VectorBuffer<pGPUPoint4d> watergrid, dihparam_v, dihparam_t, dihparam_m;
    VectorBuffer<float> eps, sig, q;
    VectorBuffer<uint32_t> excludedmask;
    VectorBuffer<pGPU4u> dihdefs;

    ForcefieldParams();
    void copyToDevice(opencl_setup& ocl);
    void defaultParams(unsigned int s);
    void setCharges(float* , unsigned int s);
    void initConstants();
};

} /* namespace powerbornGpu */

#endif /* FORCEFIELDPARAMS_H_ */
