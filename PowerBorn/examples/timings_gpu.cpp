/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#include "Atoms.h"
#include "gpu/ForcefieldParams.h"
#include "gpu/Task.h"
#include "gpu/PowerBornGPU.h"
#include "OpenMpUtil.h"

int main(int argc, char* argv[])
{
    std::string infile = argv[1];
    powerborn::BaseCoordArray atoms;
    powerborn::Array charges;
    powerbornGpu::PowerBornGpu pbrg(CL_DEVICE_TYPE_GPU);
    atoms.parsePqr(infile, &charges);
    pbrg.initDefaultForceField(atoms.size(), (float*) &charges[0]);

    for(unsigned int i=8; i<=32; i+=8)
    {
        std::cout << "Tasks: " << i << std::endl;
        for(unsigned int k=0; k<10; ++k)
        {
            pbrg.resetTasks();
            for(unsigned int j=0; j<i; ++j)
            {
                pbrg.addToTask(atoms);
            }
            {
                pbrg.prepareComputation();
                {
                    openmp::ScopedTimer t2("without transfer");
                    pbrg.doComputation();
                    pbrg.finishComputation();
                }
            }
        }
    }
    return 0;
}
