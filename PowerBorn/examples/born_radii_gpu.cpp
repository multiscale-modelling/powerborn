/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#include "Atoms.h"
#include "gpu/ForcefieldParams.h"
#include "gpu/Task.h"
#include "gpu/PowerBornGPU.h"

int main(int argc, char* argv[])
{
    if(argc != 2) {
        std::cout << "Usage:\n ./born_radii_gpu <input.pqr>\n\nComputes and prints Born radii and GB energies for ACC and FAST.\n" << std::endl;
        abort();
    }
    std::string infile = argv[1];
    std::cout << "pGPUconstants size: " << sizeof(powerbornGpu::pGPUConstants) << std::endl;
    std::cout << "pGPUEnergyConstants size: " << sizeof(powerbornGpu::pGPUEnergyConstants) << std::endl;
    std::cout << "pPGUOctNode64 size: " << sizeof(powerbornGpu::pPGUOctNode64) << std::endl;

    powerborn::BaseCoordArray atoms;
    powerborn::Array charges;
    atoms.parsePqr(infile, &charges);

    powerbornGpu::PowerBornGpu pbrg(CL_DEVICE_TYPE_GPU);
    pbrg.initDefaultForceField(atoms.size(), (float*) &charges[0]);

    pbrg.resetTasks();
    unsigned int s_count = 8;
    for(unsigned int i=0; i<s_count; ++i)
    {
        pbrg.addToTask(atoms);
    }


    pbrg.prepareComputation();
    pbrg.doComputation();
    pbrg.finishComputation();
    for(unsigned int i=0; i<s_count; ++i)
    {
        const float* energies = pbrg.getEnergiesPtr(i);
        for(unsigned int j=0; j<8; ++j)
        {
            std::cout << energies[j] << " ";
        }
        std::cout << std::endl;
    }
    // copy ff to gpu
    return 0;
}
