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
    float param_a = atof(argv[1]);
    float param_b = atof(argv[2]);
    powerborn::BaseCoordArray atoms;
    powerborn::Array charges;
    powerbornGpu::PowerBornGpu pbrg(CL_DEVICE_TYPE_GPU);
    pbrg.setFitParams(param_a, param_b);
    for(int i=3; i<argc; ++i)
    {
        std::string infile(argv[i]);
        atoms.parsePqr(infile, &charges);
        pbrg.initDefaultForceField(atoms.size(), (float*) &charges[0]);

        pbrg.resetTasks();
        pbrg.addToTask(atoms);

        pbrg.prepareComputation();
        pbrg.doComputation();
        pbrg.finishComputation();

        float gb_gpu = pbrg.getEnergiesPtr(0)[2];
        std::cout << "GB ENERGY " << infile << " " << gb_gpu << std::endl;

    }
    return 0;
}
