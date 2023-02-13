/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#include "Atoms.h"
#include "gpu/ForcefieldParams.h"
#include "gpu/Task.h"
#include "gpu/PowerBornGPU.h"
#include "PowerBornMS2.h"

int main(int argc, char* argv[])
{
    std::string infile = argv[1];
    powerborn::BaseCoordArray atoms;
    powerborn::Array charges;
    powerborn::PowerBornMS2 pbr = powerborn::PowerBornMS2::getAcc();
    powerbornGpu::PowerBornGpu pbrg(CL_DEVICE_TYPE_GPU);
    for(int i=1; i<argc; ++i)
    {
        std::string infile(argv[i]);
        atoms.parsePqr(infile, &charges);
        pbr.update(atoms);
        pbrg.initDefaultForceField(atoms.size(), (float*) &charges[0]);

        pbrg.resetTasks();
        pbrg.addToTask(atoms);

        pbrg.prepareComputation();
        pbrg.doComputation();
        pbrg.finishComputation();

        float sasa_gpu = pbrg.getEnergiesPtr(0)[3];
        float sasa_cpu = pbr.getSasa().sum();
        std::cout << infile << " " << sasa_gpu << " " << sasa_cpu << std::endl;

    }
    return 0;
}
