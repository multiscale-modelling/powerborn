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
    powerbornGpu::PowerBornGpu pbrg(CL_DEVICE_TYPE_GPU);
    powerborn::PowerBornMS2 pbr = powerborn::PowerBornMS2::getAcc();
    atoms.parsePqr(infile, &charges);
    pbrg.initDefaultForceField(atoms.size(), (float*) &charges[0]);
    double sum = 0.0;
    double sumsq = 0.0;
    double cpu_sum = 0.0;
    double cpu_sumsq = 0.0;
    unsigned int count = 0;
    for(unsigned int i=0; i<=1000; i+=32)
    {
        std::cout << i << std::endl;
        pbrg.resetTasks();
        for(unsigned int j=0; j<32; ++j)
        {
            pbrg.addToTask(atoms);
            pbr.update(atoms);
            double cpu_energy = pbr.gbEnergy(atoms, charges);
            cpu_sum += cpu_energy;
            cpu_sumsq += cpu_energy * cpu_energy;
            atoms.randomRotate();
        }
        pbrg.prepareComputation();
        pbrg.doComputation();
        pbrg.finishComputation();
        for(unsigned int j=0; j<32; ++j)
        {
            double gb_gpu = pbrg.getEnergiesPtr(j)[2];
            std::cout << "gb energy " << gb_gpu << std::endl;
            sum += gb_gpu;
            sumsq += gb_gpu * gb_gpu;
        }
        count += 32;
    }
    {
        double avg = sum / double(count);
        double stddev = sqrt(1.0 / (double(count) - 1.0) * (sumsq - sum * sum / double(count)));
        std::cout << infile << " " << count << " " << avg << " " << stddev << " " << stddev /avg << std::endl;
    }
    {
        double avg = cpu_sum / double(count);
        double stddev = sqrt(1.0 / (double(count) - 1.0) * (cpu_sumsq - cpu_sum * cpu_sum / double(count)));
        std::cout << infile << " " << count << " " << avg << " " << stddev << " " << stddev /avg << std::endl;
    }
    return 0;
}
