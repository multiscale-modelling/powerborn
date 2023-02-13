/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */


#include "gpu/Test_Integrator.h"
#include "gpu/PowerBornGPU.h"


int main(int argc, char* argv[])
{
    std::string infile = argv[1];
    std::string options;
    for(int i=2; i<argc; ++i)
    {
        options += std::string(argv[i]);
    }

    powerborn::BaseCoordArray atoms;
    atoms.parsePqr(infile);

    powerbornGpu::TestIntegrator ti;
    ti.runTest(atoms, options);
    std::cout << std::endl;
    ti.runTest2(options);

}
