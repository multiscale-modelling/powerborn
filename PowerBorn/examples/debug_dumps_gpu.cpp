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
    if(argc != 3) {
        std::cout << "Usage:\n ./debug_dumps_gpu <input.pqr> <input.dump>\n" << std::endl;
        abort();
    }

    // parse input
    const std::string infile1 = argv[1];
    const std::string infile2 = argv[2];
    powerborn::BaseCoordArray atoms, structures;
    powerborn::Array charges;
    atoms.parsePqr(infile1, &charges);
    structures.parseDump(infile2);

    powerbornGpu::PowerBornGpu pbrg1(CL_DEVICE_TYPE_GPU);
    powerbornGpu::PowerBornGpu pbrg2(CL_DEVICE_TYPE_GPU);
    pbrg1.initDefaultForceField(atoms.size(), (float*) &charges[0]);
    pbrg2.initDefaultForceField(atoms.size(), (float*) &charges[0]);

    if(structures.size() % atoms.size() != 0)
    {
        std::cerr << "Atom count in input.dump is not a multiple of atom count in input.pqr!" << std::endl;
        abort();
    }

    // now split structures and add to tasks;
    pbrg1.resetTasks();
    std::vector<powerborn::BaseCoordArray, Eigen::aligned_allocator<powerborn::BaseCoordArray> > split_atoms;
    unsigned int atom_count = 0;
    while(atom_count < structures.size())
    {
        split_atoms.push_back(powerborn::BaseCoordArray());
        powerborn::BaseCoordArray& s = split_atoms.back();
        s.setSize(atoms.size());
        for(unsigned int i=0; i<atoms.size(); ++i)
        {
            s.x()[i] = structures.x()[atom_count + i];
            s.y()[i] = structures.y()[atom_count + i];
            s.z()[i] = structures.z()[atom_count + i];
            s.r()[i] = structures.r()[atom_count + i];
        }
        pbrg1.addToTask(s);
        atom_count += atoms.size();
    }

    // now compute all structures together
    pbrg1.prepareComputation();
    pbrg1.doComputation();
    pbrg1.finishComputation();

    // now compute single structures and compare result
    for(unsigned int i=0; i<split_atoms.size(); ++i)
    {
        pbrg2.resetTasks();
        pbrg2.addToTask(split_atoms[i]);
        pbrg2.prepareComputation();
        pbrg2.doComputation();
        pbrg2.finishComputation();
        const float* res1 = pbrg1.getEnergiesPtr(i);
        const float* res2 = pbrg2.getEnergiesPtr(0);
        bool all_ok = 1;
        
        for(unsigned int j=0; j<8; ++j)
        {
            std::cout << res1[j] - res2[j] << " ";
            all_ok &= (res1[j] == res2[j]);
        }
        std::cout << std::endl;
        if(!all_ok)
        {
            std::cerr << "comparison for structure " << i << " failed!" << std::endl;
        }
    }

    return 0;
}
