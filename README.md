# PowerBorn GPU implementation 
PowerBorn calculates the born radii based on the PowerBorn approximation. This library contains the GPU implementation using OpenCL.
The CPU algorithm is explained in the following paper:
http://pubs.acs.org/doi/abs/10.1021/ct300870s
The GPU algorithm is explained in its own paper, pending publishing.

See LICENSE for for licensing information.

## Requirements
1. OpenCL 1.2
2. Boost Libraries > 1.55
3. Recent G++
4. scons
5. Eigen
6. OpenCL headers including CL/opencl.hpp from here: https://raw.githubusercontent.com/KhronosGroup/OpenCL-CLHPP/main/include/CL/opencl.hpp


## Compiling the Library
Compile the library and a test program using scons:
    scons -j4
    ./bin/test_program ./testfiles/1shg.pqr

The test program prepares 64 parallel computations of the protein defined in 1SHG.pqr and calculates them in parallel. The library is placed in lib/PowerBorn_gpu.a (or lib/PowerBorn_gpu.a.lib on Mingw/Windows).

## Linking against the library

It's best to refer to the compilation of the test_program. Short explanation:
1. The Includepath has to point to the top directory of the PowerBorn Library (i.e., your install looks like this INCLUDEDIR/PowerBorn/gpu/files).
2. Link against PowerBorn_gpu.a

## Example code

The following is an excerpt from the file test_program/test_program.cpp.
    powerborn::BaseCoordArray atoms;
    powerborn::Array charges;
    //Initialize the Object looking for a GPU device. Use the OpenCL sources compiled into the static library (true).
    powerbornGpu::PowerBornGpu pbrg(CL_DEVICE_TYPE_GPU,true);

    //A standard PQR parser:
    atoms.parsePqr(infile, &charges);

    //Simple standard forcefield parameters are set (charges)
    pbrg.initDefaultForceField(atoms.size(), (float*)&charges[0]);

    // 64 calculations are prepared. The more calculations are prepared before running the better.
    // Not all calculations are run in parallel, but the scheduler works better if the queue is big.
    int task_count = 64;
    for (unsigned int i = 0; i <= 10; i += 1)
    {
        double sum = 0.0;
        double sumsq = 0.0;
        double cpu_sum = 0.0;
        double cpu_sumsq = 0.0;
        std::cout << i;
        pbrg.resetTasks();
        for (unsigned int j = 0; j<task_count; ++j)
        {
            //The conformation is added to the task queue:
            pbrg.addToTask(atoms);
            // and rotated so that different tasks are calculated 
            // This is a rigid rotation, the born radii or energies 
            // should not change due to it, but the calculation does due to a different octree.
            atoms.randomRotate();
        }
        //Copies everything to gpu
        pbrg.prepareComputation();
        //starts the queue, but does not guarantee that it ended.
        pbrg.doComputation();
        //guarantees the completion (blocking call)
        pbrg.finishComputation();
        //The following code is self-explanatory and returns the individiual born radii per atom:
        float * br = pbrg.get_br(); //get_area also exists with the same api.
        if (br)
        {
            std::cout << "BR values for individual atoms" << std::endl;
            for(int task_num = 0; task_num < task_count; ++task_num)
            {
                for (int atomid = 0; atomid < atoms.size(); atomid++)
                {
                    std::cout << "BR(Task,Atomid): BR(" << task_num << "," << atomid << ") = " << " " << br[task_num*atoms.size() + atomid] << std::endl;
                }
            }
            std::cout << std::endl;
        }
        // The structure getEnergiesPtr(j) has the following contents:
        // ptr[0] == Coulomb
        // ptr[1] == Lennard Jones
        // ptr[2] == Generalized Born
        // ptr[3] == Np Solvation
        // ptr[4] == Dihedral Energy
        // We only output the GB energy, because this is what the library is mainly for:
        for (unsigned int j = 0; j<task_count; ++j)
        {
            double gb_gpu = pbrg.getEnergiesPtr(j)[2];
            sum += gb_gpu;
            sumsq += gb_gpu * gb_gpu;
        }
        double avg = sum / double(task_count);
        double stddev = 0.0; 
        if(task_count > 1) 
        {
            stddev = sqrt(1.0 / (double(task_count) - 1.0) * (sumsq - sum * sum / double(task_count)));
        }
        std::cout << "Result for " << infile << " for " << task_count << " tasks. Energies: (Average,standard deviation,stddev/avg):  (" 
                  << avg << " , " << stddev << " , " << stddev / avg << ")." << std::endl;
    }
