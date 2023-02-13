/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 *    
 *    As this is reference code for your own program,
 *    the code contained in this file test_program.cpp
 *    and only in this file can be copied and modified freely
 *    and can be considered under public domain.
 */

#include <openclsetup.h>

#include <PowerBorn/BaseCoordArray.h>
#include <PowerBorn/gpu/ForcefieldParams.h>
#include <PowerBorn/gpu/Task.h>
#include <PowerBorn/gpu/PowerBornGPU.h>

#ifdef WITH_BENCHMARKS
#include <boost/timer/timer.hpp>
#endif

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cout << "PowerBorn GPU testprogram. Usage: ./test_program infile.pqr. Example files are provided in the PowerBorn distribution inside the testfiles directory." << std::endl;
        exit(1);
    }

    #ifdef WITH_BENCHMARKS
    boost::timer::cpu_timer cputimer;
    #endif

    std::string infile = argv[1];
    powerborn::BaseCoordArray atoms;
    powerborn::Array charges;
    //Initialize the Object looking for a GPU device. Use the OpenCL sources compiled into the static library (true).
    powerbornGpu::PowerBornGpu pbrg(CL_DEVICE_TYPE_GPU,true);

    //A standard PQR parser:
    atoms.parsePqr(infile, &charges);

    if(atoms.size() <= 0)
    {
        std::cout << "Could not parse any atoms from the PQR file. Please compare your PQR with the test PQRs provided in the directory testfiles of the PowerBorn distribution." << std::endl;
        exit(1);
    }
    
    //Simple standard forcefield parameters are set (charges)
    pbrg.initDefaultForceField(atoms.size(), (float*)&charges[0]);

    // 64 calculations are prepared. The more calculations are prepared before running the better.
    // Not all calculations are run in parallel, but the scheduler works better if the queue is big.
    int task_count = 64;
    double cpu_runtime_sum = 0.0;
    double cpu_runtime_sumsq = 0.0;
    for (unsigned int i = 0; i <= 10; i += 1)
    {
        double sum = 0.0;
        double sumsq = 0.0;
        double cpu_sum = 0.0;
        double cpu_sumsq = 0.0;
        double cpu_runtime = 0.0;
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
        #ifdef WITH_BENCHMARKS
        cputimer.start();
        #endif
        //Copies everything to gpu
        pbrg.prepareComputation();
        //starts the queue, but does not guarantee that it ended.
        pbrg.doComputation();
        //guarantees the completion (blocking call)
        pbrg.finishComputation();
        #ifdef WITH_BENCHMARKS
        cputimer.stop();
        cpu_runtime = 1.e-9 * (double)cputimer.elapsed().wall;
        cpu_runtime_sum += cpu_runtime;
        cpu_runtime_sumsq += cpu_runtime * cpu_runtime;
        #endif
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
                  << avg << " , " << stddev << " , " << stddev / avg << ").";
        #ifdef WITH_BENCHMARKS
        std::cout << " Runtime: " << cpu_runtime << " seconds.";
        #endif
        std::cout << std::endl;
    }
    #ifdef WITH_BENCHMARKS
    double avg_runtime = cpu_runtime_sum / 11.0;
    double avg_runtime_stddev = sqrt( ( cpu_runtime_sumsq - cpu_runtime_sum * cpu_runtime_sum / 11.0 ) / 10.0 );
    std::cout << "Average runtime: " << avg_runtime << " +- " << avg_runtime_stddev << " seconds." << std::endl;
    #endif
    return 0;
}
