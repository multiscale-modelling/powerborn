/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef TEST_INTEGRATOR_H_
#define TEST_INTEGRATOR_H_

#include <PowerBorn/gpu/Task.h>
#include <PowerBorn/gpu/Integrator.h>

namespace powerbornGpu
{

class TestIntegrator: public Task
{
    opencl_setup* ocl;
    boost::shared_ptr<opencl_setup::KernelMap> kernelmap;
    VectorBuffer<float> result;
    powerborn::Array compare_results;
    std::vector<cl::Event> events;

    void initKernel(const std::string& options)
    {
        ocl = opencl_setup::get_opencl_access(CL_DEVICE_TYPE_GPU);
        std::string kernel_src;
        kernel_src += "#define GPU_BLOCK_SIZE 256\n";
        // kernel src
        kernel_src += ocl->program_source_to_string("../gpu/Powerborn_amd_bb.cl");
        kernelmap = ocl->build_kernel_map(kernel_src, options);
        std::cout << "Open CL built " << kernelmap->size() << " kernels from PowerBorn_amd_bb.cl:" << std::endl;
        for(opencl_setup::KernelMap::iterator it= kernelmap->begin(); it!=kernelmap->end(); ++it)
        {
            std::cout << it->first << std::endl;
        }
        std::cout << "Finished generating Kernels" << std::endl;

        cl_mem_flags f = CL_MEM_READ_WRITE;
        result.setFlags(f);
    }
    void initTest(const powerborn::BaseCoordArray& coord_array)
    {
        events.clear();
        this->resetTasks();
        this->addToTask(coord_array);
        result.resize(coord_array.size(), 0.0f);
        result.reallocDevice(*ocl);

        events.push_back(constants.writeToDevice(*ocl));
        events.push_back(atoms.writeToDevice(*ocl));

    }
    void runBoundingBoxKernel()
    {
        cl::Kernel& k = (*kernelmap)["BoundingBoxKernel"];
        k.setArg(0, atoms.getBuffer());
        k.setArg(1, constants.getBuffer());

        cl::NDRange offset(0);
        cl::NDRange global(256);
        cl::NDRange local(256);
        cl::Event event;

        ocl->getqueue().enqueueNDRangeKernel(k, offset, global, local, &events, &event);
        events.push_back(event);

    }
    void runPowerBornKernel()
    {
        cl::Kernel k = (*kernelmap)["PowerBornBBKernel"];
        k.setArg(0, atoms.getBuffer());
        k.setArg(1, constants.getBuffer());
        k.setArg(2, result.getBuffer());

        cl::NDRange offset(0);
        cl::NDRange global(256);
        cl::NDRange local(256);
        cl::Event event;

        ocl->getqueue().enqueueNDRangeKernel(k, offset, global, local, &events, &event);
        events.push_back(event);
        events.push_back(result.readFromDevice(*ocl, &events));
        events.push_back(constants.readFromDevice(*ocl, &events));
        cl::WaitForEvents(events);
    }
    void runCpuVersion(const powerborn::BaseCoordArray& atoms)
    {
        compare_results.resize(atoms.size());
        powerborn::Integrator integrator(10.0f, 1.0f, 0.0f, 1.0f, 0.0f);
        powerborn::Coord3 center;
        center[0] = constants[0].root.x;
        center[1] = constants[0].root.y;
        center[2] = constants[0].root.z;
        float size = constants[0].root.w;
        for(unsigned int i=0; i<atoms.size(); ++i)
        {
            compare_results[i] = integrator.cubeIntegral(atoms.getPos(i),center, size);
        }
    }
    void compareResults(const powerborn::BaseCoordArray& atoms)
    {
        powerborn::Coord3 center; 
        center[0] = constants[0].root.x;
        center[1] = constants[0].root.y;
        center[2] = constants[0].root.z;

        for(unsigned int i=0; i<compare_results.size(); ++i)
        {
            if(fabs((compare_results[i] - result[i]) / compare_results[i]) > 1e-6)
            {
                std::cout << i << " " << atoms.x()[i] - center[0] << " "
                    << atoms.y()[i] - center[1] << " "
                    << atoms.z()[i] - center[2] << " "
                    << result[i] << " " << compare_results[i] << " "
                    << (compare_results[i] - result[i]) / compare_results[i] << std::endl;
            }
        }
        std::cout << "center x " << constants[0].root.x << " "
                << "center y " << constants[0].root.y << " "
                << "center z " << constants[0].root.z << " "
                << "center size " << constants[0].root.w << std::endl;
    }
public:
    void runTest(powerborn::BaseCoordArray& atoms, const std::string& options)
    {
        this->initKernel(options);
        this->initTest(atoms);
        this->runBoundingBoxKernel();
        this->runPowerBornKernel();
        this->runCpuVersion(atoms);
        this->compareResults(atoms);
    }
    void runTest2(const std::string& options)
    {
        powerborn::BaseCoordArray atoms;
        for(unsigned int i=0; i<100; ++i)
        {
            atoms.insert(0.0f, 0.0f, i * 0.1f, 2.0f);
        }
        for(unsigned int i=0; i<100; ++i)
        {
            atoms.insert(0.0f, i * 0.1f, i * 0.1f, 2.0f);
        }
        for(unsigned int i=0; i<100; ++i)
        {
            atoms.insert(i * 0.1f, i * 0.1f, i * 0.1f, 2.0f);
        }
        this->initKernel(options);
        this->initTest(atoms);
        constants[0].root.x = 0.0f;
        constants[0].root.y = 0.0f;
        constants[0].root.z = 0.0f;
        constants[0].root.w = 10.1f;
        events.push_back(constants.writeToDevice(*ocl));
        this->runPowerBornKernel();
        this->runCpuVersion(atoms);
        this->compareResults(atoms);

    }
};

} // end of namespace





#endif /* TEST_INTEGRATOR_H_ */
