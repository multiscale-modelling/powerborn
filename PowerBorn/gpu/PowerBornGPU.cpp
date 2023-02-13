/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#include <PowerBorn/gpu/PowerBornGPU.h>
#include <opencl_critical.h>
#include <PowerBorn/gpu/ForcefieldParams.h>
#include <PowerBorn/BaseCoordArray.h>
#include <cl_code.h>

#include <sstream>
#include <cstring>
#include <iostream>
#include <string>

#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/transform_width.hpp>
#include <boost/archive/iterators/insert_linebreaks.hpp>
#include <boost/archive/iterators/remove_whitespace.hpp>

namespace powerbornGpu
{
// fit_params taken from fits as described in PowerBorn paper: 0.97734722  0.06491036
PowerBornGpu::PowerBornGpu(unsigned int device_type,bool init_sources):
        ocl(opencl_setup::get_opencl_access(device_type)),
        fit_param_a(0.97734722), fit_param_b(0.06491036),
        do_dump(1), debug_mode(0), do_profile(0),
        max_kernel_time(0.0)
{
    SourcesMap sources_map;
    if(init_sources)
    {
        sources_map.insert(std::make_pair(std::string("gpu_nv_amd"),b64decode(SIMONA::kernel_PowerBorn_gpu_nv_amd)));
        sources_map.insert(std::make_pair(std::string("gpu_intel"),b64decode(SIMONA::kernel_PowerBorn_gpu_intel)));
        sources_map.insert(std::make_pair(std::string("cpu"),b64decode(SIMONA::kernel_PowerBorn_cpu)));
    }

    this->initKernel(sources_map);
    cl_mem_flags f = CL_MEM_READ_WRITE;
    waters.setFlags(f);
    watermask.setFlags(f);
    atomIdBuffer.setFlags(f);
    stack.setFlags(f);
    result.setFlags(f);
    area.setFlags(f);
    energy.setFlags(f);
    buffersizes.setFlags(f);
}

PowerBornGpu::PowerBornGpu(opencl_setup* os, SourcesMap const& kernel_srcs):
        ocl(os), fit_param_a(0.97734722), fit_param_b(0.06491036),
        do_dump(1), debug_mode(0), do_profile(0),
        max_kernel_time(0.0)
{
    this->initKernel(kernel_srcs);
    cl_mem_flags f = CL_MEM_READ_WRITE;
    waters.setFlags(f);
    watermask.setFlags(f);
    atomIdBuffer.setFlags(f);
    stack.setFlags(f);
    result.setFlags(f);
    area.setFlags(f);
    energy.setFlags(f);
    buffersizes.setFlags(f);
}

PowerBornGpu::~PowerBornGpu()
{
    // wait for queue to finish to prevent writing to already deallocated memory when reading from device
    OpenClCriticalLock lock;
    ocl->getqueue().finish();
    lock.release_check();
}

void PowerBornGpu::setFitParams(float a, float b)
{
    fit_param_a = a;
    fit_param_b = b;
    this->initKernel(SourcesMap());
}

void PowerBornGpu::setForcefield(boost::shared_ptr<ForcefieldParams>& ff)
{
    ffparams = ff;
    ffparams->copyToDevice(*ocl);
}

void PowerBornGpu::initKernel(SourcesMap const & ext_kernel_srcs)
{
    std::string kernel_src;
    std::string kernel_src_mod;
    if(ocl->get_device_type() == CL_DEVICE_TYPE_CPU)
    {
#ifdef WITH_INTEL_GPU_CPU_EXPERIMENTAL
        std::cout << "PowerBornGPU using CPU" << std::endl;
        kernel_src += "#define GPU_BLOCK_SIZE 256u\n";
        kernel_src_mod = "cpu";
#else
        std::cout << "Intel GPU and CPU support is experimental. Try it by enabling it in the SConstruct and setting WITH_INTEL_GPU_CPU_EXPERIMENTAL=True. Exiting." << std::endl;
        exit(5);
#endif
    }
    else if(ocl->get_device_vendor().compare(0,6,"NVIDIA",0,6) == 0)
    {
        std::cout << "PowerBornGPU using NVIDIA Platform" << std::endl;
        kernel_src += "#define GPU_BLOCK_SIZE 256u\n";
        kernel_src_mod = "gpu_nv_amd";
    }
    else if( (ocl->get_device_vendor().compare(0,25,"Advanced Micro Devices, Inc.",0,25) == 0) ||
        (ocl->get_device_vendor().compare(0,3,"AMD",0,3) == 0) )
    {
        //We found some implementations, which put a \n at the end of the AMD, Inc. Therefore we only compare the first 25 chars.
        std::cout << "PowerBornGPU using AMD Platform" << std::endl;
        kernel_src += "#define GPU_BLOCK_SIZE 256u\n";
        kernel_src_mod = "gpu_nv_amd";
    }
    else if(ocl->get_device_vendor().compare(0,5,"Intel",0,5) == 0 )
    {
#ifdef WITH_INTEL_GPU_CPU_EXPERIMENTAL
        std::cout << "PowerBornGPU using Intel Platform" << std::endl;
        kernel_src += "#define GPU_BLOCK_SIZE 256u\n";
        kernel_src_mod = "gpu_intel";
#else
        std::cout << "Intel GPU and CPU support is experimental. Try it by enabling it in the SConstruct and setting WITH_INTEL_GPU_CPU_EXPERIMENTAL=True. Exiting." << std::endl;
        exit(5);
#endif
    }
    else
    {
        std::cerr << "Unknown or Unsupported OpenCL Device" << std::endl;
        throw PowerBornGpuException();
    }
    if(ocl->get_max_workgroup_size() < 256)
    {
        std::cerr << "Maximum work group size of the chosen device is too small!" << std::endl;
        throw PowerBornGpuException();
    }

    // fit parameters
    std::stringstream s;
    s << "#define FIT_PARAM_A " << fit_param_a << std::endl;
    s << "#define FIT_PARAM_B " << fit_param_b << std::endl;
    kernel_src += s.str();

    // kernel src, in case of an empty map:
    if(ext_kernel_srcs.empty())
    {
        kernel_src += ocl->program_source_to_string(   (std::string("../gpu/PowerBorn") + kernel_src_mod + std::string(".cl")).c_str()  );
    }
    else
    {
        kernel_src += ext_kernel_srcs.find(kernel_src_mod)->second;
    }
    if(ocl->get_platform_vendor() == "NVIDIA Corporation")
    {
        kernelmap = ocl->build_kernel_map(kernel_src, "-cl-fast-relaxed-math -cl-single-precision-constant -cl-nv-maxrregcount=48");
    }
    else
    {
        kernelmap = ocl->build_kernel_map(kernel_src, "-cl-fast-relaxed-math -cl-single-precision-constant");
    }
    std::cout << "Open CL built " << kernelmap->size() << " kernels from PowerBorn_" << kernel_src_mod << ".cl:" << std::endl;
    for(opencl_setup::KernelMap::iterator it= kernelmap->begin(); it!=kernelmap->end(); ++it)
    {
        std::cout << it->first << std::endl;
    }
    std::cout << "Finished generating Kernels" << std::endl;
}

void PowerBornGpu::initDefaultForceField(unsigned int s, float* charges)
{
    OpenClCriticalLock lock;
    if(!ffparams)
    {
        ffparams.reset(new ForcefieldParams());
    }
    if(charges)
    {
        ffparams->setCharges(charges, s);
    }else
    {
        ffparams->defaultParams(s);
    }
    ffparams->copyToDevice(*ocl);
    lock.release_check();
}

unsigned int PowerBornGpu::getWaterSizeY(unsigned int atomCnt)
{
    return (atomCnt + (this->getLocalSize()>>5))/(this->getLocalSize()>>5);
}

unsigned int PowerBornGpu::getEnergySizeY(unsigned int atomCnt)
{
    return (1 + (atomCnt + this->getLocalSize() - 1) / this->getLocalSize());
}

void PowerBornGpu::resizeDeviceBuffers()
{
    unsigned int s = task_count * constants[0].atomPitch;
    pGPUPoint4d p;
    p.x = p.y = p.z = p.w = 0.0f;
    watermask.resize(*ocl, s * (512 / 32));

    atomIdBuffer.resize(*ocl, task_count * constants[0].atomPitch * this->getIntegrationBlocks() * this->getIntegrationBlockScale());
    buffersizes.resize(task_count * this->getWaterSizeY(constants[0].atomCnt));
    for(unsigned int i=0; i<buffersizes.getHost().size(); ++i)
    {
        buffersizes[i].s[0] = -1;
        buffersizes[i].s[1] = -1;
    }
    result.resize(s*this->getIntegrationBlocks(), 2.0f);
    area.resize(s, 1.0f);
    energy.resize(task_count * 4 * (1 + (constants[0].atomCnt + this->getLocalSize() - 1) / this->getLocalSize()), 0.0f);
}

void PowerBornGpu::updateDeviceBuffers()
{
    this->resizeDeviceBuffers();
    this->prepareConstants();
    constants.writeToDevice(*ocl);
    atoms.writeToDevice(*ocl);
    buffersizes.writeToDevice(*ocl);
    // debug: copy ff params again to make sure they have not been modified by last run or whatever
    ffparams->copyToDevice(*ocl);

    result.reallocDevice(*ocl);
    area.reallocDevice(*ocl);
    energy.reallocDevice(*ocl);
    ocl->getqueue().finish();
}

void PowerBornGpu::prepareConstants()
{
    for(unsigned int j=0; j<task_count; ++j)
    {
        unsigned int offset = constants[j].atomPitch*j;
        unsigned int atomCnt = constants[j].atomCnt;
        float xmax, ymax, zmax, rmax, xmin, ymin, zmin;
        xmax = ymax = zmax = rmax = -1e10f;
        xmin = ymin = zmin = 1e10f;
        for(unsigned int i=0; i<atomCnt; ++i)
        {
            float x1 = atoms[offset + i].x;
            float y1 = atoms[offset + i].y;
            float z1 = atoms[offset + i].z;
            float r1 = atoms[offset + i].w;
            xmax = x1 > xmax ? x1 : xmax;
            ymax = y1 > ymax ? y1 : ymax;
            zmax = z1 > zmax ? z1 : zmax;
            rmax = r1 > rmax ? r1 : rmax;
            xmin = x1 < xmin ? x1 : xmin;
            ymin = y1 < ymin ? y1 : ymin;
            zmin = z1 < zmin ? z1 : zmin;
        }
        constants[j].root.x = 0.5f*(xmax + xmin);
        constants[j].root.y = 0.5f*(ymax + ymin);
        constants[j].root.z = 0.5f*(zmax + zmin);
        constants[j].root.w = 1.0f + 2.0f * rmax + std::max(xmax - xmin, std::max(ymax - ymin, zmax - zmin));
        float box_size = resolution * 128.0f;
        int levels = 7;
        while (box_size < constants[j].root.w && levels <= 10)
        {
            box_size *= 2.0f;
            levels++;
        }
        if (levels <= 10)
        {
            constants[j].root.w = box_size;
        }
        else
        {
            levels = 10;
        }
        constants[j].maxLevel = levels;
        constants[j].maincelllevel = 3 - (levels & 1);
        constants[j].maincellcount = 1 << (3 * (3 - (levels & 1)));
    }
}

void PowerBornGpu::prepareComputation()
{
    OpenClCriticalLock lock;

    //std::cout << "preparing computaiton" << std::endl;
    if(debug_mode & 1)
    {
        std::cout << "PowerBornGpu dump step: " << do_dump << std::endl;
    }
    if((do_dump == 0 && (debug_mode & 1)) || (debug_mode & 2))
    {
        this->dumpStructures("powerborn_gpu_step.dump");
    }
    do_dump++;

    this->updateDeviceBuffers(); // allocate device and/or host memory
    lock.release_check();
}

void PowerBornGpu::doComputation()
{
    OpenClCriticalLock lock;

    //std::cout << "running computaiton for " << task_count << " tasks" << std::endl;
    // now run the kernels
    //
    event_runtime_sum = 0.0;
    this->runWaterKernel();
    this->runBornRadiiKernel();
    this->runEnergyKernel();
    lock.release_check();
}

void PowerBornGpu::finishComputation()
{
    this->finishComputation(final_result);
}

void PowerBornGpu::finishComputation(std::vector<float>& result, bool is_recompute)
{
    OpenClCriticalLock lock;

    // debug readback of energy, make sure kernel finished
    unsigned int eBad = 1;
    while (eBad)
    {
        ocl->getqueue().finish();
        buffersizes.readFromDevice(*ocl);
        ocl->getqueue().finish();
        eBad = 0;
        for(unsigned int j=0; j<task_count * getEnergySizeY(constants[0].atomCnt); ++j)
        {
            eBad += (buffersizes[j].s[0] == -1) ? 1 : 0;
        }
        if (eBad)
        {
            std::cerr << "eBad " << eBad << std::endl;
            throw eBadException();
        }
    }

    // enque data read backs
    try
    {
        //constants.readFromDevice(*ocl);
        //result.readFromDevice(*ocl);
        //area.readFromDevice(*ocl);
        energy.readFromDevice(*ocl);
        ocl->getqueue().finish();
    } catch(cl::Error& err)
    {
        std::cerr << "PowerBornGPU OpenCL error while reading energy from device!" << std::endl
                << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
        this->dumpStructures();
        throw err;
    }

    try
    {
        this->computeFinalResult(result);
        if((debug_mode & 4) && !is_recompute)
        {
            this->recomputeEnergies();
        }
    } catch(cl::Error& err)
    {
        std::cerr << "PowerBornGPU OpenCL error while finishing computation!" << std::endl
                << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
        this->dumpStructures();
        throw err;
    }
    //this->debug();
    //std::cout << "finished computaiton" << std::endl;
    lock.release_check();
}

const float* PowerBornGpu::getEnergiesPtr(unsigned int task)
{
    if(task >= task_count)
    {
        std::cerr << "PowerBornGpu  task energy too large " << task << " " <<  task_count << std::endl;
        throw PowerBornGpuException();
    }
    return &final_result[8 * task];
}

void PowerBornGpu::runWaterKernel()
{
    /* __kernel void
WaterKernel( __global float4* atoms,
           __global pGPUConstants* constants,
           __global const float4* waterSeeds,
           __global int* waterMask)
     */
    cl::Kernel& k = (*kernelmap)["WaterKernel"];
    k.setArg(0, atoms.getBuffer());
    k.setArg(1, constants.getBuffer());
    k.setArg(2, ffparams->watergrid.getBuffer());
    k.setArg(3, watermask.getBuffer());
    k.setArg(4, buffersizes.getBuffer());

    cl::NDRange offset(0, 0);
    cl::NDRange global(this->getGlobalSize(), this->getWaterSizeY(constants[0].atomCnt));
    cl::NDRange local(this->getLocalSize(), 1);
    //std::cout << "running water kernel" << std::endl;
    try
    {
        cl::Event event;
        ocl->getqueue().enqueueNDRangeKernel(k, offset, global, local, 0, &event);
        ocl->getqueue().finish();
        if(do_profile)
        {
            double kernel_time = this->getEventRuntime(event);
            event_runtime_sum += kernel_time;
            max_kernel_time = std::max(max_kernel_time, kernel_time);
        }
    } catch(cl::Error& err)
    {
        std::cerr << "Error while running water kernel!" << std::endl
                << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
        throw err;
    }
}

void PowerBornGpu::runBornRadiiKernel()
{
/*
 * PowerBornKernel( __global float4* atoms,
           __global float4* water,
           __global pGPUConstants* constants,
           __global float* result,
           __global float* resultArea,
           __global const float4* waterSeeds,
           __global const int* waterMask,
           __global int* stack,
           __global int* atomIdBuffer - second only
        )
 */
    int integrationBlocks = this->getIntegrationBlocks();
    int integrationBlocksScale = this->getIntegrationBlockScale();

    // debug readback
    unsigned int cBad = 1;
    while (cBad)
    {
        ocl->getqueue().finish();
        buffersizes.readFromDevice(*ocl);
        ocl->getqueue().finish();
        cBad = 0;
        for(unsigned int j=0; j<task_count * getWaterSizeY(constants[0].atomCnt); ++j)
        {
            cBad += (buffersizes[j].s[0] == -1) ? 1 : 0;
            cBad += (buffersizes[j].s[1] == -1) ? 1 : 0;
        }
        if (cBad)
        {
            std::cerr << "cBad " << cBad << std::endl;
            throw cBadException();
        }
    }

    // stack and water reallocation
    unsigned int waterTotal = 0;
    unsigned int stackTotal = 0;
    for(unsigned int j=0; j<task_count; ++j)
    {
        constants[j].octTreeMaxSize = 0; // reset for debug in energy kernel
        constants[j].waterOffset = waterTotal;
        unsigned int waterCnt = 0;
        unsigned int stackCnt = 0;
        unsigned int start = j * this->getWaterSizeY(constants[0].atomCnt);
        unsigned int stop = (j + 1) * this->getWaterSizeY(constants[0].atomCnt);
        for(unsigned int i=start; i<stop; i++)
        {
            waterCnt += buffersizes[i].s[0];
            stackCnt += buffersizes[i].s[1];
        }
        constants[j].waterCnt = waterCnt;
        unsigned int waterPitch = (waterCnt + 31) & 0xffffffe0;
        waterTotal += waterPitch;

        constants[j].stackOffset = stackTotal + 512*4;
        constants[j].stackSize = 2 * stackCnt + 512 + 64 * constants[j].atomCnt;
        stackTotal += constants[j].stackSize + 512*4;
    }

    constants.writeToDevice(*ocl);
    ocl->getqueue().finish();

    if (waterTotal > waters.size())
    {
        waterTotal = ((unsigned int)(waterTotal*1.2f) + 16383) & (~16383); // add 20% for security, round up to 16384
        std::cout << "Reallocating water buffer from " << waters.size() << " to " << waterTotal << std::endl;
        waters.resize(*ocl, (unsigned int) waterTotal);
    }
    if (stackTotal > stack.size())
    {
        stackTotal = ((unsigned int)(stackTotal*1.2f) + 16383) & (~16383); // add 20% for security, round up to 16384
        std::cout << "Reallocating stack buffer from " << stack.size() << " to " << stackTotal << std::endl;
        stack.resize(*ocl, stackTotal); // resize host buffer
    }

    cl::Kernel& k = (*kernelmap)["PowerBornKernel1"];
    k.setArg(0, atoms.getBuffer());
    k.setArg(1, waters.getBuffer());
    k.setArg(2, constants.getBuffer());
    k.setArg(3, result.getBuffer());
    k.setArg(4, area.getBuffer());
    k.setArg(5, ffparams->watergrid.getBuffer());
    k.setArg(6, watermask.getBuffer());
    k.setArg(7, stack.getBuffer());
    k.setArg(8, integrationBlocks);

    cl::NDRange offset(0, 0);
    cl::NDRange global(this->getGlobalSize(), 1);
    cl::NDRange local(this->getLocalSize(), 1);

    try
    {
        cl::Event event;
        ocl->getqueue().enqueueNDRangeKernel(k, offset, global, local, 0, &event);
        ocl->getqueue().finish();
        if(do_profile)
        {
            double kernel_time = this->getEventRuntime(event);
            event_runtime_sum += kernel_time;
            max_kernel_time = std::max(max_kernel_time, kernel_time);
        }
    } catch(cl::Error& err)
    {
        std::cerr << "Error while running born radii 1 kernel!" << std::endl
                << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
        throw err;
    }

    cl::Kernel& kwi = (*kernelmap)["PowerBornTopLevelWaterIntegral"];
    kwi.setArg(0, atoms.getBuffer());
    kwi.setArg(1, constants.getBuffer());
    kwi.setArg(2, stack.getBuffer());
    kwi.setArg(3, result.getBuffer());
    kwi.setArg(4, atomIdBuffer.getBuffer());
    kwi.setArg(5, integrationBlocksScale);


    cl::NDRange offsetwi(0, 0);
    cl::NDRange globalwi(this->getGlobalSize(), integrationBlocks);
    cl::NDRange localwi(this->getLocalSize(), 1);
    try
    {
        cl::Event event;
        ocl->getqueue().enqueueNDRangeKernel(kwi, offsetwi, globalwi, localwi, 0, &event);
        ocl->getqueue().finish();
        if(do_profile)
        {
            double kernel_time = this->getEventRuntime(event);
            event_runtime_sum += kernel_time;
            max_kernel_time = std::max(max_kernel_time, kernel_time);
        }
    } catch(cl::Error& err)
    {
        std::cerr << "Error while running born radii top level water integral kernel!" << std::endl
                << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
        throw err;
    }

    cl::Kernel& k2 = (*kernelmap)["PowerBornKernel2"];
    k2.setArg(0, atoms.getBuffer());
    k2.setArg(1, waters.getBuffer());
    k2.setArg(2, constants.getBuffer());
    k2.setArg(3, result.getBuffer());
    k2.setArg(4, stack.getBuffer());
    k2.setArg(5, atomIdBuffer.getBuffer());
    k2.setArg(6, integrationBlocksScale);

    cl::NDRange offset2(0,0);
    cl::NDRange global2(this->getGlobalSize() / 2, integrationBlocks);
    cl::NDRange local2(this->getLocalSize() / 2, 1);

    try
    {
        cl::Event event;
        ocl->getqueue().enqueueNDRangeKernel(k2, offset2, global2, local2, 0, &event);
        ocl->getqueue().finish();
        if(do_profile)
        {
            double kernel_time = this->getEventRuntime(event);
            event_runtime_sum += kernel_time;
            max_kernel_time = std::max(max_kernel_time, kernel_time);
        }
    } catch(cl::Error& err)
    {
        std::cerr << "Error while running born radii 2 kernel!" << std::endl
                << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
        throw err;
    }

    cl::Kernel& k3 = (*kernelmap)["PowerBornKernelAggregate"];
    k3.setArg(0, atoms.getBuffer());
    k3.setArg(1, constants.getBuffer());
    k3.setArg(2, result.getBuffer());
    k3.setArg(3, integrationBlocks);

    cl::NDRange offset3(0,0);
    cl::NDRange global3(this->getGlobalSize(),(constants[0].atomCnt+this->getLocalSize()-1)/this->getLocalSize());
    cl::NDRange local3(this->getLocalSize(),1);
    //std::cout << "running born radii kernel" << std::endl;
    try
    {
        cl::Event event;
        ocl->getqueue().enqueueNDRangeKernel(k3, offset3, global3, local3, 0, &event);
        ocl->getqueue().finish();
        if(do_profile)
        {
            double kernel_time = this->getEventRuntime(event);
            event_runtime_sum += kernel_time;
            max_kernel_time = std::max(max_kernel_time, kernel_time);
        }
    } catch(cl::Error& err)
    {
        std::cerr << "Error while running born radii aggregate kernel!" << std::endl
                << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
        throw err;
    }
}

float * PowerBornGpu::get_br()
{
    unsigned int const &  offset = constants[0].atomPitch;
    unsigned int const & atomCnt = constants[0].atomCnt;
    result.readFromDevice(*ocl);
    ocl->getqueue().finish();
    Bradii.resize(constants[0].atomCnt*task_count);
    for(unsigned int j=0; j<task_count; ++j)
    {
        for(unsigned int i=0; i<atomCnt; ++i)
        {
            Bradii[j*atomCnt + i] = result[getIntegrationBlocks()*offset*j + i];
        }
    }
    return &(Bradii[0]);
}

float * PowerBornGpu::get_area()
{
    area.readFromDevice(*ocl);
    ocl->getqueue().finish();
    unsigned int offset = constants[0].atomPitch;
    unsigned int atomCnt = constants[0].atomCnt;
    Surfaces.resize(constants[0].atomCnt*task_count);

    for(unsigned int j=0; j<task_count; ++j)
    {
        for(unsigned int i=0; i<atomCnt; ++i)
        {
            Surfaces[j*atomCnt + i] = area[offset * j + i];
        }
    }
    return &(Surfaces[0]);
}

void PowerBornGpu::runEnergyKernel()
{
    /*EnergyKernel(__global const float4* atoms,
           __global const pGPUEnergyConstants* energyconstants,
           __global const pGPUConstants* constants,
           __global const float* result,
           __global const float* resultArea,
           __global float* resultEnergy,
           __global const float* eps,
           __global const float* sig,
           __global const float* charge,
           __global const int* excludedmask,
           __global const int4* dihdefs,
           __global const float4* dihparams_v,
           __global const float4* dihparams_t,
           __global const float4* dihparams_m
        )
     */

    unsigned int finished_size = task_count * this->getEnergySizeY(constants[0].atomCnt);
    if(buffersizes.getHost().size() < finished_size)
    {
        buffersizes.resize(finished_size);
    }
    for(unsigned int i=0; i<buffersizes.getHost().size(); ++i)
    {
        buffersizes[i].s[0] = -1;
        buffersizes[i].s[1] = -1;
    }
    buffersizes.writeToDevice(*ocl);
    ocl->getqueue().finish();

    cl::Kernel& k = (*kernelmap)["EnergyKernel"];
    k.setArg(0, atoms.getBuffer());
    k.setArg(1, ffparams->energyconstants.getBuffer());
    k.setArg(2, constants.getBuffer());
    k.setArg(3, this->getIntegrationBlocks());
    k.setArg(4, result.getBuffer());
    k.setArg(5, area.getBuffer());
    k.setArg(6, energy.getBuffer());
    k.setArg(7, ffparams->eps.getBuffer());
    k.setArg(8, ffparams->sig.getBuffer());
    k.setArg(9, ffparams->q.getBuffer());
    k.setArg(10, ffparams->excludedmask.getBuffer());
    k.setArg(11, ffparams->dihdefs.getBuffer());
    k.setArg(12, ffparams->dihparam_v.getBuffer());
    k.setArg(13, ffparams->dihparam_t.getBuffer());
    k.setArg(14, ffparams->dihparam_m.getBuffer());
    k.setArg(15, buffersizes.getBuffer());

    cl::NDRange offset(0,0);
    cl::NDRange global(this->getGlobalSize(), this->getEnergySizeY(constants[0].atomCnt));
    cl::NDRange local(this->getLocalSize(),1);
    //std::cout << "running energy kernel" << std::endl;
    try
    {
        cl::Event event;
        ocl->getqueue().enqueueNDRangeKernel(k, offset, global, local, 0, &event);
        ocl->getqueue().finish();
        if(do_profile)
        {
            double kernel_time = this->getEventRuntime(event);
            event_runtime_sum += kernel_time;
            max_kernel_time = std::max(max_kernel_time, kernel_time);
        }
    } catch(cl::Error& err)
    {
        std::cerr << "Error while running energy kernel!" << std::endl
                << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
        throw err;
    }
}

void PowerBornGpu::computeFinalResult(std::vector<float>& final_result)
{
    final_result.resize(8 * task_count);
    unsigned int blocks_per_structure = (constants[0].atomCnt + this->getLocalSize() - 1) / this->getLocalSize();
    unsigned int blocksize = 4 + blocks_per_structure * 4;
    for(unsigned int task=0; task<task_count; ++task)
    {
        unsigned int offset = blocksize * task;
        final_result[task*8 + 0] = 0.0f;
        final_result[task*8 + 1] = 0.0f;
        final_result[task*8 + 2] = 0.0f;
        final_result[task*8 + 3] = 0.0f;
        final_result[task*8 + 4] = energy[offset];
        final_result[task*8 + 5] = 0.0f;
        final_result[task*8 + 6] = 0.0f;
        final_result[task*8 + 7] = 0.0f;
        offset += 4;
        for(unsigned int i=0; i<blocks_per_structure; ++i)
        {
            final_result[task*8 + 0] += energy[offset + i * 4 + 0];
            final_result[task*8 + 1] += energy[offset + i * 4 + 1];
            final_result[task*8 + 2] += energy[offset + i * 4 + 2];
            final_result[task*8 + 3] += energy[offset + i * 4 + 3];
        }
    }
}

double PowerBornGpu::getEventRuntime(cl::Event& event)
{
    cl_ulong ev_start_time = event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
    cl_ulong ev_end_time   = event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
    return 1e-9 * (double)(ev_end_time - ev_start_time);
}

void PowerBornGpu::debug()
{
    constants.readFromDevice(*ocl);
    result.readFromDevice(*ocl);
    area.readFromDevice(*ocl);
    energy.readFromDevice(*ocl);
    ocl->getqueue().finish();
    unsigned int offset = constants[0].atomPitch;
    unsigned int atomCnt = constants[0].atomCnt;
    std::cout << "born radii" << std::endl;
    for(unsigned int j=0; j<task_count; ++j)
    {
        for(unsigned int i=0; i<atomCnt; ++i)
        {
            std::cout << result[getIntegrationBlocks()*offset*j + i] << " "; 
        }
        std::cout << std::endl;
    }
    std::cout << "area" << std::endl;
    for(unsigned int i=0; i<atomCnt; ++i)
    {
        for(unsigned int j=0; j<task_count; ++j)
        {
            std::cout << area[offset * j + i] << " ";
        }
        std::cout << std::endl;
    }
}

void PowerBornGpu::dumpStructures(const std::string& fname)
{
    powerborn::BaseCoordArray coords;
    for(unsigned int i=0; i<atoms.getHost().size(); ++i)
    {
        coords.insert(atoms[i].x, atoms[i].y, atoms[i].z, atoms[i].w);
    }
    coords.dump(fname);
}

void PowerBornGpu::recomputeEnergies()
{
    this->resetConstants();
    this->prepareComputation();
    this->doComputation();
    std::vector<float> comparison;
    this->finishComputation(comparison, true);
    bool all_ok = true;
    for(unsigned int i=0; i<task_count; i++)
    {
        bool structure_ok =true;
        for(unsigned int j=i*8; j<i*8+8; ++j)
        {
            structure_ok &= comparison[j] == final_result[j];
        }
        if(!structure_ok)
        {
            std::cout << "energy mismatch for structure " << i << std::endl;
            for(unsigned int j=i*8; j<i*8+8; ++j)
            {
                std::cout <<  comparison[j] << " " << final_result[j] << std::endl;
            }
            all_ok = false;
            std::cout << "----------------------------------------------" << std::endl;
        }
    }
    if(!all_ok)
    {
        static unsigned int dump_id = 0;
        std::stringstream s;
        s << "powerborngpu_energy_mismatch_" << dump_id << ".dump";
        this->dumpStructures(s.str());
        dump_id++;
        std::cout << "saved dump of structures to " << s.str() << std::endl;
        throw EnergyMismatchException();
    }
}

void PowerBornGpu::setDump(int d)
{
    // count up to 0 and then dump
    do_dump = -d;
}

void PowerBornGpu::setDebugMode(int d)
{
    debug_mode = d;
    std::cout << "PowerBornGpu debug mode set: " << d << std::endl;
}

void PowerBornGpu::setProfiling(int s)
{
    do_profile = s;
    if(do_profile)
    {
        ocl->enable_profiling();
    }
    else
    {
        ocl->disable_profiling();
    }
}

double PowerBornGpu::getMaxKernelTime()
{
    return max_kernel_time;
}


std::string PowerBornGpu::b64decode(const std::string& s)
{
    std::stringstream os;
    namespace bai = boost::archive::iterators;
    typedef bai::transform_width< bai::binary_from_base64<bai::remove_whitespace<const char*> >, 8, 6 > base64_dec;
    unsigned int size = s.size();
    if (size == 0)
        return std::string();
 
    std::copy(base64_dec(s.data()), base64_dec(s.data() + size),
              std::ostream_iterator<char>(os));
    
    return os.str();
}

} /* namespace powerborn */
