/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef POWERBORNGPU_H_
#define POWERBORNGPU_H_

#include <openclsetup.h>
#include <PowerBorn/gpu/Task.h>

#include <exception>
#include <map>

namespace powerbornGpu
{

class PowerBornGpuException: public std::exception
{
};

class EnergyMismatchException: public std::exception
{
};

class cBadException: public std::exception
{
};

class eBadException: public std::exception
{
};

class ForcefieldParams;

class PowerBornGpu: public Task
{
public:
    typedef std::map<std::string const,std::string const> SourcesMap;
protected:
    // OpenCL stuff
    opencl_setup* ocl;
    boost::shared_ptr<opencl_setup::KernelMap> kernelmap;
    boost::shared_ptr<ForcefieldParams> ffparams;
    //std::vector<cl::Event> events;

    double event_runtime_sum;
    
    DeviceBuffer<pGPUPoint4d> waters; // position and radii
    DeviceBuffer<cl_uint> stack, watermask, atomIdBuffer;
    VectorBuffer<float> result, area, energy;
    VectorBuffer<cl_int2> buffersizes;
    std::vector<float> final_result;
    float fit_param_a, fit_param_b;
    int do_dump, debug_mode, do_profile;
    double max_kernel_time;

    void debug();
    void dumpStructures(const std::string& fname="powerborngpu.dump");

    void initKernel(SourcesMap const& ext_kernel_srcs);
    void resizeDeviceBuffers();
    void updateDeviceBuffers();
    void runWaterKernel();
    void runBornRadiiKernel();
    void runEnergyKernel();
    void computeFinalResult(std::vector<float>&);
    void recomputeEnergies();
    void finishComputation(std::vector<float>&, bool is_recompute=false);
    double getEventRuntime(cl::Event&);
    unsigned int getWaterSizeY(unsigned int);
    unsigned int getEnergySizeY(unsigned int);
    std::vector<float> Bradii;
    std::vector<float> Surfaces;

public:
    float * get_br();
    float * get_area();
    double get_event_runtime_sum() { return event_runtime_sum; }
    /* This map typedef should contain the sources with a understood string 
     * Something like:
     * maps["intel-cpu"] = intel cpu code
     * It is only used for initKernel
     */
    
    PowerBornGpu(unsigned int device_type = CL_DEVICE_TYPE_GPU,bool init_sources = false);
    PowerBornGpu(opencl_setup*, SourcesMap const& ext_kernel_srcs);
    virtual ~PowerBornGpu();

    void initDefaultForceField(unsigned int s, float* charges = NULL);
    void setForcefield(boost::shared_ptr<ForcefieldParams>&);
    void prepareConstants();
    void prepareComputation();
    void doComputation();
    void finishComputation();
    void setFitParams(float, float);
    void setDump(int);
    void setDebugMode(int);
    void setProfiling(int);

    const float* getEnergiesPtr(unsigned int task);
    double getMaxKernelTime();

    static std::string b64decode(std::string const & sources);
};

} /* namespace powerbornGpu */

#endif /* POWERBORNGPU_H_ */
