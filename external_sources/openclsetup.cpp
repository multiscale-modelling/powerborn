/*
 *    This file is part of the PowerBorn Library and SIMONA
 *    Please refer to the file LICENSE in either the PowerBorn package (LGPL)
 *    for SIMONA (proprietary) for License and Copyright information.
 */

#include <openclsetup.h>
#include <opencl_critical.h>

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>

void opencl_setup::output_platform_info(cl::Platform &platform) {
  OpenClCriticalLock lock;
  std::cout << "OpenCL Platform Info" << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "CL_PLATFORM_VENDOR            " << platform.getInfo<CL_PLATFORM_VENDOR>() << std::endl;
  std::cout << "CL_PLATFORM_PROFILE           " << platform.getInfo<CL_PLATFORM_PROFILE>() << std::endl;
  std::cout << "CL_PLATFORM_VERSION           " << platform.getInfo<CL_PLATFORM_VERSION>() << std::endl;
  std::cout << "CL_PLATFORM_NAME              " << platform.getInfo<CL_PLATFORM_NAME>() << std::endl;
  std::cout << "CL_PLATFORM_EXTENSIONS        " << platform.getInfo<CL_PLATFORM_EXTENSIONS>() << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
}

void opencl_setup::output_device_info(cl::Device &device) {
  OpenClCriticalLock lock;
  std::cout << "OpenCL Device Info" << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "CL_DEVICE_NAME                " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
  std::cout << "CL_DEVICE_VENDOR              " << device.getInfo<CL_DEVICE_VENDOR>() << std::endl;
  std::cout << "CL_DEVICE_AVAILABLE           " << device.getInfo<CL_DEVICE_AVAILABLE>() << std::endl;
  std::cout << "CL_DEVICE_GLOBAL_MEM_SIZE     " << device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() << std::endl;
  std::cout << "CL_DEVICE_LOCAL_MEM_SIZE      " << device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() << std::endl;

  std::cout << "CL_DEVICE_MAX_CONSTANT_ARGS   " << device.getInfo<CL_DEVICE_MAX_CONSTANT_ARGS>() << std::endl;
  std::cout << "CL_DEVICE_MAX_MEM_ALLOC_SIZE  " << device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>() << std::endl;
  std::cout << "CL_DEVICE_MAX_WORK_GROUP_SIZE " << device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() << std::endl;
  std::cout << "CL_DEVICE_MAX_WORK_ITEM_SIZES " << device.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>()[0] << " "
            << device.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>()[1] << " "
            << device.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>()[2] << std::endl;
  std::cout << "CL_DRIVER_VERSION             " << device.getInfo<CL_DRIVER_VERSION>() << std::endl;
  std::cout << "CL_DEVICE_MAX_COMPUTE_UNITS   " << device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << std::endl;
  std::cout << "CL_DEVICE_MAX_CLOCK_FREQUENCY " << device.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << std::endl;
  std::cout << "CL_DEVICE_EXTENSIONS          " << device.getInfo<CL_DEVICE_EXTENSIONS>() << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
}

void opencl_setup::init_standalone(unsigned int device_type) {
  OpenClCriticalLock lock;
  cl_uint first_usable_platform_index = 0;
  available_devices.clear();
  chosen_devices.clear();
  cl::Platform::get(&platforms);

  std::cout << "Searching devices in " << platforms.size() << " OpenCL platforms" << std::endl;

  for (cl_uint i = 0; i < platforms.size(); ++i) {
    try {
      std::cout << "OpenCL platform " << i << " contains: " << std::endl;
      this->output_platform_info(platforms[i]);
    } catch (cl::Error &err) {
      std::cerr << "Error while getting Platform info for platform " << i << ": " << std::endl
                << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl
                << "Trying to get devices for this platform nevertheless..." << std::endl;
    }
    try {
      available_devices.clear();
      platforms[i].getDevices(device_type, &available_devices);
      for (cl_uint dev_id = 0; dev_id < available_devices.size(); ++dev_id) {
        try {
          this->output_device_info(available_devices[dev_id]);
#ifdef WITH_OPENCL_TEST
          this->test_device(available_devices[dev_id]);
#endif // WITH_OPENCL_TEST
        } catch (cl::Error &err) {
          std::cerr << "Error in getting device info for platform " << i << " device nr. " << dev_id << std::endl
                    << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
        }
      }
    } catch (cl::Error &err) {
      std::cerr << "Error while getting devices for platform " << i << ": " << std::endl
                << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
      // Do not choose this platform
      if (i == first_usable_platform_index) {
        if (++first_usable_platform_index >= platforms.size()) {
          std::cerr << "No usable Device found. Exiting." << std::endl;
          throw err;
        }
      }
    }
  }

  chosen_platform = platforms[first_usable_platform_index];
  std::cout << "Chosen platform index is: " << first_usable_platform_index << std::endl;
  try {
    chosen_platform.getDevices(device_type, &available_devices);
    if (available_devices.size()) {
      // just use first device
      chosen_device = available_devices[0];
      chosen_devices.push_back(chosen_device);
      context = cl::Context(chosen_devices);
      queues.push_back(cl::CommandQueue(context, chosen_device));
      std::cout << "Command Queue constructed" << std::endl;
      this->output_device_info(chosen_device);
    } else {
      throw cl::Error(CL_DEVICE_NOT_FOUND);
    }
  } catch (cl::Error &err) {
    std::cerr << "Error while getting device of chosen platform and initializing context and queue: " << std::endl
              << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
    throw err;
  }
  lock.release_check();
}

opencl_setup::opencl_setup(unsigned int device_type)
    : err(CL_SUCCESS), workgroup_size(64), max_queue_count(1), active_queue_id(0) {
  this->init_standalone(device_type);
}

std::string opencl_setup::program_source_to_string(const char *filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cout << "Program Source " << filename << " not found." << std::endl;
  }
  return std::string(std::istreambuf_iterator<char>(file), (std::istreambuf_iterator<char>()));
}

void opencl_setup::build_program(const char *filename, std::vector<cl::Kernel> &outputvector) {
  std::string localsource = program_source_to_string(filename);
  build_program_string(localsource, outputvector);
}

boost::shared_ptr<opencl_setup::KernelMap> opencl_setup::build_kernel_map(std::string const &localsource,
                                                                          const std::string options) {
  OpenClCriticalLock lock;
  std::cout << "OpenCl build options: " << options << std::endl;
  programs.push_back(cl::Program(context, localsource));
  try {
    programs.back().build(chosen_devices, options.c_str());
  } catch (cl::Error &err) {
    std::cerr << localsource << std::endl;
    std::string buildstring;
    programs.back().getBuildInfo(chosen_device, static_cast<cl_program_build_info>(CL_PROGRAM_BUILD_LOG), &buildstring);
    std::cerr << buildstring.c_str() << std::endl;
    std::cerr << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
    throw err;
  } catch (std::exception &err) {
    std::cout << "unknown exception" << std::endl;
    throw err;
  }
  std::vector<cl::Kernel> Kernels;
  try {
    programs.back().createKernels(&Kernels);
  } catch (cl::Error &err) {
    std::cerr << "Error while OpenCL program was creating Kernels!" << std::endl
              << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
    throw err;
  }

  boost::shared_ptr<KernelMap> kernelmap;
  kernelmap.reset(new KernelMap());
  for (unsigned int i = 0; i < Kernels.size(); ++i) {
    std::string kernel_info_string;
    Kernels[i].getInfo(static_cast<cl_kernel_info>(CL_KERNEL_FUNCTION_NAME), &kernel_info_string);
    kernelmap->insert(std::make_pair(std::string(kernel_info_string.c_str()), Kernels[i]));
  }
  lock.release_check();
  return kernelmap;
}

void opencl_setup::build_program_string(std::string const &localsource, std::vector<cl::Kernel> &outputvector,
                                        const std::string options) {
  OpenClCriticalLock lock;
  programs.push_back(cl::Program(context, localsource));
  try {
    programs.back().build(chosen_devices, options.c_str());
  } catch (cl::Error &err) {
    std::string buildstring;
    programs.back().getBuildInfo(chosen_device, static_cast<cl_program_build_info>(CL_PROGRAM_BUILD_LOG), &buildstring);
    std::cerr << buildstring.c_str() << std::endl;
    std::cerr << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
    exit(err.err());
  } catch (...) {
    std::cout << "unknown exception" << std::endl;
  }
  std::vector<cl::Kernel> Kernels;
  try {
    programs.back().createKernels(&Kernels);
  } catch (cl::Error &err) {
    std::cerr << "Error while OpenCL program was creating Kernels!" << std::endl
              << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
    throw err;
  }
  KernelMap kernelmap;

  // sort lexicographically
  for (unsigned int i = 0; i < Kernels.size(); ++i) {
    std::string kernel_info_string;
    Kernels[i].getInfo(static_cast<cl_kernel_info>(CL_KERNEL_FUNCTION_NAME), &kernel_info_string);
    kernelmap.insert(std::make_pair(std::string(kernel_info_string.c_str()), Kernels[i]));
  }
  for (KernelMap::const_iterator it = kernelmap.begin(); it != kernelmap.end(); ++it) {
    outputvector.push_back(it->second);
  }
  lock.release_check();
}

unsigned int opencl_setup::preferred_workgroup_size(cl::Kernel &kernel) {
  OpenClCriticalLock lock;
  int error = 0;
  size_t workgroup = 0;
  error = kernel.getWorkGroupInfo(chosen_device, CL_KERNEL_WORK_GROUP_SIZE, &workgroup);
  if (workgroup == 0 || error != 0) {
    std::cout << error << " <-ERROR WORKGROUP_SIZE-> " << workgroup << std::endl;
  } else {
    std::cout << "Set Workgroup size to: " << workgroup << std::endl;
  }
  lock.release_check();
  return workgroup;
}

std::string opencl_setup::get_platform_vendor() {
  OpenClCriticalLock lock;
  std::string tmp = chosen_platform.getInfo<CL_PLATFORM_VENDOR>();
  lock.release_check();
  return tmp;
}

std::string opencl_setup::get_device_vendor() {
  OpenClCriticalLock lock;
  std::string tmp = chosen_device.getInfo<CL_DEVICE_VENDOR>();
  lock.release_check();
  return tmp;
}

unsigned int opencl_setup::get_device_type() {
  OpenClCriticalLock lock;
  unsigned int tmp = chosen_device.getInfo<CL_DEVICE_TYPE>();
  lock.release_check();
  return tmp;
}

unsigned int opencl_setup::get_compute_units() {
  OpenClCriticalLock lock;
  unsigned int tmp = chosen_device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
  lock.release_check();
  return tmp;
}

unsigned int opencl_setup::get_max_workgroup_size() {
  OpenClCriticalLock lock;
  unsigned int tmp = chosen_device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
  lock.release_check();
  return tmp;
}

void opencl_setup::enable_profiling() {
  OpenClCriticalLock lock;
  queues[0].finish();
  queues[0] = cl::CommandQueue(context, chosen_device, CL_QUEUE_PROFILING_ENABLE, &err);
  std::cout << "OpenCL CommandQueue profiling enabled" << std::endl;
  lock.release_check();
}

void opencl_setup::disable_profiling() {
  OpenClCriticalLock lock;
  queues[0].finish();
  queues[0] = cl::CommandQueue(context, chosen_device, 0, &err);
  std::cout << "OpenCL CommandQueue profiling disabled" << std::endl;
  lock.release_check();
}

cl::CommandQueue &opencl_setup::getqueue() { return queues[0]; }

cl::Context &opencl_setup::getcontext() { return context; }

cl::Buffer *opencl_setup::new_buffer_object(cl_mem_flags flags, size_t size, void *host_ptr, cl_int *err) {
  OpenClCriticalLock lock;
  cl::Buffer *newbo = new cl::Buffer(this->context, flags, size, host_ptr, err);
  lock.release_check();
  return newbo;
}

opencl_setup *opencl_setup::get_opencl_access(unsigned int device_type) {
  static const unsigned int init_device_type(device_type);
  static opencl_setup ocl(device_type);
  if (device_type != init_device_type) {
    std::cerr << "Opencl_setup was requested with a different device type than it was initialized! Different device "
                 "types in one run are not supported."
              << std::endl;
    throw std::exception(); // for some reason compiler complains about other exception types being not defined!
  }
  return &ocl;
}

const char *opencl_setup::getErrorString(cl_int err) {
  switch (err) {
  case 0:
    return "CL_SUCCESS";
  case -1:
    return "CL_DEVICE_NOT_FOUND";
  case -2:
    return "CL_DEVICE_NOT_AVAILABLE";
  case -3:
    return "CL_COMPILER_NOT_AVAILABLE";
  case -4:
    return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
  case -5:
    return "CL_OUT_OF_RESOURCES";
  case -6:
    return "CL_OUT_OF_HOST_MEMORY";
  case -7:
    return "CL_PROFILING_INFO_NOT_AVAILABLE";
  case -8:
    return "CL_MEM_COPY_OVERLAP";
  case -9:
    return "CL_IMAGE_FORMAT_MISMATCH";
  case -10:
    return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
  case -11:
    return "CL_BUILD_PROGRAM_FAILURE";
  case -12:
    return "CL_MAP_FAILURE";
  case -30:
    return "CL_INVALID_VALUE";
  case -31:
    return "CL_INVALID_DEVICE_TYPE";
  case -32:
    return "CL_INVALID_PLATFORM";
  case -33:
    return "CL_INVALID_DEVICE";
  case -34:
    return "CL_INVALID_CONTEXT";
  case -35:
    return "CL_INVALID_QUEUE_PROPERTIES";
  case -36:
    return "CL_INVALID_COMMAND_QUEUE";
  case -37:
    return "CL_INVALID_HOST_PTR";
  case -38:
    return "CL_INVALID_MEM_OBJECT";
  case -39:
    return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
  case -40:
    return "CL_INVALID_IMAGE_SIZE";
  case -41:
    return "CL_INVALID_SAMPLER";
  case -42:
    return "CL_INVALID_BINARY";
  case -43:
    return "CL_INVALID_BUILD_OPTIONS";
  case -44:
    return "CL_INVALID_PROGRAM";
  case -45:
    return "CL_INVALID_PROGRAM_EXECUTABLE";
  case -46:
    return "CL_INVALID_KERNEL_NAME";
  case -47:
    return "CL_INVALID_KERNEL_DEFINITION";
  case -48:
    return "CL_INVALID_KERNEL";
  case -49:
    return "CL_INVALID_ARG_INDEX";
  case -50:
    return "CL_INVALID_ARG_VALUE";
  case -51:
    return "CL_INVALID_ARG_SIZE";
  case -52:
    return "CL_INVALID_KERNEL_ARGS";
  case -53:
    return "CL_INVALID_WORK_DIMENSION";
  case -54:
    return "CL_INVALID_WORK_GROUP_SIZE";
  case -55:
    return "CL_INVALID_WORK_ITEM_SIZE";
  case -56:
    return "CL_INVALID_GLOBAL_OFFSET";
  case -57:
    return "CL_INVALID_EVENT_WAIT_LIST";
  case -58:
    return "CL_INVALID_EVENT";
  case -59:
    return "CL_INVALID_OPERATION";
  case -60:
    return "CL_INVALID_GL_OBJECT";
  case -61:
    return "CL_INVALID_BUFFER_SIZE";
  case -62:
    return "CL_INVALID_MIP_LEVEL";
  case -63:
    return "CL_INVALID_GLOBAL_WORK_SIZE";
  default:
    return "Unknown OpenCL error";
  }
}

void opencl_setup::test_device(cl::Device &device) {
  OpenClCriticalLock lock;
  std::string kernel_src = "__kernel void \
compute1(__global const float4* a, \
        __global float4* w, \
        __global float* r1, \
        __global float* r2, \
        __global float* r3, \
        const int integrationBlockCount \
)\
{\
    __local int ss[16];\
    __local int vv[16];\
    barrier(CLK_LOCAL_MEM_FENCE| CLK_GLOBAL_MEM_FENCE);\
}\
\
__kernel void \
compute2(__global const float4* a, \
        __global float4* w, \
        __global float* r1, \
        const int integrationBlockCount \
)\
{\
}";
  queues.clear();
  try {
    chosen_device = device;
    chosen_devices.clear();
    chosen_devices.push_back(chosen_device);
    context = cl::Context(chosen_devices);
    queues.push_back(cl::CommandQueue(context, chosen_device));
    std::cout << "TEST Command Queue constructed" << std::endl;
    boost::shared_ptr<KernelMap> km = this->build_kernel_map(kernel_src, "");
    if (km->size() == 2) {
      std::cout << "TEST successfull: " << km->size() << std::endl;
    } else {
      std::cout << "TEST failed: " << km->size() << std::endl;
    }
  } catch (cl::Error &err) {
    std::cerr << "TEST Error:" << std::endl
              << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
  }
  queues.clear();
  lock.release_check();
}
