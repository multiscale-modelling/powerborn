/*
 *    This file is part of the PowerBorn Library and SIMONA
 *    Please refer to the file LICENSE in either the PowerBorn package (LGPL)
 *    for SIMONA (proprietary) for License and Copyright information.
 */

#ifndef OPENCLSETUP_H_
#define OPENCLSETUP_H_

#define __CL_ENABLE_EXCEPTIONS
//#define __NO_STD_VECTOR
//#define __NO_STD_STRING

#include <CL/opencl.hpp>

#include <string>
#include <map>
#include <vector>

#include "boost/shared_ptr.hpp"
#include "boost/noncopyable.hpp"

/*
 * Singleton Pattern for opencl access through method get_opencl_access
 */
class opencl_setup : public boost::noncopyable {
private:
  cl_int err;
  cl::Platform chosen_platform;
  cl::Device chosen_device;
  std::vector<cl::Platform> platforms;
  cl::Context context;
  std::vector<cl::Device> chosen_devices, available_devices;
  std::vector<cl::Program> programs;
  std::vector<cl::CommandQueue> queues;
  unsigned int workgroup_size, max_queue_count, active_queue_id;

  void output_platform_info(cl::Platform &platform);
  void output_device_info(cl::Device &device);
  void init_standalone(unsigned int device_type);
  opencl_setup(unsigned int device_type = CL_DEVICE_TYPE_CPU);

public:
  typedef std::map<std::string, cl::Kernel> KernelMap;
  typedef boost::shared_ptr<KernelMap> KernelMapPtr;

  std::string program_source_to_string(const char *filename);
  void build_program(const char *filename, std::vector<cl::Kernel> &outputvector);
  boost::shared_ptr<KernelMap> build_kernel_map(std::string const &localsource, const std::string options = "");
  void build_program_string(std::string const &localsource, std::vector<cl::Kernel> &outputvector,
                            const std::string options = "-cl-single-precision-constant");
  unsigned int preferred_workgroup_size(cl::Kernel &kernel);
  std::string get_platform_vendor();
  std::string get_device_vendor();
  unsigned int get_device_type();
  unsigned int get_compute_units();
  unsigned int get_max_workgroup_size();
  void enable_profiling();
  void disable_profiling();
  cl::CommandQueue &getqueue();
  cl::Context &getcontext();
  cl::Buffer *new_buffer_object(cl_mem_flags flags, size_t size, void *host_ptr = NULL, cl_int *err = NULL);
#ifdef OPENCL_CPU
  static opencl_setup *get_opencl_access(unsigned int device_type = CL_DEVICE_TYPE_CPU);
#else
  static opencl_setup *get_opencl_access(unsigned int device_type = CL_DEVICE_TYPE_GPU);
#endif
  static const char *getErrorString(cl_int err);
  void test_device(cl::Device &device);
};

#endif /* OPENCLSETUP_H_ */
