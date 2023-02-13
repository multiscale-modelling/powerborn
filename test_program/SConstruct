# -*- mode: python -*-
# vim: syntax=python
from __future__ import print_function
import os
import sys


outputdir = "lib"
bindir = "bin"
libname = "PowerBorn_gpu"

WITH_BENCHMARKS=False
CROSS_COMPILE=False


#If you want to cross compile, please setup all paths here:
if CROSS_COMPILE:
    compiler_prefix = 'i686-w64-mingw32'
    #env.Append(CCFLAGS="-m64")
    #env.Append(LINKFLAGS="-m64")
    env.Append(EXE_NAME="_32bit.exe")

    env['CC'] = '%s-gcc' %(compiler_prefix)
    env['CXX'] = '%s-g++' %(compiler_prefix)
    env['RANLIB'] = '%s-ranlib' %(compiler_prefix)
    env['AR'] = '%s-ar' %(compiler_prefix)



##############################
####### OpenCL options #######
##############################
USE_CPU_FOR_OPENCL = False
WITH_OPENCL_DEPRECATED = False
USE_GLOBAL_OPENCL_LIBS = False
SMOG_NO_BONDED = True
WITH_INTEL_GPU_CPU_EXPERIMENTAL = False
##############################


##############################
######### Environment ########
##############################
env = None

if CROSS_COMPILE :
    env = Environment(ENV= os.environ, tools = ["mingw"],toolpath=".")
else:
    env = Environment(ENV= os.environ, tools = ["default"],toolpath=".")

##############################
### Absolute Include paths ###
##############################
env.Append(CPPPATH= ['#/'])
env.Append(CPPPATH= ['#/opencl'])
env.Append(CPPPATH= ['#/PowerBorn/external'])
env.Append(CPPPATH= ['#/external_includes'])
env.Append(LIBPATH= ['#/'+outputdir])
##############################

env.Append(LIBS=['OpenCL'])
env.Append(LIBS=[libname])
##############################

##############################
## Preparing OpenCL kernels ##
##############################
sys.path.append(os.getcwd() + os.sep + "opencl")

from opencl_scons import convert_kernels
convert_kernels(os.getcwd() + os.sep + "opencl")
##############################

##############################
## Add. Libs for Benchmarks ##
##############################

if WITH_BENCHMARKS:
    env.Append(LIBS=['boost_timer'])
    env.Append(LIBS=['boost_chrono'])
    env.Append(LIBS=['boost_system'])
    env.Append(LIBS=['rt'])
    env.Append(CCFLAGS=['-DWITH_BENCHMARKS'])

if WITH_OPENCL_DEPRECATED == True:
    env.Append(CCFLAGS=['-DCL_USE_DEPRECATED_OPENCL_1_1_APIS'])

if WITH_INTEL_GPU_CPU_EXPERIMENTAL == True:
    env.Append(CCFLAGS=['-DWITH_INTEL_GPU_CPU_EXPERIMENTAL'])

libfiles=[ "#/PowerBorn/gpu/PowerBornGPU.cpp",
           "#/PowerBorn/gpu/Task.cpp",
           "#/PowerBorn/gpu/ForcefieldParams.cpp",
           "#/opencl/cl_code.cpp",
           "#/external_sources/openclsetup.cpp" ,
           "#/external_sources/opencl_critical.cpp"]

binfiles = [ "#/test_program/test_program.cpp" ]
try:
    os.makedirs(outputdir)
except:
    pass

try:
    os.makedirs(bindir)
except:
    pass

env.Library(outputdir + "/" + libname,libfiles)
env.Program(bindir + "/test_program",binfiles)
