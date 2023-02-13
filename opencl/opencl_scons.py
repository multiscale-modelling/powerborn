import os
import glob
import base64
import shutil

def convert_kernels(dir):

    # the current powerborn opencl kernel is located in a different dir, so first copy it here
    pbrg_kernel = os.path.join(dir, "../external/PowerBorn/gpu/PowerBorn_amd_mod.cl")
    if(os.path.isfile(pbrg_kernel)):
        shutil.copy(pbrg_kernel, dir)

    # now convert all kernels
    openclfiles = glob.glob(dir + os.sep + "*.cl")
    filedata = {}
   
    for clfile in openclfiles:
        kernel = os.path.basename(clfile)[:-3]
        f = open(clfile, 'rb')
        filedata[kernel] = base64.b64encode(f.read()).decode("utf-8")
        f.close()

    cpp = open(dir + os.sep + "cl_code.cpp", "wt")
    hpp = open(dir + os.sep + "cl_code.h", "wt")

    cpp.write("""
//This is a auto-generated file; changes will be lost !!

#include <cl_code.h>

namespace SIMONA
{
""")

    hpp.write("""
//This is a auto-generated file; changes will be lost !!
#ifndef __CLCODEINC__
#define __CLCODEINC__

#include <string>

namespace SIMONA
{
""")
    
    for filename, base64data in filedata.items():
        cpp.write('extern const std::string kernel_%s = "%s";' % (filename, base64data))
        hpp.write('extern const std::string kernel_%s;\n' % (filename))

    cpp.write("}\n")
    hpp.write("}\n")
    hpp.write("#endif\n")

if __name__ == "__main__":
    convert_kernels(os.getcwd())
