/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#include <PowerBorn/gpu/Task.h>
#include <PowerBorn/BaseCoordArray.h>

#define PRINTVARF(X) (printf("HOST %s: %f\n", #X, X));
#define PRINTVARI(X) (printf("HOST %s: %d\n", #X, X));

namespace powerbornGpu
{

Task::Task(): task_count(0), water_ratio(5)
{
    resolution = 0.25f;
    cl_mem_flags f = CL_MEM_READ_ONLY;
    atoms.setFlags(f);
    f = CL_MEM_READ_WRITE;
    constants.setFlags(f);
}

Task::~Task()
{
}

void Task::printTask()
{
    unsigned int i=0;
    while(i < task_count)
    {
        printf("Constants: i=%d\n",i);
        PRINTVARF(constants[i].root.x)
        PRINTVARF(constants[i].root.y)
        PRINTVARF(constants[i].root.z)
        PRINTVARF(constants[i].root.w)
        PRINTVARI(constants[i].atomCnt)
        PRINTVARI(constants[i].atomPitch)
        PRINTVARI(constants[i].waterCnt)
        PRINTVARI(constants[i].atomsOffset)
        PRINTVARI(constants[i].waterOffset)
        PRINTVARI(constants[i].stackOffset)
        PRINTVARI(constants[i].stackSize)
        PRINTVARI(constants[i].octTreeMaxSize)
        PRINTVARF(constants[i].integration_factor)
        PRINTVARI(constants[i].maxLevel)
        PRINTVARI(constants[i].maincelllevel)
        PRINTVARI(constants[i].maincellcount)
        i++;
    }
}

unsigned int Task::paddedSize(unsigned int s, unsigned int p)
{
    unsigned int p1 = p > 0 ? p - 1 : 0;
    return (s + p1) & (~p1);
}

unsigned int Task::stackSize(unsigned int s, unsigned int water_ratio)
{
    unsigned int psize = Task::paddedSize(s);
    return 16 * (psize + psize * water_ratio) + 64 * psize;
}

void Task::resetTasks()
{
    task_count = 0;
}

void Task::setWaterRatio(unsigned int w)
{
    if(task_count)
    {
        throw SetWaterRatioException();
    }
    water_ratio = w;
}

void Task::initConstants(unsigned int s)
{
    constants.resize(task_count + 1);
    pGPUConstants& c = constants.getHost().back();
    if(task_count > 0)
    {
        c = constants[0];
        c.atomsOffset = task_count * constants[0].atomPitch;
        c.waterOffset = task_count * constants[0].atomPitch * water_ratio;
        c.stackOffset = task_count * this->stackSize(constants[0].atomCnt, water_ratio);
    }
    else
    {
        c.root.x = 0.0f;
        c.root.y = 0.0f;
        c.root.z = 0.0f;
        c.root.w = resolution;

        c.atomCnt = s;
        c.atomPitch = this->paddedSize(s);
        c.atomsOffset = 0;
        c.waterOffset = 0;
        c.stackOffset = 0;
        c.octTreeMaxSize = 0;
        c.integration_factor = 10.0f;

        c.waterCnt = 0; // set in water_kernel
        c.stackSize = 0; // set in water_kernel
        c.maxLevel = 0; // set in water_kernel
        c.maincelllevel = 0; // set in water_kernel
        c.maincellcount = 0; // set in water_kernel
    }
}

void Task::resetConstants()
{
    for(unsigned int i=0; i<task_count; ++i)
    {
        pGPUConstants& c = constants[i];
        c = constants[0];

        c.atomsOffset = i * constants[0].atomPitch;
        c.waterOffset = i * constants[0].atomPitch * water_ratio;
        c.stackOffset = i * this->stackSize(constants[0].atomCnt, water_ratio);

        c.root.x = 0.0f;
        c.root.y = 0.0f;
        c.root.z = 0.0f;
        c.root.w = resolution;

        c.octTreeMaxSize = 0;
        c.integration_factor = 10.0f;

        c.waterCnt = 0; // set in water_kernel
        c.stackSize = 0; // set in water_kernel
        c.maxLevel = 0; // set in water_kernel
        c.maincelllevel = 0; // set in water_kernel
        c.maincellcount = 0; // set in water_kernel
    }
}

void Task::initAtoms(const powerborn::BaseCoordArray& a)
{
    unsigned int padded_size = this->paddedSize(a.size());
    unsigned int offset = task_count * padded_size;
    atoms.resize(padded_size + offset);
    for(unsigned int i=0; i<a.size(); ++i)
    {
        atoms[i + offset].x = a.x()[i];
        atoms[i + offset].y = a.y()[i];
        atoms[i + offset].z = a.z()[i];
        atoms[i + offset].w = a.r()[i];
    }
}

void Task::addToTask(const powerborn::BaseCoordArray& atoms)
{
    if(task_count > 0 && atoms.size() != constants[0].atomCnt)
    {
        throw TaskAtomSizeException();
    }
    this->initConstants(atoms.size());
    this->initAtoms(atoms);
    task_count++;
}

unsigned int Task::taskCount()
{
    return task_count;
}

unsigned int Task::getGlobalSize()
{
    return 256u * task_count;
}

unsigned int Task::getLocalSize()
{
    return 256u;
}

unsigned int Task::getIntegrationBlocks()
{
    return 8u;
}

unsigned int Task::getIntegrationBlockScale()
{
    int ib = getIntegrationBlocks();
    int scale = (64 + ib - 1) / ib;
    scale = std::max(scale, 19); //19 is minimal threshold
    return scale;
}


} /* namespace powerbornGpu */
