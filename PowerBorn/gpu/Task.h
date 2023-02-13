/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef TASK_H_
#define TASK_H_

#include <PowerBorn/gpu/Buffer.h>
#include <PowerBorn/gpu/pGPU_defs.h>

namespace powerborn
{
class BaseCoordArray;
}

namespace powerbornGpu
{

class TaskAtomSizeException: public std::exception
{
public:
    virtual const char* what() const throw()
    {
        return "Tried to atom group to task with wrong group size!";
    }
};

class SetWaterRatioException: public std::exception
{
public:
    virtual const char* what() const throw()
    {
        return "Tried to set water ratio in task with non-zero task_count! Use Task::reset() first!";
    }
};

/*
 * Data and Memory for each evaluation, host and device.
 */
class Task
{
protected:
    VectorBuffer<pGPUConstants> constants;
    VectorBuffer<pGPUPoint4d> atoms;
    float resolution;
    unsigned int task_count;
    unsigned int water_ratio;

    void printTask();

    void initConstants(unsigned int s);
    void initAtoms(const powerborn::BaseCoordArray& atoms);
    void resetConstants();

    static unsigned int paddedSize(unsigned int s, unsigned int padding=32);
    static unsigned int stackSize(unsigned int s, unsigned int water_ratio);
public:
    Task();
    virtual ~Task();
    void resetTasks();
    void setWaterRatio(unsigned int w);
    void addToTask(const powerborn::BaseCoordArray& atoms);
    unsigned int taskCount();
    virtual unsigned int getGlobalSize();
    virtual unsigned int getLocalSize();
    virtual unsigned int getIntegrationBlocks();
    virtual unsigned int getIntegrationBlockScale();
};

} /* namespace powerbornGpu */

#endif /* TASK_H_ */
