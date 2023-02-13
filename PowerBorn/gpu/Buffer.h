/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef BUFFER_H_
#define BUFFER_H_

#include <openclsetup.h>
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>
#include <vector>
#include <iostream>
#include <exception>
#include <cstring>

namespace powerbornGpu
{

class BufferFlagsException: public std::exception
{
public:
    virtual const char* what() const throw()
    {
        return "Tried to enque buffer operation with no OpenCL memory flags set!";
    }
};

class BufferNonExistentException: public std::exception
{
public:
    virtual const char* what() const throw()
    {
        return "Tried to read from a OpenCL buffer that does not exist!";
    }
};

class BufferSizeMismatchException: public std::exception
{
public:
    virtual const char* what() const throw()
    {
        return "Tried to copy between host and device buffer with different sizes!";
    }
};

class BufferSizeZeroException: public std::exception
{
public:
    virtual const char* what() const throw()
    {
        return "Tried to create a buffer with 0 size!";
    }
};

class BufferInvalidHostPtr: public std::exception
{
public:
    virtual const char* what() const throw()
    {
        return "There is no valid Buffer host pointer available";
    }
};

class BufferWriteToDeviceFailure: public std::exception
{
public:
    virtual const char* what() const throw()
    {
        return "Buffer data on device is incorrect";
    }
};

class BaseBuffer: public boost::noncopyable // noncopyable to prevent unintended sharing of cl::Buffers
{
protected:
    boost::shared_ptr<cl::Buffer> devicebuffer;
    cl_mem_flags flags;
    bool flags_set;
    unsigned int size_; // size of memory on the device in bytes

    virtual unsigned int getSize()
    {
        return size_;
    }
public:
    typedef boost::shared_ptr<cl::Buffer> BufferPtr;

    BaseBuffer()
    {
        flags = 0;
        flags_set = false;
        size_ = 0;
    }
    virtual ~BaseBuffer()
    {
    }
    cl::Buffer& getBuffer()
    {
        return *devicebuffer;
    }
    const cl::Buffer& getBuffer() const
    {
        return *devicebuffer;
    }
    void setFlags(cl_mem_flags f)
    {
        flags = f;
        flags_set = true;
    }
};

template<class T>
class DeviceBuffer: public BaseBuffer
{
public:
    void resize(opencl_setup& ocl, unsigned int s)
    {
        if(!s)
        {
            throw BufferSizeZeroException();
        }
        unsigned int tmp = s * sizeof(T);
        if(size_ < tmp)
        {
            devicebuffer.reset(ocl.new_buffer_object(flags, tmp));
            size_ = tmp;
        }
    }
    unsigned int size()
    {
        return this->getSize() / sizeof(T);
    }
};

template<class T>
class Buffer : public BaseBuffer
{
protected:
    T host;
    virtual unsigned int getSize()
    {
        return sizeof(T);
    }
    virtual void* getHostLocation()
    {
        return (void*) &host;
    }
    virtual void checkWriteToDevice(opencl_setup& ocl)
    {
        HostType host_copy = host;
        ocl.getqueue().finish();
        this->readFromDevice(ocl);
        ocl.getqueue().finish();
        int comparison = memcmp(&host_copy, this->getHostLocation(), this->getSize());
        if(comparison != 0)
        {
            throw BufferWriteToDeviceFailure();
        }
    }
public:
    typedef T HostType;
    typedef boost::shared_ptr<cl::Buffer> BufferPtr;
    Buffer(): BaseBuffer()
    {
        size_ = sizeof(T);
    }
    Buffer(const T& t): BaseBuffer()
    {
        host = t;
        size_ = sizeof(T);
    }
    virtual ~Buffer()
    {
    }
    HostType& getHost()
    {
        return host;
    }
    const HostType& getHost() const
    {
        return host;
    }
    virtual void writeToDevice(opencl_setup& ocl, const std::vector<cl::Event>* events = NULL)
    {
        if(!flags_set)
        {
            throw BufferFlagsException();
        }
        if(!this->getSize())
        {
            throw BufferSizeZeroException();
        }
        try
        {
            if(!devicebuffer)
            {
                devicebuffer.reset(ocl.new_buffer_object(flags, this->getSize()));
                size_ = this->getSize();
            }
            if(size_ < this->getSize())
            {
                devicebuffer.reset(ocl.new_buffer_object(flags, this->getSize()));
                size_ = this->getSize();
            }
            ocl.getqueue().enqueueWriteBuffer(*devicebuffer, CL_TRUE, 0, this->getSize(), this->getHostLocation(), events);

            // debug check to make sure data arrived at device
            this->checkWriteToDevice(ocl);
        } catch(cl::Error& err)
        {
            std::cerr << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
            throw err;
        }
    }
    virtual void readFromDevice(opencl_setup& ocl, const std::vector<cl::Event>* events = NULL)
    {
        try
        {
            if(!devicebuffer)
            {
                throw BufferNonExistentException();
            }
            if(size_ < this->getSize())
            {
                throw BufferSizeMismatchException();
            }
            ocl.getqueue().enqueueReadBuffer(*devicebuffer, CL_TRUE, 0, this->getSize(), (void*) this->getHostLocation(), events);
        } catch(cl::Error& err)
        {
            std::cerr << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
            throw err;
        }
    }
};

template<class T>
class VectorBuffer: public Buffer<std::vector<T> >
{
    typedef Buffer<std::vector<T> > BaseType;
    virtual unsigned int getSize()
    {
        // size of host memory
        return sizeof(T) * BaseType::host.size();
    }
    virtual void* getHostLocation()
    {
        if(! (BaseType::host.size()))
        {
            throw BufferInvalidHostPtr();
        }
        return (void*)(&(BaseType::host)[0]);
    }
    virtual void checkWriteToDevice(opencl_setup& ocl)
    {
        typename BaseType::HostType host_copy = BaseType::host;
        ocl.getqueue().finish();
        this->readFromDevice(ocl);
        ocl.getqueue().finish();
        int comparison = memcmp(host_copy.data(), this->getHostLocation(), this->getSize());
        if(comparison != 0)
        {
            throw BufferWriteToDeviceFailure();
        }
    }
public:
    virtual ~VectorBuffer()
    {
    }
    T& operator[](unsigned int i)
    {
        return BaseType::host[i];
    }
    const T& operator[](unsigned int i) const
    {
        return BaseType::host[i];
    }
    VectorBuffer()
    {
        BaseType::flags = 0;
        BaseType::flags_set = false;
        BaseType::size_ = 0;
    }
    VectorBuffer(const typename BaseType::HostType& t)
    {
        BaseType::host = t;
        BaseType::flags = 0;
        BaseType::flags_set = false;
        BaseType::size_ = 0;
    }
    void resize(unsigned int s)
    {
        BaseType::host.resize(s);
    }
    void resize(unsigned int s, const T& t)
    {
        BaseType::host.resize(s, t);
    }
    void reallocDevice(opencl_setup& ocl)
    {
        if(!BaseType::flags_set)
        {
            throw BufferFlagsException();
        }
        if(!this->getSize())
        {
            throw BufferSizeZeroException();
        }
        try
        {
            if(!BaseType::devicebuffer || (BaseType::size_ < this->getSize()))
            {
                BaseType::devicebuffer.reset(ocl.new_buffer_object(BaseType::flags, this->getSize()));
                BaseType::size_ = this->getSize();
            }
        } catch(cl::Error& err)
        {
            std::cerr << err.what() << " " << err.err() << " " << opencl_setup::getErrorString(err.err()) << std::endl;
            throw err;
        }
    }
};

} // end of namespace powerbornGpu



#endif /* BUFFER_H_ */
