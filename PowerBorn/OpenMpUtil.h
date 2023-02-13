/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */
 
// This file is based on code from
// http://bisqwit.iki.fi/story/howto/openmp/

#ifndef OPENMPUTIL_H_
#define OPENMPUTIL_H_

#ifdef _OPENMP
# include <omp.h>
#endif

#include <string.h>
#include <iostream>
#include <ctime>

namespace openmp
{


#ifdef _OPENMP
static double getTime()
{
    return omp_get_wtime();
}

class Mutex
{
public:
	Mutex() { omp_init_lock(&lock); }
	~Mutex() { omp_destroy_lock(&lock); }
	void Lock() { omp_set_lock(&lock); }
	void Unlock() { omp_unset_lock(&lock); }

	Mutex(const Mutex& ) { omp_init_lock(&lock); }
	Mutex& operator= (const Mutex& ) { return *this; }
private:
	omp_lock_t lock;
};

#else

static double getTime()
{
    const double factor = 1.0 / static_cast<double>(CLOCKS_PER_SEC);
    return static_cast<double>(clock()) * factor;
}

/* A dummy mutex that doesn't actually exclude anything,
 * but as there is no parallelism either, no worries. */
class Mutex
{
public:
	void Lock() {}
	void Unlock() {}
};

#endif

class ScopedTimer
{
    double time_;
    std::string name_;

public:
    ScopedTimer(const std::string& name=""): time_(getTime()), name_(name)
    {
    }
    virtual ~ScopedTimer()
    {
        time_ = getTime() - time_;
        std::cout << "Timer " << name_ << ": " << time_ << " seconds." << std::endl;
    }
};

class Timer
{
    double time_, tmp_;
    std::string name_;

public:
    void setName(const std::string& name)
    {
        name_ = name;
    }
    void start()
    {
        tmp_ = getTime();
    }
    void stop()
    {
        time_ += getTime() - tmp_;
    }
    void reset()
    {
        time_ = 0.0;
    }
    double time()
    {
        return time_;
    }
    void print()
    {
        std::cout << "Timer " << name_ << ": " << time_ << std::endl;
    }
    Timer(const std::string& name=""): time_(0.0), tmp_(0.0), name_(name)
    {
    }
};

/* An exception-safe scoped lock-keeper. */
class ScopedLock
{
public:
	explicit ScopedLock(Mutex& m) : mut(m), locked(true) { mut.Lock(); }
	~ScopedLock() { Unlock(); }
	void Unlock() { if(!locked) return; locked=false; mut.Unlock(); }
	void LockAgain() { if(locked) return; mut.Lock(); locked=true; }
private:
	Mutex& mut;
	bool locked;

	// prevent copying the scoped lock.
	void operator=(const ScopedLock&);
	ScopedLock(const ScopedLock&);
};

class OpenMpInterface
{
protected:
	inline int getThreadId() const
	{
#ifdef _OPENMP
		return omp_get_thread_num();
#else
		return 0;
#endif
	}
	inline int getMaxThreads()
	{
#ifdef _OPENMP
		return omp_get_max_threads();
#else
		return 1;
#endif
	}
	virtual ~OpenMpInterface(){}
};

template<class T>
class ThreadPrivateObject: public OpenMpInterface
{
public:
	T& get()
	{
		return *objects_[this->getThreadId()];
	}
	const T& get() const
	{
		return *objects_[this->getThreadId()];
	}
	void setAll(const T& t)
	{
		for(unsigned int i=0; i<objects_.size(); ++i)
		{
			*(objects_[i]) = t;
		}
	}
	unsigned int size() const
	{
		return objects_.size();
	}
	T& getFromThread(int thread_id)
	{
		return *objects_[thread_id];
	}
	const T& getFromThread(int thread_id) const
	{
		return *objects_[thread_id];
	}
	ThreadPrivateObject()
	{
		unsigned int has_exception = 0;
	    objects_.resize(this->getMaxThreads(), NULL);
#pragma omp parallel
		{
	        try
	        {
	            objects_[this->getThreadId()] = new T(); // this should take care of local memory allocation!
	        }
	        catch(...)
	        {
#pragma omp atomic
	            has_exception++;
	        }
		}
		if(has_exception)
		{
		    std::cout << "Exception in ThreadPrivateObject ctor " << std::endl;
		    abort();
		}
	}
	ThreadPrivateObject(const ThreadPrivateObject<T>& other)
	{
	    *this = other;
	}
    ThreadPrivateObject(const T& object)
    {
        unsigned int has_exception = 0;
        objects_.resize(this->getMaxThreads(), NULL);
#pragma omp parallel
        {
            try
            {
                objects_[this->getThreadId()] = new T(object); // this should take care of local memory allocation!
            }
            catch(...)
            {
#pragma omp atomic
                has_exception++;
            }
        }
        if(has_exception)
        {
            std::cout << "Exception in ThreadPrivateObject copy ctor " << std::endl;
            abort();
        }
    }
	ThreadPrivateObject<T>& operator=(const ThreadPrivateObject<T>& other)
	{
		if(this == &other) return *this;
        for(unsigned int i=0; i<objects_.size(); ++i)
        {
            if(objects_[i])
            {
                delete objects_[i];
                objects_[i] = NULL;
            }
        }
        unsigned int has_exception = 0;
		objects_.resize(this->getMaxThreads(), NULL);
#pragma omp parallel
		{
            try
            {
                objects_[this->getThreadId()] = new T(*(other.objects_[this->getThreadId()])); // this should take care of local memory allocation!
            }
            catch(...)
            {
#pragma omp atomic
                has_exception++;
            }
        }
        if(has_exception)
        {
            std::cout << "Exception in ThreadPrivateObject assignment op " << std::endl;
            abort();
        }
		return *this;
	}
	~ThreadPrivateObject()
	{
		for(size_t i=0; i<objects_.size(); ++i)
		{
			if(objects_[i])
            {
                delete objects_[i];
                objects_[i] = NULL;
            }
		}
	}
private:
    std::vector<T*> objects_;
};

} // end of namespace openmp

#endif /* OPENMPUTIL_H_ */
