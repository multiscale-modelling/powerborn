/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef OBJECTPOOL_H_
#define OBJECTPOOL_H_

#include "Eigen/StdDeque"
#include "Eigen/StdVector"
#include "Eigen/Core"

namespace powerborn
{

template <class T, unsigned int S=8>
class BaseStorage
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    typedef Eigen::Array<T, S, 1> Storage;

    unsigned int size() const {return size_;}
    void increaseStorage()
    {
        storage_.push_back(new Storage());
        size_ += S;
    }
    T& operator[](unsigned int i){return (*storage_[i/S])[i%S];}
    const T& operator[](unsigned int i) const {return (*storage_[i/S])[i%S];}
    BaseStorage& operator=(const BaseStorage& other)
    {
        if(this==&other) return *this;
        for(unsigned int i=0; i<size_; ++i)
        {
            if(storage_[i]) delete storage_[i];
        }
        storage_.clear();
        for(unsigned int i=0; i<other.storage_.size(); ++i)
        {
            storage_.push_back(new Storage(*(other.storage_[i])));
        }
        size_ = other.size_;
        return *this;
    }
    BaseStorage(const BaseStorage& other): size_(0)
    {
        *this = other;
    }
    BaseStorage(): size_(0) {}
    ~BaseStorage()
    {
        for(unsigned int i=0; i<storage_.size(); ++i)
        {
            if(storage_[i]) delete storage_[i];
        }
    }
private:
    typename AlignedVector<Storage*>::Type storage_;
    unsigned int size_;
};

template<class T, unsigned int S=64>
class Pool
{
	BaseStorage<T,S> pool_;
	unsigned int counter_;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	T& getObject()
	{
		if(counter_ >= pool_.size())
		{
			pool_.increaseStorage();
		}
		T& obj = pool_[counter_];
		counter_++;
		return obj;
	}
	void reset()
	{
		counter_ = 0;
	}
	inline T* getObjectPtr()
	{
		return &(this->getObject());
	}
	Pool(): counter_(0)
	{
	}
};

template<class T>
class ReusablePool
{
	std::deque<T, Eigen::aligned_allocator<T> > pool_;
	std::deque<T*,  Eigen::aligned_allocator<T*> > unused_;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    // TODO copy ctor and assignment ?
	T& getObject()
	{
		if (unused_.size())
		{
			T& obj = *(unused_.back());
			unused_.pop_back();
			return obj;
		}
		else
		{
			pool_.push_back(T());
			return pool_.back();
		}
	}
	T* getObjectPtr()
	{
		return &(this->getObject());
	}
	void releaseObject(T* obj)
	{
		unused_.push_back(obj);
	}
    void releaseObject(T& obj)
    {
        this->releaseObject(&obj);
    }
    ReusablePool()
    {
    }
	ReusablePool& operator=(const ReusablePool& other)
    {
        // do not copy any data! This is just a pool that provides access to uninitialized objects!!!
        if(this == &other) return *this;
        unused_.clear();
        pool_.clear();
        return *this;
    }
    ReusablePool(const ReusablePool& other)
    {
        *this = other;
    }

};

} // end of namespace

#endif /* OBJECTPOOL_H_ */
