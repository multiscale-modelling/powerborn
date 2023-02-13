/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef ATOMS_DEFAULT_H_
#define ATOMS_DEFAULT_H_

#include "BaseCoordArray.h"

namespace powerborn {

class CoveredAtomException: public std::exception
{
};

class CoordArray: public BaseCoordArray
{
protected:
    ArrayI tmp_;
    ArrayU ids1_, ids2_;
    unsigned int true_count_, false_count_;
    void collectTrue(const CoordArray& other)
    {
        unsigned int stop=true_count_;
        for(unsigned int a=0; a < stop; ++a)
        {
            unsigned int id = ids1_[a];
            x_[a]=other.x_[id];
            y_[a]=other.y_[id];
            z_[a]=other.z_[id];
            r_[a]=other.r_[id];
        }
        this->fixPadding(stop);
    }
    void collectFalse(const CoordArray& other)
    {
        unsigned int stop=false_count_;
        for(unsigned int a=0; a < stop; ++a)
        {
            unsigned int id=ids2_[a];
            x_[a]=other.x_[id];
            y_[a]=other.y_[id];
            z_[a]=other.z_[id];
            r_[a]=other.r_[id];
        }
        this->fixPadding(stop);
    }
    void collectFixedRadTrue(const CoordArray& other)
    {
        unsigned int stop=true_count_;
        for(unsigned int a=0; a < stop; ++a)
        {
            unsigned int id = ids1_[a];
            x_[a]=other.x_[id];
            y_[a]=other.y_[id];
            z_[a]=other.z_[id];
        }
        this->fixPadding(stop);
        this->setRad(other.r_[0]);
    }
    void compressTrue()
    {
        if(false_count_==0) return;
        if(true_count_==0)
        {
            this->fixPadding(0);
            return;
        }
        int index1=true_count_-1;
        int stop = false_count_;
        int index2=0;
        while((index1 >= 0) & (index2 < stop)) // TODO check this
        {
            unsigned int empty = ids2_[index2];
            unsigned int full = ids1_[index1];
            if (empty >= full) break;
            x_[empty]=x_[full];
            y_[empty]=y_[full];
            z_[empty]=z_[full];
            r_[empty]=r_[full];
            index1--;
            index2++;
        }
        this->fixPadding(true_count_); // TODO check this!
    }
    void compressFalse()
    {
        if(true_count_==0) return;
        if(false_count_==0)
        {
            this->fixPadding(0);
            return;  // no false values
        }
        int index1 = 0;
        int stop = true_count_;
        int index2 = false_count_ - 1;
        while( (index2 >= 0) & (index1 < stop))
        {
            unsigned int empty=ids2_[index2];
            unsigned int full=ids1_[index1];
            if (full >= empty) break;
            x_[full]=x_[empty];
            y_[full]=y_[empty];
            z_[full]=z_[empty];
            r_[full]=r_[empty];
            index1++;
            index2--;
        }
        this->fixPadding(false_count_);  // TODO check this line!
    }
    void allIndices(unsigned int stop)
    {
        int32_t index1=0;
        int32_t index2=0;
        for(unsigned int a=0; a < stop; ++a)
        {
            ids1_[index1]=a; // truth indices
            ids2_[index2]=a; // false indices
            int32_t tmp= *(int32_t*) &tmp_[a];
            index2 -= ~tmp;
            index1 -= tmp;
        }
        true_count_=index1;  // array length of _ids1
        false_count_=index2;  // array length of _ids2
    }
    void trueIndices(unsigned int stop)
    {
        int32_t index=0;
        for(unsigned int a=0; a < stop; ++a)
        {
            ids1_[index]=a;
            index -= *(int32_t*) &tmp_[a];
        }
        true_count_ = index;
    }
    void falseIndices(unsigned int stop)
    {
        int32_t index=0;
        for(unsigned int a=0; a < stop; ++a)
        {
            ids2_[index]=a;
            index -= ~(*(int32_t*) &tmp_[a]);
        }
        false_count_ = index;
    }
    void computeFixedRad(float x,float y,float z,float thres, const CoordArray& other)
    {
        tmp_.head(other.size_) = -((other.x()-x).square() + (other.y()-y).square() < -(other.z()-z).square() + thres).cast<int32_t>();
    }
    void compute(float x,float y,float z,float s, const CoordArray& other)
    {
        tmp_.head(other.size_) = -((other.x()-x).square() + (other.y()-y).square() < (other.r() + s).square() - (other.z()-z).square()).cast<int32_t>();
    }
    bool computeRealNeighbours(float x,float y,float z,float s, const CoordArray& other)
    {

        bool zero_nb_dist = false;
        unsigned int stop=other.p_size_;
        for (unsigned int a=0; a< stop; a++)
        {
            float dx = other.x_[a] - x;
            float dy = other.y_[a] - y;
            float dz = other.z_[a] - z;
            float dr = other.r_[a] + s;
            dx *= dx;
            dy *= dy;
            dz *= dz;
            dr *= dr;
            dx += dy;
            dx += dz;
            zero_nb_dist |= dx == 0.0f;
            tmp_[a] = -(int32_t( dx < dr));
        }
        return zero_nb_dist;
    }
    void checkPower(float x,float y,float z,float r, const CoordArray& other)
    {
        tmp_.head(other.size_) = -((other.x()-x).square() + (other.y()-y).square() + (other.z()-z).square() < other.r() + (r*r)).cast<int32_t>();
    }
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    ArrayU::SegmentReturnType trueIds() {return ids1_.head(true_count_);}
    ArrayU::ConstSegmentReturnType trueIds() const {return ids1_.head(true_count_);}
    ArrayU::SegmentReturnType falseIds() {return ids2_.head(false_count_);}
    ArrayU::ConstSegmentReturnType falseIds() const {return ids2_.head(false_count_);}
    ArrayI::SegmentReturnType tmp() {return tmp_.head(size_);}
    ArrayI::ConstSegmentReturnType tmp() const {return tmp_.head(size_);}
    unsigned int getRealNeighbours(float x, float y, float z, float r)
    {
        if(tmp_.size() < p_size_)
        {
            tmp_.resize(p_size_);
            ids1_.resize(p_size_);
            ids2_.resize(p_size_);
        }
        unsigned int same_radii_count = 0;
        bool zero_dist = this->computeRealNeighbours(x,y,z,r,*this);
        if(zero_dist)
        {
            //std::cout << "size before duplicate removal " << size_ << " ";
            this->computeFixedRad(x,y,z,0.0f,*this);
            this->allIndices(size_);
            //std::cout << "analysis " << tmp_.sum() << " " << true_count_ << " " << false_count_ << " ";
            bool is_covered = true_count_ > 0;
            for(unsigned int i=0; i<true_count_; ++ i)
            {
                unsigned int nbid = ids1_[i];
                is_covered &=  r <= r_[nbid];
                same_radii_count += (r == r_[nbid]);
            }
            if(is_covered && same_radii_count==0)
            {
                throw CoveredAtomException();
            }
            this->compressFalseOrdered();
            //std::cout << "size after removal of duplicates " << size_ << std::endl;
            this->computeRealNeighbours(x,y,z,r,*this);
        }
        this->trueIndices(size_);
        this->collectTrue(*this);
        return same_radii_count;
    }
    void update(float x,float y,float z,float s, const CoordArray& other)
    {
        if(tmp_.size() < other.p_size_) tmp_.resize(other.p_size_);
        if (ids1_.size() < other.p_size_) ids1_.resize(other.p_size_);
        s *= float(sqrt(float(3.0)));

        this->compute(x,y,z,s,other);
        this->trueIndices(other.size_);
        this->resize(true_count_);
        this->collectTrue(other);
    }
    void updateFixedRad(float x,float y,float z,float s, const CoordArray& other)
    {
        if(other.size() == 0)
        {
            this->clear();
            return;
        }
        if(tmp_.size() < other.p_size_) tmp_.resize(other.p_size_);
        if (ids1_.size() < other.p_size_) ids1_.resize(other.p_size_);
        float thres = sqrt(float(3.0))*s;
        thres+=other.r_[0];
        thres*=thres;
        this->computeFixedRad(x,y,z,thres,other);
        this->trueIndices(other.size_);
        this->resize(true_count_);
        this->collectFixedRadTrue(other);
    }
    void checkSphere(float x, float y, float z, float r)
    {
        if(tmp_.size() < p_size_) tmp_.resize(p_size_);
        if (ids1_.size() < p_size_) ids1_.resize(p_size_);
        if (ids2_.size() < p_size_) ids2_.resize(p_size_);
        this->computeFixedRad(x,y,z,r*r,*this);
        this->allIndices(size_);
        this->compressFalse();
    }
    void checkVertices(float x, float y, float z, float r)
    {
        if(tmp_.size() < p_size_) tmp_.resize(p_size_);
        if (ids1_.size() < p_size_) ids1_.resize(p_size_);
        if (ids2_.size() < p_size_) ids2_.resize(p_size_);
        this->checkPower(x,y,z,r,*this);
        this->allIndices(size_);
        //this->compressFalse();
    }
    bool insideSphere(float xi, float yi, float zi) const
    {
        return ((this->x()-xi).square() + (this->y()-yi).square() < (this->r()).square() - (this->z()-zi).square()).any();
    }
    unsigned int cubeInsideSphere(float xc, float yc, float zc, float sc) const
    {
         sc *= float(sqrt(float(3.0)));
        return ((this->x()-xc).square() + (this->y()-yc).square() < (this->r() - sc).square() - (this->z()-zc).square()).any();
    }
    void compressFalseOrdered()
    {
        if(true_count_==0) return;
        if(false_count_==0)
        {
            this->fixPadding(0);
            return;  // no false values
        }
        unsigned int stop = false_count_;
        for(unsigned int i=0; i<stop; ++i)
        {
            unsigned int id2 = ids2_[i];
            x_[i]=x_[id2];
            y_[i]=y_[id2];
            z_[i]=z_[id2];
            r_[i]=r_[id2];
        }
        this->fixPadding(false_count_);  // TODO check this line!
    }
    void compressTrueOrdered()
    {
        if(true_count_==0) return;
        if(false_count_==0)
        {
            this->fixPadding(0);
            return;  // no false values
        }
        unsigned int stop = true_count_;
        for(unsigned int i=0; i<stop; ++i)
        {
            unsigned int id1 = ids1_[i];
            x_[i]=x_[id1];
            y_[i]=y_[id1];
            z_[i]=z_[id1];
            r_[i]=r_[id1];
        }
        this->fixPadding(true_count_);  // TODO check this line!
    }
    CoordArray& operator=(const CoordArray& other)
    {
        /*
         * For performance reasons does not copy arrays tmp_, ids1_, ids2_, false_count, true_count !!!
         * This is just an optimized base class assignment op, that copies the coords and radii in one loop.
         */
        this->setSize(other.size_);
        for(unsigned int i=0; i<other.size_; ++i)
        {
            x_[i] = other.x_[i];
            y_[i] = other.y_[i];
            z_[i] = other.z_[i];
            r_[i] = other.r_[i];
        }
        return *this;
    }
    CoordArray(const CoordArray& other) = default;
    CoordArray(const BaseCoordArray& base): BaseCoordArray(base)
    {
        this->fixPadding();
    }
    CoordArray()
    {
    }
    CoordArray(unsigned int s)
    {
        this->reserve(s);
    }
    virtual ~CoordArray()
    {
    }
};

} // end of namespace

#endif /* ATOMS_DEFAULT_H_ */
