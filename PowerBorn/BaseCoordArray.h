/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef BASECOORDARRAY_H_
#define BASECOORDARRAY_H_

#include "Typedefs.h"
#include <string>
#include <sstream>
#include <fstream>
#include <cstring>

namespace powerborn {

class BaseCoordArray
{

protected:
    Array x_, y_, z_, r_;
    unsigned int size_, capacity_, p_size_;
    void fixPadding(unsigned int test)
    {
        size_ = test;
        p_size_ = this->round4(test);
        //printf("padding: %d %d\n",test, p_size_);
        for(unsigned int a=size_; a<p_size_; ++a)
        {
            x_[a]=float(0.0);
            y_[a]=float(0.0);
            z_[a]=float(0.0);
            r_[a]=float(0.0);
        }
    }
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    void fixPadding()
    {
        this->fixPadding(size_);
    }
    unsigned int round4(unsigned int n) const
    {
        unsigned int t = (3 + n) & (~3);
        //printf("input: %d   rounded: %d\n",n,t);
        return t;
    }
    Float4 boundingBox() const
    {
        float xmax = this->x().maxCoeff();
        float xmin = this->x().minCoeff();
        float ymax = this->y().maxCoeff();
        float ymin = this->y().minCoeff();
        float zmax = this->z().maxCoeff();
        float zmin = this->z().minCoeff();
        float rmax = this->r().maxCoeff();

        float size = xmax - xmin;
        size = std::max(size, (ymax - ymin));
        size = std::max(size, (zmax - zmin));
        size += 2.0f * rmax;
        Float4 size_center(size, (xmax + xmin), (ymax + ymin), (zmax + zmin));
        return size_center * 0.5f;
    }
    void reserve(const unsigned int s)
    {
        if(s <= capacity_) return;
        capacity_=this->round4(s);
        x_.conservativeResize(capacity_);  // make sure there is enough space to always have multiple of 4 alignment
        y_.conservativeResize(capacity_);  // make sure there is enough space to always have multiple of 4 alignment
        z_.conservativeResize(capacity_);  // make sure there is enough space to always have multiple of 4 alignment
        r_.conservativeResize(capacity_);  // make sure there is enough space to always have multiple of 4 alignment
    }
    void resize(const unsigned int r)
    {
        size_= 0;
        p_size_= 0;
        if(r <= capacity_) return;
        capacity_ = this->round4(r);
        x_.resize(capacity_);  // make sure there is enough space to always have multiple of 4 alignment
        y_.resize(capacity_);  // make sure there is enough space to always have multiple of 4 alignment
        z_.resize(capacity_);  // make sure there is enough space to always have multiple of 4 alignment
        r_.resize(capacity_);  // make sure there is enough space to always have multiple of 4 alignment
    }
    void setSize(unsigned int s)
    {
        size_= s;
        p_size_ = this->round4(s);
        if(s <= capacity_) return;
        capacity_ = this->round4(s);
        x_.resize(capacity_);  // make sure there is enough space to always have multiple of 4 alignment
        y_.resize(capacity_);  // make sure there is enough space to always have multiple of 4 alignment
        z_.resize(capacity_);  // make sure there is enough space to always have multiple of 4 alignment
        r_.resize(capacity_);  // make sure there is enough space to always have multiple of 4 alignment
        this->fixPadding();
    }
    void clear()
    {
        size_ = 0;
        p_size_ = 0;
    }
    unsigned int size() const
    {
        return size_;
    }
    unsigned int pSize() const
    {
        return p_size_;
    }
    void insert(float x, float y, float z)
    {
        if(capacity_ == size_) this -> reserve(size_ + 128);
        x_[size_]=x;
        y_[size_]=y;
        z_[size_]=z;
        size_++;
        p_size_ = this->round4(size_);
    }
    void insert(float x, float y, float z, float r)
    {
        if(capacity_ == size_) this -> reserve(size_ + 128);
        x_[size_]=x;
        y_[size_]=y;
        z_[size_]=z;
        r_[size_]=r;
        size_++;
        p_size_ = this->round4(size_);
    }
    void append(const BaseCoordArray& other)
    {
        this->reserve(size_ + other.size_);
        for(unsigned int i=0; i< other.size_; ++i)
        {
            x_[size_+i] = other.x_[i];
            y_[size_+i] = other.y_[i];
            z_[size_+i] = other.z_[i];
            r_[size_+i] = other.r_[i];
        }
        this->fixPadding(size_ + other.size_);
    }
    void random(const unsigned int s)
    {
        size_=s;
        p_size_=round4(s);
        capacity_=p_size_;
        x_=Eigen::Array<float,Eigen::Dynamic,1>::Random(p_size_)*10.0f;
        y_=Eigen::Array<float,Eigen::Dynamic,1>::Random(p_size_)*10.0f;
        z_=Eigen::Array<float,Eigen::Dynamic,1>::Random(p_size_)*10.0f;
        r_=Eigen::Array<float,Eigen::Dynamic,1>::Random(p_size_).abs()*2.0f + 1.3f;
        this->fixPadding(s);
    }
    void translate(float x, float y, float z)
    {
        x_.head(p_size_) += x;
        y_.head(p_size_) += y;
        z_.head(p_size_) += z;
    }
    void getScaledCoords(const BaseCoordArray& other, float scale)
    {
        this->resize(other.size_);
        for(unsigned int i=0; i<other.size_; ++i)
        {
            x_[i] = other.x_[i] * scale;
            y_[i] = other.y_[i] * scale;
            z_[i] = other.z_[i] * scale;
            r_[i] = other.r_[i];
        }
        this->fixPadding(other.size_);
    }
    void appendTransformed(const BaseCoordArray& other, const Coord3& trans,
            float scale, float r_scale=1.0f)
    {
        unsigned int old_size = size_;
        unsigned int new_size = size_ + other.size_;
        float p0 = trans[0];
        float p1 = trans[1];
        float p2 = trans[2];
        this->reserve(new_size);
        for(unsigned int i=0; i<other.size_; ++i)
        {
            x_[old_size + i] = (other.x_[i] - p0) * scale;
            y_[old_size + i] = (other.y_[i] - p1) * scale;
            z_[old_size + i] = (other.z_[i] - p2) * scale;
            r_[old_size + i] = other.r_[i] * r_scale;
        }
        this->fixPadding(new_size);
    }
    void setRad(float r)
    {
        this->r().setConstant(r);
        for(unsigned int a=size_; a<p_size_; ++a)
        {
            r_[a] = float(0.0f);
        }
    }
    void writePqr(const std::string& outfile, float probe_radius = 1.4f) const
    {
        FILE * pFile;
        pFile = fopen (outfile.c_str(),"w");
        for(unsigned int a=0; a < size_; ++a)
        {
            float x=x_[a];
            float y=y_[a];
            float z=z_[a];
            float r=r_[a] - probe_radius;
            fprintf(pFile,"%6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n","ATOM  ",a+1,"ZN","WAT","A",1,x,y,z,0.0,r);
        }
        fclose(pFile);
        std::cout << "Wrote atoms to pqr file " << outfile << std::endl;
    }
    BaseCoordArray(): size_(0), capacity_(0), p_size_(0)
    {
    }
    BaseCoordArray(unsigned int s): size_(0), p_size_(0)
    {
        this->reserve(s);
    }
    virtual ~BaseCoordArray()
    {
    }
    Array::SegmentReturnType x()
    {
        return x_.head(size_);
    }
    Array::SegmentReturnType y()
    {
        return y_.head(size_);
    }
    Array::SegmentReturnType z()
    {
        return z_.head(size_);
    }
    Array::SegmentReturnType r()
    {
        return r_.head(size_);
    }
    Array::ConstSegmentReturnType x() const
    {
        return x_.head(size_);
    }
    Array::ConstSegmentReturnType y() const
    {
        return y_.head(size_);
    }
    Array::ConstSegmentReturnType z() const
    {
        return z_.head(size_);
    }
    Array::ConstSegmentReturnType r() const
    {
        return r_.head(size_);
    }
    const Coord3 getPos(unsigned int i) const
    {
        return Coord3(x_[i], y_[i], z_[i]);
    }
    void randomRotate()
    {
        Coord3 angle = Coord3::Random();
        const float twopi = 2.0f * std::acos(-1.0);
        angle *= twopi;
        Eigen::Matrix3f m;
        m = Eigen::AngleAxisf(angle[0], Eigen::Vector3f::UnitZ())
            * Eigen::AngleAxisf(angle[1], Eigen::Vector3f::UnitY())
            * Eigen::AngleAxisf(angle[2], Eigen::Vector3f::UnitZ());
        for(unsigned int i=0; i<size_; ++i)
        {
            float tx = x_[i];
            float ty = y_[i];
            float tz = z_[i];
            x_[i] = tx * m(0,0) + ty * m(0,1) + tz * m(0,2);
            y_[i] = tx * m(1,0) + ty * m(1,1) + tz * m(1,2);
            z_[i] = tx * m(2,0) + ty * m(2,1) + tz * m(2,2);
        }
    }
    void getIds(const BaseCoordArray& other, const AlignedVector<unsigned int>::Type& ids)
    {
        unsigned int size = ids.size();
        this->setSize(size);
        for(unsigned int i=0; i<size; ++i)
        {
            unsigned int id = ids[i];
            x_[i] = other.x_[id];
            y_[i] = other.y_[id];
            z_[i] = other.z_[id];
            r_[i] = other.r_[id];
        }
    }
    void appendIds(const BaseCoordArray& other, const AlignedVector<unsigned int>::Type& ids)
    {
        unsigned int size = ids.size();
        unsigned int os = size_;
        this->setSize(os + size);
        for(unsigned int i=0; i<size; ++i)
        {
            unsigned int id = ids[i];
            x_[os + i] = other.x_[id];
            y_[os + i] = other.y_[id];
            z_[os + i] = other.z_[id];
            r_[os + i] = other.r_[id];
        }
    }
    void getIds(const BaseCoordArray& other, const AlignedVector<unsigned int>::Type& ids, unsigned int myid)
    {
        unsigned int size = ids.size();
        this->resize(size);
        unsigned int count = 0;
        for(unsigned int i=0; i<size; ++i)
        {
            unsigned int id = ids[i];
            x_[count] = other.x_[id];
            y_[count] = other.y_[id];
            z_[count] = other.z_[id];
            r_[count] = other.r_[id];
            count += id != myid;
        }
        this->fixPadding(count);
    }
    void parsePqr(const std::string& fname, Array* charges = 0, float probe_radius = 1.4f)
    {
        this->clear();
        std::string atomstr("ATOM");
        std::string temp;
        std::ifstream input(fname.c_str());
        std::vector<float> c;
        while (getline(input,temp))
        {
            try
            {
                if(temp.size() >= 4 && temp.substr(0,4) == atomstr)
                {
                    std::stringstream instream(temp);
                    std::string line=instream.str();
                    std::stringstream data(line.substr(30,line.size()-30));
                    float x,y,z,q,r;
                    data >> x >> y >> z >> q >> r;
                    this->insert(x, y, z, r + probe_radius);
                    c.push_back(q);
                }
            }
            catch(std::out_of_range& e)
            {
                printf("wrong line:  %s",temp.c_str());
            }
        }
        this->fixPadding();
        if(charges)
        {
            charges->resize(c.size());
            for(unsigned int i=0; i<c.size(); ++i)
            {
                (*charges)[i] = c[i];
            }
        }
        if(this->size() == 0)
        {
            std::cout << "no atoms found in file: "<< fname << std::endl;
        }
    }
    void dump(const std::string& outfile) const
    {
        std::ofstream out;
        out.open(outfile.c_str());
        unsigned int *xi, *yi, *zi, *ri;
        if(size_ > 0)
        {
            xi = (unsigned int*) &x_[0];
            yi = (unsigned int*) &y_[0];
            zi = (unsigned int*) &z_[0];
            ri = (unsigned int*) &r_[0];
            for(unsigned int a=0; a < size_; ++a)
            {
                out << xi[a] << " " << yi[a] << " " << zi[a] << " " << ri[a] << std::endl;
            }
        }
        out.close();
    }
    void parseDump(const std::string& fname)
    {
        std::string temp;
        std::ifstream input(fname.c_str());
        uint32_t xi, yi, zi, ri;
        float *x, *y, *z, *r;
        x = (float*) &xi;
        y = (float*) &yi;
        z = (float*) &zi;
        r = (float*) &ri;
        while (getline(input,temp))
        {
            std::stringstream instream(temp);
            instream >> xi >> yi >> zi >> ri;
            this->insert(*x, *y, *z, *r);
        }
        if(this->size() == 0)
        {
            std::cout << "no atoms found in dump: "<< fname << std::endl;
        }
    }
    bool operator==(const BaseCoordArray& other) const
    {
        if(size_ != other.size_) return false;
        if(size_ == 0) return true;
        bool equal = (this->x() == other.x()).all();
        equal &= (this->y() == other.y()).all();
        equal &= (this->z() == other.z()).all();
        equal &= (this->r() == other.r()).all();
        return equal;
    }
};

} // end of namespace

#endif /* BASECOORDARRAY_H_ */
