/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef TYPEDEFS_H_
#define TYPEDEFS_H_

#include "Eigen/StdVector"
#include "Eigen/StdDeque"
#include "Eigen/Core"
#include "Eigen/Geometry"

#include <stdint.h>
#include <exception>
#include <stdexcept>
#include <iostream>

namespace powerborn {

#ifndef NULL
#define NULL 0;
#endif

typedef Eigen::ArrayXf Array;
typedef Eigen::Array<int32_t, Eigen::Dynamic, 1> ArrayI;
typedef Eigen::Array<uint32_t, Eigen::Dynamic, 1> ArrayU;
typedef Eigen::Array<unsigned char, Eigen::Dynamic, 1> ArrayUc;
typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayB;

typedef Eigen::Array<float, 8, 1> Float8;
typedef Eigen::Matrix<float, 8, 1> Vector8;

typedef Eigen::Array<uint32_t, 8, 1> UInt8;
typedef Eigen::Array<uint32_t, 4, 1> UInt4;
typedef Eigen::Array<uint32_t, 2, 1> UInt2;
typedef Eigen::Array<int32_t, 4, 1> Int4;
typedef Eigen::Array<int32_t, 2, 1> Int2;
typedef Eigen::Array4f Float4;
typedef Eigen::Vector4f Vector4;
typedef Eigen::Vector3f Coord3;
typedef Eigen::Vector3f Coord;
typedef Eigen::Vector3d Coord3d;
typedef Eigen::Array<Float4, Eigen::Dynamic, 1> Coords;
typedef Eigen::Array<Float4, 8, 1> Coord8;
typedef Eigen::Array<float, 32, 1> Float32;
typedef Eigen::Array<double, 24, 1> Double24;

template<class T, int Size>
struct AlignedArray {
    typedef Eigen::Array<T, Size, 1> Type;
};

template<class T>
struct AlignedVector {
    typedef std::vector<T, Eigen::aligned_allocator<T> > Type;
private:
    AlignedVector() {}
};

template<class T>
struct AlignedDeque {
    typedef std::deque<T, Eigen::aligned_allocator<T> > Type;
};

template<class T>
void prefetch_nta(T* p)
{
#ifdef __SSE__
    _mm_prefetch((const char*) p, _MM_HINT_NTA);
#endif
}

template<class T>
void prefetch(T* p)
{
#ifdef __SSE__
    _mm_prefetch((const char*) p, _MM_HINT_T0);
#endif
}

} // end of namespace

#endif /* TYPEDEFS_H_ */
