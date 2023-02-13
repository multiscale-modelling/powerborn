/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef OCTREECELL_H_
#define OCTREECELL_H_

#include "Typedefs.h"

#include <stdio.h>

namespace powerborn {

class Octree;
class OctreeCell;
class CubicGrid;

class CubicGrid
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	CubicGrid()
	{
		dirx_ << 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0;
		diry_ << 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0;
		dirz_ << 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0;
	}
	virtual ~CubicGrid(){}
protected:
	Float8 dirx_, diry_, dirz_;
};


class OctreeCell {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Float8 x, y, z, s;
	UInt8 dtype;
	void print() const {
		for (int i = 0; i < 8; ++i) {
			printf("%7.3f %7.3f %7.3f %7.3f %2d\n", s[i], x[i], y[i], z[i],
					dtype[i]);
		}
	}
	OctreeCell& getChild(unsigned int i)
	{
		return *children[i];
	}
	const OctreeCell& getChild(unsigned int i) const
	{
		return *children[i];
	}
	void getDataFromChildren() {
		for (unsigned int i = 0; i < 8; ++i) {
			unsigned int type = dtype[i];
			if (type == 16)
			{
				OctreeCell& child = getChild(i);
				unsigned int sum = child.dtype.sum();
				type = sum == 8 ? 1 : 16;
				type = sum == 0 ? 0 : type;
				dtype[i] = type;
			}
			if(type == 0) {
				s[i] = 0.0f;
			}
			else if(type == 1) {
				s[i] = (8.0f* s[i]) * (s[i] * s[i]);
			}
			else {
				const OctreeCell& child = getChild(i);
				Float4 sv = *(Float4*) &child.s[0];
				Float4 sx = *(Float4*) &child.x[0] * sv;
				Float4 sy = *(Float4*) &child.y[0] * sv;
				Float4 sz = *(Float4*) &child.z[0] * sv;
				Float4 tv = *(Float4*) &child.s[4];
				sv += tv;
				sx +=*(Float4*) &child.x[4] * tv;
				sy +=*(Float4*) &child.y[4] * tv;
				sz +=*(Float4*) &child.z[4] * tv;
				sv[0] = sv.sum();
				sx[0] = sx.sum();
				sy[0] = sy.sum();
				sz[0] = sz.sum();
				s[i] = sv[0];
				float ivol = 1.0f/sv[0];
				x[i] = sx[0] * ivol;
				y[i] = sy[0] * ivol;
				z[i] = sz[0] * ivol;
			}
		}
	}
	OctreeCell()
	{
	}
	void setChild(OctreeCell& child, unsigned int i)
	{
		children[i] = &child;
	}
private:
	OctreeCell* children[8];
};

}

#endif /* OCTREECELL_H_ */
