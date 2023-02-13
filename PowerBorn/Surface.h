/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef SURFACE_H_
#define SURFACE_H_

#include "Typedefs.h"
#include "Octree.h"
#include "Atoms.h"
#include "ObjectPool.h"

namespace powerborn
{

class Surface : public Octree
{

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Surface(float resolution=0.25f): Octree(resolution)
	{
	}
	virtual ~Surface() {}
	void init(const CoordArray& atoms)
	{
#ifdef _OPENMP

		this->initOmp(atoms);
#else
		this->initSerial(atoms);
#endif
	}
	void build(const CoordArray& atoms, const CoordArray& waters)
	{
#ifdef _OPENMP

		this->buildOmp(atoms,waters);
#else
		this->buildSerial(atoms,waters);
#endif
	}
    void finalize()
    {
#ifdef _OPENMP
#pragma omp single nowait
        {
            for(unsigned int i=0; i<64; ++i)
            {
                level2_[i]->getDataFromChildren();
            }
            for(unsigned int i=0; i<8; ++i)
            {
                for(unsigned int j=0; j<8; ++j)
                {
                    level1_[i]->setChild(*level2_[i*8+j], j); // dtype is already set for level1 cells!
                }
                level1_[i]->getDataFromChildren();
            }
            root_->getDataFromChildren();
        }
#endif
    }

private:
	openmp::ThreadPrivateObject<ReusablePool<CoordArray> > coord_array_pools_;
	AlignedArray<OctreeCell*, 8>::Type level1_;
	AlignedArray<OctreeCell*, 64>::Type level2_;


	void initSerial(const CoordArray& atoms)
	{
		this->initSizeCenter(atoms);
		this->initOctree();
	}
	void initOmp(const CoordArray& atoms)
	{
#pragma omp single nowait
		{
            this->initSizeCenter(atoms);
            this->initOctree();
			Pool<OctreeCell>& cell_pool = cell_pools_.get();
			for(unsigned int i=0; i<8; ++i)
			{
				level1_[i] = cell_pool.getObjectPtr();
				this->constructChildCell(*root_, *level1_[i], i);
				level1_[i]->dtype.setConstant(16); // they will get children in the next step!
			}
			for(unsigned int i=0; i<64; ++i)
			{
				OctreeCell& parent = *level1_[i/8];
				unsigned int j = i%8;
				OctreeCell* child = cell_pool.getObjectPtr();
				level2_[i] = child;
				this->constructChildCell(*child, parent.s[j], parent.x[j], parent.y[j], parent.z[j]);
			}
		}
	}
	void buildSerial(const CoordArray& atoms, const CoordArray& waters)
	{
		this->buildRecursive(*root_, 0, atoms, waters, coord_array_pools_.get(), cell_pools_.get());
	}
	void buildOmp(const CoordArray& atoms, const CoordArray& waters)
	{
		Pool<OctreeCell>& cell_pool = cell_pools_.get();
		ReusablePool<CoordArray>& coord_pool = coord_array_pools_.get();
#pragma omp for schedule(dynamic) nowait
		for(unsigned int i=0; i<512; ++i)
		{
			this->recursiveBuildStep1(*level2_[i/8], 2, atoms, waters, i%8, coord_pool, cell_pool);
		}
	}
	void recursiveBuildStep1(OctreeCell& c, const int level,
			const CoordArray& atoms, const CoordArray& waters, unsigned int i,
			ReusablePool<CoordArray>& coord_pool, Pool<OctreeCell>& cell_pool)
	{
		CoordArray& my_atoms = coord_pool.getObject();
		my_atoms.update(c.x[i], c.y[i], c.z[i], c.s[i], atoms);
		if (my_atoms.size())
		{
			CoordArray& my_waters = coord_pool.getObject();
			my_waters.updateFixedRad(c.x[i], c.y[i], c.z[i], c.s[i], waters);
			unsigned int cell_type = this->getCellType(my_atoms.size(), my_waters.size());
			c.dtype[i] = cell_type;
			if(cell_type == 16)
			{
				OctreeCell& child = cell_pool.getObject();
				this->constructChildCell(c, child, i);
				this->buildRecursive(child, level + 1, my_atoms, my_waters, coord_pool, cell_pool);
			}
			coord_pool.releaseObject(my_waters);
		}
		else
		{
			c.dtype[i] = 1;
		}
		coord_pool.releaseObject(my_atoms);
	}
	void recursiveBuildStep2(OctreeCell& c, const int level,
			const CoordArray& atoms, const CoordArray& waters, unsigned int i,
			ReusablePool<CoordArray>& coord_pool, Pool<OctreeCell>& cell_pool)
	{
		CoordArray& my_atoms = coord_pool.getObject();
		my_atoms.update(c.x[i], c.y[i], c.z[i], 0.75f * c.s[i], atoms);
		if (my_atoms.size())
		{
			CoordArray& my_waters = coord_pool.getObject();
			my_waters.updateFixedRad(c.x[i], c.y[i], c.z[i], 0.75f * c.s[i], waters);
			unsigned int cell_type = this->getCellType(my_atoms.size(), my_waters.size());
			c.dtype[i] = cell_type;
			if(cell_type == 16)
			{
				unsigned int is_buried = this->isBuriedInWater(c, i, my_waters);
				if(! is_buried)
				{
					OctreeCell& child = cell_pool.getObject();
					this->constructChildCell(c, child, i);
					this->buildRecursive(child, level + 1, my_atoms, my_waters, coord_pool, cell_pool);
				}
				else
				{
					c.dtype[i] = 1;
				}
			}
			coord_pool.releaseObject(my_waters);
		}
		else
		{
			c.dtype[i] = 1;
		}
		coord_pool.releaseObject(my_atoms);
	}
	void buildRecursive(OctreeCell& c, const int level,
			const CoordArray& atoms, const CoordArray& waters,
			ReusablePool<CoordArray>& coord_pool, Pool<OctreeCell>& cell_pool)
	{
		if(level < end_level_ - 1)
		{
			for(unsigned int i=0; i<8; ++i)
			{
				this->recursiveBuildStep1(c, level, atoms, waters, i, coord_pool, cell_pool);
			}
		}
		else if(level < end_level_)
		{
			for(unsigned int i=0; i<8; ++i)
			{
				this->recursiveBuildStep2(c, level, atoms, waters, i, coord_pool, cell_pool);
			}
		}
		else
		{
			for(unsigned int i=0; i<8; ++i)
			{
				unsigned int cell_type = this->isBuriedInWater(c, i, waters, 0.5f);
				c.dtype[i] = cell_type;
				if (cell_type == 0)
				{
					this->decideSurface(c, i, atoms, waters, cell_pool);
				}
			}
		}
        prefetch_nta((const char*) &c);
		c.getDataFromChildren();
	}
	unsigned int isBuriedInWater(OctreeCell& c, unsigned int i, const CoordArray& waters, float size_scale = 1.0f)
	{
		return waters.cubeInsideSphere(c.x[i], c.y[i], c.z[i], c.s[i] * size_scale);
	}
	inline void decideSurface(OctreeCell& c, unsigned int blockid, const CoordArray& atoms,
			const CoordArray& waters, Pool<OctreeCell>& cell_pool)
	{
		Float8 cx, cy, cz, cs = Float8::Zero();
		UInt8 dtype;
		float s = c.s[blockid] * 0.5f;
		//printf("smallest size: %f\n", s);
		cx = dirx_ * s + c.x[blockid];
		cy = diry_ * s + c.y[blockid];
		cz = dirz_ * s + c.z[blockid];
		float volume = s * s * s * 8.0f;
		for(unsigned int i=0; i<8; ++i)
		{
			bool cell_type = !atoms.insideSphere(cx[i], cy[i], cz[i]);
			if (!cell_type)
			{
				cell_type = waters.insideSphere(cx[i], cy[i], cz[i]);
			}
			dtype[i] = cell_type; //? 1: 0;
			cs[i] = cell_type? volume : cs[i];
		}
		unsigned int sum = dtype.sum();
		unsigned int type = sum == 8 ? 1 : 16;
		type = sum == 0 ? 0 : type;
		c.dtype[blockid] = type;
		if (type == 16)
		{
			OctreeCell& child = cell_pool.getObject();
			c.setChild(child, blockid);
			child.x = cx;
			child.y = cy;
			child.z = cz;
			child.s = cs;
			child.dtype = dtype;
		}
	}
	inline unsigned int getCellType(unsigned int atomcount, unsigned int watercount)
	{
		unsigned int type = (atomcount != 0 && watercount != 0)? (unsigned int) 16 : (unsigned int) 0;
		type = (atomcount == 0)? (unsigned int)1 : type;
		return type;
	}
};

} // end of namespace

#endif /* SURFACE_H_ */
