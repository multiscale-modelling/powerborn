/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef NEIGHBOURGRID_H_
#define NEIGHBOURGRID_H_

#include <stdint.h>
/*
 * the results of the PowerDiagram depend on the order of the atoms being added.
 * Unordered Maps result in different orders of atoms being added, if the nb grid is reconstructed with the same atom
 * therefore results are not binary identical anymore
 * thuse the unordered maps are replaced with std ones
 */
//#include <boost/unordered_map.hpp>
#include <map>

#include "Typedefs.h"
#include "Atoms.h"
#include "OpenMpUtil.h"

namespace powerborn
{

class NeighbourGrid
{
public:
	typedef uint64_t Key;
	typedef Eigen::Array<uint64_t, Eigen::Dynamic, 1> KeyArray;
	typedef Eigen::Array<uint64_t, 27, 1> KeyArray27;
    typedef std::map<uint64_t, unsigned int > Map;
	typedef std::multimap<uint64_t, unsigned int> MultiMap;
	typedef MultiMap::iterator Iterator;
	typedef MultiMap::const_iterator ConstIterator;
	typedef std::pair<Iterator, Iterator> IteratorPair;
	typedef std::pair<ConstIterator, ConstIterator> ConstIteratorPair;
	typedef AlignedVector<unsigned int>::Type IdVec;

private:
	Map key_to_id_;
	MultiMap map_;
	AlignedVector<IdVec>::Type vmap_;
	KeyArray keys_;
	float grid_space_, inv_grid_space_, xmin_, ymin_, zmin_;
	IdVec empty_;

	uint64_t intsToKey(uint64_t xi, uint64_t yi, uint64_t zi) const
	{
		uint64_t r = ( xi << 42) | (yi << 21) | zi;
		return r;
	}
	int64_t pointToKey(float x, float y, float z)
	{
		const uint64_t mask = ((1 << 21) - 1);
		uint64_t xi = (uint64_t(inv_grid_space_ * x) & mask) + 1;
		uint64_t yi = (uint64_t(inv_grid_space_ * y) & mask) + 1;
		uint64_t zi = (uint64_t(inv_grid_space_ * z) & mask) + 1;
		return this->intsToKey(xi, yi, zi);
	}
    void getKeys(KeyArray27& keys, uint64_t start) const
    {
        uint64_t xi, yi, zi;
        const uint64_t mask = ((1u << 21) - 1);
        zi = (mask & start);
        yi = (mask & (start >> 21));
        xi = (mask & (start >> 42));

        keys[0] = this->intsToKey(xi, yi, zi);
        keys[1] = this->intsToKey(xi-1, yi, zi);
        keys[2] = this->intsToKey(xi+1, yi, zi);
        keys[3] = this->intsToKey(xi, yi-1, zi);
        keys[4] = this->intsToKey(xi, yi+1, zi);
        keys[5] = this->intsToKey(xi, yi, zi-1);
        keys[6] = this->intsToKey(xi, yi, zi+1);
        keys[7] = this->intsToKey(xi+1, yi+1, zi);
        keys[8] = this->intsToKey(xi+1, yi-1, zi);
        keys[9] = this->intsToKey(xi-1, yi+1, zi);
        keys[10] = this->intsToKey(xi-1, yi-1, zi);
        keys[11] = this->intsToKey(xi-1, yi, zi+1);
        keys[12] = this->intsToKey(xi+1, yi, zi+1);
        keys[13] = this->intsToKey(xi, yi-1, zi+1);
        keys[14] = this->intsToKey(xi, yi+1, zi+1);
        keys[15] = this->intsToKey(xi-1, yi, zi-1);
        keys[16] = this->intsToKey(xi+1, yi, zi-1);
        keys[17] = this->intsToKey(xi, yi-1, zi-1);
        keys[18] = this->intsToKey(xi, yi+1, zi-1);
        keys[19] = this->intsToKey(xi-1, yi-1, zi-1);
        keys[20] = this->intsToKey(xi-1, yi-1, zi+1);
        keys[21] = this->intsToKey(xi-1, yi+1, zi-1);
        keys[22] = this->intsToKey(xi-1, yi+1, zi+1);
        keys[23] = this->intsToKey(xi+1, yi-1, zi-1);
        keys[24] = this->intsToKey(xi+1, yi-1, zi+1);
        keys[25] = this->intsToKey(xi+1, yi+1, zi-1);
        keys[26] = this->intsToKey(xi+1, yi+1, zi+1);
    }
    void init(const BaseCoordArray& atoms)
    {
#pragma omp single nowait
        if(keys_.size() != atoms.size()) keys_.resize(atoms.size());
#pragma omp single nowait
        {
            grid_space_ = 2.0f * atoms.r().maxCoeff();
            inv_grid_space_ = 1.0f / grid_space_;
        }
#pragma omp single nowait
        xmin_ = atoms.x().minCoeff();
#pragma omp single nowait
        ymin_ = atoms.y().minCoeff();
#pragma omp single nowait
        zmin_ = atoms.z().minCoeff();
#pragma omp barrier
        //std::cout << "Neighbour Grid init " << grid_space_ << " " << xmin_ << " " << ymin_ << " " << zmin_ << std::endl;
    }
    void computeKeys(const BaseCoordArray& atoms)
    {
        unsigned int stop = atoms.size();
#pragma omp for
        for(unsigned int i=0; i<stop; ++i)
        {
            keys_[i] = pointToKey(atoms.x()[i] - xmin_, atoms.y()[i] - ymin_, atoms.z()[i] - zmin_);
        }
    }
    void insertAtoms()
    {
        map_.clear();
        for(unsigned int i=0; i<keys_.size(); ++i)
        {
            map_.insert(std::pair<Key, unsigned int>(keys_[i],i));
        }
    }
    void keyToVecId()
    {
        // determine number of different keys;
        key_to_id_.clear();
        unsigned int count=0;
        for(unsigned int i=0; i<keys_.size(); ++i)
        {
            std::pair<Map::iterator, bool> it = key_to_id_.insert(std::pair<Key, unsigned int>(keys_[i], count));
            count += it.second;
        }
        if(vmap_.size() < key_to_id_.size()) vmap_.resize(key_to_id_.size());
    }
    void remap()
    {
        NeighbourGrid::KeyArray27 keys;
        unsigned int key_count = key_to_id_.size();
#pragma omp for schedule(dynamic) nowait
        for(unsigned int i=0; i<key_count; ++i)
        {
            Map::const_iterator it = key_to_id_.begin();
            std::advance(it, i);
            this->getKeys(keys, it->first);
            IdVec& nbids = vmap_[it->second];
            nbids.clear();
            for(unsigned int j=0; j<27; ++j)
            {
                std::pair<MultiMap::const_iterator, MultiMap::const_iterator> mit = map_.equal_range(keys[j]);
                for(; mit.first!=mit.second; ++mit.first)
                {
                    nbids.push_back((mit.first)->second);
                }
            }
        }
    }
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	void createNeighbourLists(const BaseCoordArray& atoms) __attribute__((noinline))
	{
        this->init(atoms);
        this->computeKeys(atoms);
#pragma omp single nowait
        this->insertAtoms();
#pragma omp single nowait
        this->keyToVecId();
#pragma omp barrier
        this->remap();
	}
	const IdVec& getNeighbourIds(unsigned int aid) const
	{
	    Map::const_iterator it = key_to_id_.find(keys_[aid]);
	    if(it != key_to_id_.end()) return vmap_[it->second];
	    return empty_;
	}
	NeighbourGrid() {}
};

} // end of namespace

#endif /* NEIGHBOURGRID_H_ */
