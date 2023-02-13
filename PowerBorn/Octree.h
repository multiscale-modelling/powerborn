/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef OCTREE_H_
#define OCTREE_H_

#include <fstream>
#include <string>

#include "Atoms.h"
#include "ObjectPool.h"
#include "OctreeCell.h"
#include "OpenMpUtil.h"

namespace powerborn {

class Octree: public CubicGrid
{

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Octree(float minsize) : CubicGrid(), end_level_(0), size_center_(Float4::Zero()),
		minsize_(minsize), root_(0)
	{
	}
	Float4 getSizeCenter() const
	{
		return size_center_;
	}
	const OctreeCell& getRoot() const
	{
		return *root_;
	}
	void writeCgo(const std::string name) const{
		std::ofstream outfile;
		outfile.open(name.c_str());
		printf("writing CGO file %s\n",name.c_str());
		if (outfile.is_open()) {
			outfile << "10.0 2.0" << std::endl;
			outfile << "2.0 1.0" << std::endl;
			outfile << "6.0 1.0 0.0 0.0" << std::endl;
			this->writeCgo(outfile, *root_, size_center_);
		}
		outfile << "3.0" << std::endl;
		outfile.close();
	}
    void writeBoundingBoxCgo(const std::string name) const{
        std::ofstream outfile;
        outfile.open(name.c_str());
        printf("writing CGO file %s\n",name.c_str());
        if (outfile.is_open()) {
            outfile << "10.0 2.0" << std::endl;
            outfile << "2.0 1.0" << std::endl;
            outfile << "6.0 1.0 0.0 0.0" << std::endl;
            this->writeSingleCgo(size_center_, outfile);
        }
        outfile << "3.0" << std::endl;
        outfile.close();
    }
	virtual ~Octree(){}

protected:
	int end_level_;
	Float4 size_center_;
	float minsize_;
    OctreeCell* root_;
    openmp::ThreadPrivateObject<Pool<OctreeCell> > cell_pools_;
	void clear()
	{
		for(unsigned int i=0; i<cell_pools_.size(); ++i) {
			cell_pools_.getFromThread(i).reset();
		}
	}
	virtual void initSizeCenter(const BaseCoordArray& atoms){
		size_center_ = atoms.boundingBox();
        float s = minsize_;
        end_level_ = 0;
        while (s < size_center_[0]) {
            end_level_++;
            s *= 2.0f;
        }
        end_level_--;
        size_center_[0] = s;
	}
	void initOctree() {
		this->clear();
		root_ = cell_pools_.get().getObjectPtr();
		this->constructChildCell(*root_, size_center_[0], size_center_[1], size_center_[2], size_center_[3]);
	}
	void constructChildCell(OctreeCell& child, float s, float x, float y, float z) {
		s *= 0.5f;
		child.s.setConstant(s);
		child.x = dirx_ * s + x;
		child.y = diry_ * s + y;
		child.z = dirz_ * s + z;
		child.dtype.setZero();
	}
	void constructChildCell(OctreeCell& parent, OctreeCell& child, unsigned int blockid) {
		parent.setChild(child, blockid);
		parent.dtype[blockid] = 16;
		this->constructChildCell(child, parent.s[blockid], parent.x[blockid], parent.y[blockid], parent.z[blockid]);
	}
	void writeCgo(std::ofstream& outfile, const OctreeCell& c, Float4 sc) const {
		sc[0] *= 0.5f;
		Float8 dx = dirx_ * sc[0] + sc[1];
		Float8 dy = diry_ * sc[0] + sc[2];
		Float8 dz = dirz_ * sc[0] + sc[3];
		for (unsigned int i = 0; i < 8; ++i) {
			Float4 tmp = Float4(sc[0], dx[i], dy[i], dz[i]);
			if (c.dtype[i] == 0) {
				this->writeSingleCgo(tmp, outfile);
			}
			else if (c.dtype[i] > 1) {
				this->writeCgo(outfile, c.getChild(i), tmp);
			}
		}
	}
	void writeSingleCgo(Float4 c, std::ofstream& outfile) const {
		if (outfile.is_open()) {
			typedef Eigen::Vector3f Coord;
			Coord corners[8];
			corners[0] = Coord(1.0, 1.0, 1.0);
			corners[1] = Coord(-1.0, 1.0, 1.0);
			corners[2] = Coord(-1.0, -1.0, 1.0);
			corners[3] = Coord(1.0, -1.0, 1.0);
			corners[4] = Coord(1.0, 1.0, -1.0);
			corners[5] = Coord(-1.0, 1.0, -1.0);
			corners[6] = Coord(-1.0, -1.0, -1.0);
			corners[7] = Coord(1.0, -1.0, -1.0);
			Coord center(c[1], c[2], c[3]);
			float size = c[0];
			for (unsigned int a = 0; a < 8; ++a) {
				corners[a] *= size;
				corners[a] += center;
			}
			for (unsigned int a = 0; a < 4; ++a) {
				outfile << " 4.0 " << corners[a][0] << " " << corners[a][1] << " "
						<< corners[a][2];
				outfile << " 4.0 " << corners[(a + 1) % 4][0] << " "
						<< corners[(a + 1) % 4][1] << " " << corners[(a + 1) % 4][2];
				outfile << " 4.0 " << corners[a + 4][0] << " " << corners[a + 4][1] << " "
						<< corners[a + 4][2];
				outfile << " 4.0 " << corners[(a + 1) % 4 + 4][0] << " "
						<< corners[(a + 1) % 4 + 4][1] << " " << corners[(a + 1) % 4 + 4][2];
				outfile << " 4.0 " << corners[a][0] << " " << corners[a][1] << " "
						<< corners[a][2];
				outfile << " 4.0 " << corners[(a) % 4 + 4][0] << " "
						<< corners[(a) % 4 + 4][1] << " " << corners[(a) % 4 + 4][2];
			}
			outfile << std::endl;
		}
	}
};


} // end of namespace


#endif /* OCTREE_H_ */
