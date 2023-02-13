/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef POWERCELL_H_
#define POWERCELL_H_

#include "TernaryNet.h"
#include "Atoms.h"
#include "NeighbourGrid.h"

#include <set>

namespace powerborn
{

/*
 * these classes are for interfacing with old PowerSasa Code
 */

class ZeroPoint
{
    Coord pos_;
    UInt2 generators_;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    ZeroPoint(const Coord pos, const Vertex* node1, const Vertex* node2):
        pos_(pos), generators_(node1->data().commonGenerators(node2->data()))
    {
    }
    Coord& position() {return pos_;}
    const Coord& position() const {return pos_;}
    UInt2& generators() {return generators_;}
    const UInt2& generators() const {return generators_;}
    bool operator==(const ZeroPoint& other) const
    {
        bool test1 = (pos_ == other.pos_);
        bool test2 = (generators_ == other.generators_).all();
        return test1 & test2;
    }
    std::string str() const
    {
        std::stringstream s;
        s << pos_[0] << " " << pos_[1] << " " << pos_[2] << " " << generators_[0] << " " << generators_[1];
        return s.str();
    }
};

class PowerCell
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    void init(unsigned int aid, const BaseCoordArray& atoms)
    {
        real_neighbours_.clear(); // TODO performance
        powers_.resize(0); // TODO performance
        zero_points_.clear();
        rel_neighbours_.clear();
        single_circle_excluded_.setZero();
        pos_ = Coord(atoms.x()[aid], atoms.y()[aid], atoms.z()[aid]);
        rad_ = atoms.r()[aid];
        id_ = aid;
        is_covered_ = 0;
        duplicate_atoms_ = 0;
        normed_rad_ = 1.0f;
    }
    unsigned int id() const {return id_;}
    void setNeighbours(const BaseCoordArray& other) {rel_neighbours_ = other;}
    const BaseCoordArray& getNeighbours() const {return rel_neighbours_;}
    unsigned int neighbourSize() const {return rel_neighbours_.size();}
    Coord position(unsigned int i) const {return rel_neighbours_.getPos(i);}
    float radius(unsigned int i) const {return rel_neighbours_.r()[i];}

    const Coord& position() const {return pos_;}
    float radius() const {return rad_;}
    float normedRadius() const {return normed_rad_;}
    void setNormedRadius(float r) {normed_rad_ = r;}

    bool& isCovered() {return is_covered_;}
    const bool& isCovered() const {return is_covered_;}

    const float& powerVerticesPower(unsigned int i) const {return powers_[i];}
    void setPowers(const Array& p) {powers_ = p;}

    AlignedVector<unsigned int>::Type& realNeighbours() {return real_neighbours_;}
    const AlignedVector<unsigned int>::Type& realNeighbours() const {return real_neighbours_;}

    void setSingleCircleExcluded(const ArrayB& other, unsigned int s)
    {
        if(single_circle_excluded_.size() < s) single_circle_excluded_.resize(s);
        single_circle_excluded_.head(s) = other.head(s);
    }
    bool singleCircleExcluded(unsigned int i) const {return single_circle_excluded_[i];}

    void setZeroPoints(const AlignedVector<ZeroPoint>::Type& zp) {zero_points_ = zp;}
    const ZeroPoint& zeroPoint(unsigned int i) const {return zero_points_[i];}
    unsigned int zeroPointSize() const {return zero_points_.size();}

    void writeZeroPoints(const std::string name="zeropoints.pqr")
    {
        BaseCoordArray zp;
        zp.insert(0.0, 0.0, 0.0, 1.0);
        for(unsigned int i=0; i<zero_points_.size(); ++i)
        {
            zp.insert(zero_points_[i].position()[0], zero_points_[i].position()[1], zero_points_[i].position()[2], 0.05f);
        }
        zp.writePqr(name, 0.0f);
    }
    bool operator==(const PowerCell& other) const
    {
        if(&other == this) return 1;
        ArrayB test(10);
        test[0] = (rel_neighbours_ == other.rel_neighbours_);
        test[1] = (zero_points_ == other.zero_points_);
        test[2] = (powers_.size() == other.powers_.size());
        if(test[2])
        {
            test[2] = ((powers_ == other.powers_).all());
        }
        test[3] = ((single_circle_excluded_.head(rel_neighbours_.size()) == other.single_circle_excluded_.head(rel_neighbours_.size())).all());
        test[4] = (id_ == other.id_);
        test[5] = (pos_ == other.pos_);
        test[6] = (real_neighbours_ == other.real_neighbours_);
        test[7] = (rad_ == other.rad_);
        test[8] = (is_covered_ == other.is_covered_);
        test[9] = (duplicate_atoms_ == other.duplicate_atoms_);
        Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ", "", "", " ", "");
        if(!test.all()) std::cout << id_ << " TEST " << test.format(fmt) << std::endl;
        return test.all();
    }
    unsigned int& duplicateAtoms() {return duplicate_atoms_;}
    const unsigned int& duplicateAtoms() const {return duplicate_atoms_;}
private:
    BaseCoordArray rel_neighbours_;
    AlignedVector<ZeroPoint>::Type zero_points_;
    Array powers_;
    ArrayB single_circle_excluded_;
    AlignedVector<unsigned int>::Type real_neighbours_;
    Coord pos_;
    float rad_, normed_rad_;
    unsigned int id_;
    bool is_covered_;
    unsigned int duplicate_atoms_;
};

}

#endif /* POWERCELL_H_ */
