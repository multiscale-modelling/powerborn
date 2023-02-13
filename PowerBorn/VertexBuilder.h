/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef VERTEXBUILDER_H_
#define VERTEXBUILDER_H_

#include "TernaryNet.h"
#include "Atoms.h"
#include "PowerCell.h"

#include <algorithm>
#include <fstream>
#include <sstream>

#undef EIGEN_DEFAULT_IO_FORMAT
#define EI_FORMAT Eigen::IOFormat(3, 0, ", ", ", ", "", "", "", ";")

namespace powerborn
{

struct PlaneValueSorter
{
    typedef std::pair<float, unsigned int> pair;
    bool operator() (const pair a, const pair b) { return (a.first < b.first);}
};

class VertexBuilder
{
private:
    struct VertexTriple
    {
        // non-obsolete, obsolete, new
        Vertex* first, *second, *third;
        VertexTriple(Vertex* f, Vertex* s): first(f), second(s), third(NULL)
        {
        }
        VertexTriple(Vertex* f, Vertex* s, Vertex* t): first(f), second(s), third(t)
        {
        }
        VertexTriple(){}
    };

    PowerVertices net_;

    // only used for constructing, all persisting data is stored in PowerCell
    AlignedVector<Vertex*>::Type gen1_;
    AlignedVector<VertexTriple>::Type connections_;
    AlignedVector<unsigned int>::Type nb_ids_;
    AlignedVector<ZeroPoint>::Type zero_points_;
    AlignedVector<std::pair<float, unsigned int> >::Type plane_values_;

    ArrayB single_circle_excluded_, has_single_circle_;
    CoordArray neighbours_, init_neighbours_, all_neighbours_, vertices_, nv_;

    void initNet() __attribute__((noinline))
    {
        net_.reset();
        AlignedVector<Vertex*>::Type& nodes = net_.nodes();
        if(nodes.size() < 8) nodes.resize(8);
        for(unsigned int i=0; i<8; ++i)
        {
            nodes[i] = &net_.getNewNode();
        }
        // position the nodes
        nodes[0]->connects() = Vertex::Connects(nodes[4], nodes[1], nodes[3]);
        nodes[0]->data().generators = UInt4(0, 1, 2, 4);
        nodes[1]->connects() = Vertex::Connects(nodes[5], nodes[0], nodes[2]);
        nodes[1]->data().generators = UInt4(0, 1, 3, 4);
        nodes[2]->connects() = Vertex::Connects(nodes[6], nodes[1], nodes[3]);
        nodes[2]->data().generators = UInt4(0, 0, 3, 4);
        nodes[3]->connects() = Vertex::Connects(nodes[7], nodes[0], nodes[2]);
        nodes[3]->data().generators = UInt4(0, 0, 2, 4);
        nodes[4]->connects() = Vertex::Connects(nodes[0], nodes[5], nodes[7]);
        nodes[4]->data().generators = UInt4(0, 1, 2, 5);
        nodes[5]->connects() = Vertex::Connects(nodes[1], nodes[4], nodes[6]);
        nodes[5]->data().generators = UInt4(0, 1, 3, 5);
        nodes[6]->connects() = Vertex::Connects(nodes[2], nodes[5], nodes[7]);
        nodes[6]->data().generators = UInt4(0, 0, 3, 5);
        nodes[7]->connects() = Vertex::Connects(nodes[3], nodes[4], nodes[6]);
        nodes[7]->data().generators = UInt4(0, 0, 2, 5);
        net_.setSize(8);
    }
    void computeObsoleteVertices(const Coord3& apos, float rad) __attribute__((noinline))
    {
        vertices_.checkVertices(apos[0], apos[1], apos[2], rad);
    }
    unsigned int findObsoleteConnections(unsigned int count, unsigned int new_atom) __attribute__((noinline))
    {
        if(connections_.size() < count * 3) connections_.resize(count * 3);
        int ccount = 0;
        for(unsigned int i = 0; i < count; ++i)
        {
            Vertex* start = net_.node(vertices_.trueIds()[i]);
            for(int j=0; j<3; ++j)
            {
                Vertex* tmp = start->connects(j);
                connections_[ccount].first = tmp;
                connections_[ccount].second = start;
                unsigned int tmp_id = tmp->id();
                ccount++;
                ccount += vertices_.tmp()[tmp_id];
            }
        }
        for(int i=0; i<ccount; ++i)
        {
            Vertex* nv = &net_.getNewNode();
            connections_[i].third = nv;
            // first connection is the one to an old vertex
            nv->connects(0) = connections_[i].first;
            net_.reconnect(connections_[i].first->connects(), connections_[i].second, nv);
            // new vertex gets the new atom and the two common generators of the old connection
            UInt2 tmp = connections_[i].first->data().commonGenerators(connections_[i].second->data());
            UInt4& gen3 = connections_[i].third->data().generators;
            gen3[1] = new_atom;
            gen3[2] = tmp[0];
            gen3[3] = tmp[1];
        }
        return ccount;
    }
    void connectNewVertices(unsigned int count) __attribute__((noinline))
    {
        // two new vertex have a connection if they share a generator
        for(unsigned int i=0; i<count; i++)
        {
            Vertex* nv = connections_[i].third;
            UInt4& generators = nv->data().generators;
            unsigned int g1 = generators[2];
            unsigned int g2 = generators[3];
            Vertex* nv1 = gen1_[g1];
            Vertex* nv2 = gen1_[g2];
            gen1_[g1] = nv1 == NULL ? nv : NULL;
            gen1_[g2] = nv2 == NULL ? nv : NULL;
            if(nv1 != NULL)
            {
                Vertex::Connects& c1 = nv->connects();
                Vertex::Connects& c2 = nv1->connects();
                c1[2] = c1[1];
                c1[1] = nv1;
                c2[2] = c2[1];
                c2[1] = nv;
            }
            if(nv2 != NULL)
            {
                Vertex::Connects& c1 = nv->connects();
                Vertex::Connects& c2 = nv2->connects();
                c1[2] = c1[1];
                c1[1] = nv2;
                c2[2] = c2[1];
                c2[1] = nv;
            }
        }
    }
    void removeAndInsert(unsigned int add_count) __attribute__((noinline))
    {
        unsigned int stop = vertices_.falseIds().size();
        unsigned int new_size = add_count + stop;
        AlignedVector<Vertex*>::Type& nodes = net_.nodes();
        if(nodes.size() < new_size) nodes.resize(new_size);
        for(unsigned int i=0; i<stop; ++i)
        {
            nodes[i] = nodes[vertices_.falseIds()[i]];
            nodes[i]->setId(i);
        }
        for(unsigned int i=0; i<add_count; ++i)
        {
            nodes[stop+i] = connections_[i].third;
            nodes[stop+i]->setId(stop+i);
        }
        net_.setSize(new_size);
    }
    void computeNewVertexPositions(unsigned int ccount, const Coord3& apos, float rad) __attribute__((noinline))
    {
        // get new Vertex and compute position and power of new vertex
        nv_.setSize(ccount);
        Coord3 p0 = this->getPointOnPlane(apos, rad);
        float planeval = apos.dot(p0);
        for(unsigned int i=0; i<ccount; ++i)
        {
            Coord3 l0 = vertices_.getPos(connections_[i].first->id());
            Coord3 l1 = vertices_.getPos(connections_[i].second->id());
            Coord3 new_pos = this->getLinePlaneIntersection(planeval, apos,
                    l0, l1);
            nv_.x()[i] = new_pos[0];
            nv_.y()[i] = new_pos[1];
            nv_.z()[i] = new_pos[2];
            nv_.r()[i] = new_pos.squaredNorm() - 1.0f;
        }
        vertices_.compressFalseOrdered();
        vertices_.append(nv_);
    }
    Coord3 getPointOnPlane(const Coord3& apos, float rad)
    {
        // TODO numerical stability, clamp x
        float d = apos.squaredNorm();
        d = d < std::numeric_limits<float>::epsilon() ? std::numeric_limits<float>::epsilon() : d;
        float x = (1.0f - rad * rad) * (0.5f / d) + 0.5f;
        //std::cout << "point on plane " << x << " d " << d << std::endl;
        return apos * x;
    }
    Coord3 getLinePlaneIntersection(float planeval, const Coord3& n,
            const Coord3& l0, const Coord3& l1)
    {
        float dot1 = l0.dot(n);
        float dot2 = l1.dot(n);
        float d = (planeval - dot1) / (dot2 - dot1);
        d = ! (d > 0.0f) ? 0.0f : d; // also includes nan test
        d = d > 1.0f ? 1.0f : d;
        return l0 * (1.0f - d) + l1 * d;
    }
    void printGenerators(const Vertex& v) const
    {
        std::cout << v.data().generators.format(EI_FORMAT);
        for(int i=0; i<3; ++i)
        {
            const UInt4& nb_gen = v.connects(i)->data().generators;
            std::cout << " || " << nb_gen.format(EI_FORMAT);
        }
        std::cout << "  node ids: " << v.id();
        for(int i=0; i<3; ++i)
        {
            std::cout << " " << v.connects(i)->id();
        }
        std::cout << std::endl;
    }
    bool checkGenerators() const
    {
        //std::cout << "Generator check for " << net_.size() << " vertices" << std::endl;
        bool all_ok = true;

        for(unsigned int i=0; i<net_.size(); ++i)
        {
            // each Vertex has to have 2 common generators with each of its connected Vertices
            const Vertex* v = net_.node(i);
            const UInt4& my_gen = v->data().generators;
            for(int i=0; i<3; ++i)
            {
                const UInt4& nb_gen = v->connects(i)->data().generators;
                int count = 0;
                for(int j=1; j<4; ++j)
                {
                    for(int k=1; k<4; ++k)
                    {
                        count += int(my_gen[j] == nb_gen[k]);
                    }
                }
                all_ok &= count == 2;
                if(!all_ok)
                {
                    this->printGenerators(*v);
                }
            }
            //std::cout << std::endl;
        }
        if(!all_ok)
        {
            std::cout << "Generator check failed" << std::endl;
        }
        return all_ok;
    }
    void initNeighbours(const NeighbourGrid& nb, const BaseCoordArray& atoms, PowerCell& pc) __attribute__((noinline))
    {
        all_neighbours_.clear();
        all_neighbours_.append(init_neighbours_);
        unsigned int aid = pc.id();
        const AlignedVector<unsigned int>::Type& nb_ids = nb.getNeighbourIds(aid);
        neighbours_.getIds(atoms, nb_ids, aid);
        Coord3 my_pos = atoms.getPos(aid);
        pc.duplicateAtoms() = neighbours_.getRealNeighbours(my_pos[0], my_pos[1], my_pos[2], atoms.r()[aid]);
        float inv_r = 1.0f / atoms.r()[aid];

        all_neighbours_.appendTransformed(neighbours_, my_pos, inv_r, inv_r);
        if(gen1_.size() < all_neighbours_.size()) gen1_.resize(all_neighbours_.size());
        //std::cout << "New Real neighbours " << nb_ids_.size() << " are " << all_neighbours_.size() << std::endl;
    }
    void init(const BaseCoordArray& atoms)
    {
        nb_ids_.resize(atoms.size());
        zero_points_.clear();
        this->initNet();
        vertices_.setSize(8);
        vertices_.x() << -1.01f, -1.01f,  1.01f,  1.01f, -1.01f, -1.01f,  1.01f,  1.01f;
        vertices_.y() <<  1.01f, -1.01f, -1.01f,  1.01f,  1.01f, -1.01f, -1.01f,  1.01f;
        vertices_.z() <<  1.01f,  1.01f,  1.01f,  1.01f, -1.01f, -1.01f, -1.01f, -1.01f;
        vertices_.r() = vertices_.x().square() + vertices_.y().square() + vertices_.z().square() - 1.0f;
        for(unsigned int i=0; i<net_.size(); ++i)
        {
            net_.node(i)->setId(i);
        }
    }
    bool isCoveredCheck()
    {
        bool tmp = (vertices_.r() < 0.0f).all();
        //if(tmp) std::cout << "all powers smaller 0.0f;" << std::endl;
        return tmp;
    }
    bool addVerticesForAtom(const Coord3& apos, float rad, int new_atom_id)
    {
        // find all vertices that need to be replaced
        this->computeObsoleteVertices(apos, rad);
        unsigned int count = vertices_.trueIds().size();
        if(!count) return 0; // no vertices to replace

        // handle case when all vertices will be replaced
        if(count == net_.size())
        {
            //std::cout << "no vertices left" << std::endl;
            net_.setZero();
            return 1;
        }

        // find all connections that lead to obsolete vertices
        // assign generators to new vertices
        unsigned int ccount = this->findObsoleteConnections(count, new_atom_id);

        // connect the new vertices
        this->connectNewVertices(ccount);

        // compute positions for the new vertices
        this->computeNewVertexPositions(ccount, apos, rad);

        // remove old from net and add new ones
        this->removeAndInsert(ccount);

        return this->isCoveredCheck();
    }
    float stableZeroPoints()
    {
        float rad_square_change = 0.0f;
        const float min_power = 1000.0f * std::numeric_limits<float>::epsilon();
        while(((vertices_.r() + rad_square_change).abs() <= min_power).any() )
        {
            rad_square_change += min_power;
        }
        // change the powers and radii;
        if(rad_square_change != 0.0f)
        {
            //std::cout << "Changed squared radii by " << rad_square_change << std::endl;
            vertices_.r() += rad_square_change;
            all_neighbours_.r() = (all_neighbours_.r().square() - rad_square_change).sqrt();
            return sqrt(1.0f - rad_square_change);
        }
        return 1.0f;
    }
    void computeZeroPoints(PowerCell& pc) __attribute__((noinline))
    {
        zero_points_.clear();
        float new_rad = this->stableZeroPoints();
        pc.setNormedRadius(new_rad);
        float new_rad_square = new_rad * new_rad;
        float inv_new_rad = 1.0f / new_rad;
        for(unsigned int i=0; i<net_.size(); ++i)
        {
            Vertex* node1 = net_.node(i);
            unsigned int node1_id = node1->id();
            for(unsigned int j=0; j<3; ++j)
            {
                Vertex* node2 = node1->connects(j);
                unsigned int node2_id = node2->id();
                if(node1_id > node2_id) continue;
                bool is_valid1 = (vertices_.r()[node1_id] > 0.0f);
                bool is_valid2 = (vertices_.r()[node2_id] > 0.0f);
                if (!(is_valid1 | is_valid2)) continue;

                Coord3 x1 = vertices_.getPos(node1_id);
                Coord3 d = vertices_.getPos(node2_id) -x1;
                float a = d.squaredNorm();
                a = a < std::numeric_limits<float>::epsilon() ? std::numeric_limits<float>::epsilon() : a;
                float b = 2.0f * x1.dot(d);
                float c = x1.squaredNorm() - new_rad_square;
                float disc = b * b - 4.0f * a * c;
                if(disc <= 0.0f) continue;

                disc = sqrt(disc);
                float tmp = 0.5f / a;
                float s1 = (-b + disc) * tmp;
                float s2 = (-b - disc) * tmp;
                float numeric_cutoff1 = 0.0f;
                float numeric_cutoff2 = 1.0f - numeric_cutoff1;
                bool range1 = s1 >= numeric_cutoff1 && s1 <= numeric_cutoff2;
                bool range2 = s2 >= numeric_cutoff1 && s2 <= numeric_cutoff2;
                if(is_valid1 && is_valid2)
                {
                    // need to add two zeropoints or none
                    if(range1 && range2)
                    {
                        zero_points_.push_back(ZeroPoint((x1 + s1 * d) * inv_new_rad, node1, node2));
                        zero_points_.push_back(ZeroPoint((x1 + s2 * d) * inv_new_rad, node1, node2));
                    }
                    else if(range1 || range2)
                    {
                        std::cout << "Warning: Zeropoint one invalid range." << std::endl;
                        std::cout << "s1 " << s1 << " s2 " << s2-1 << " power1 " << vertices_.r()[node1_id] <<
                                " power2 " << vertices_.r()[node2_id] <<  " dist " << a << std::endl;
                        zero_points_.push_back(ZeroPoint((x1 + s1 * d) * inv_new_rad, node1, node2));
                        zero_points_.push_back(ZeroPoint((x1 + s2 * d) * inv_new_rad, node1, node2));
                    }
                }
                else
                {
                    // one zeropoint needs to be added, since one power is larger 0, the other smaller
                    if(range1)
                    {
                        zero_points_.push_back(ZeroPoint((x1 + s1 * d) * inv_new_rad, node1, node2));
                    }
                    else if(range2)
                    {
                        zero_points_.push_back(ZeroPoint((x1 + s2 * d) * inv_new_rad, node1, node2));
                    }
                    else
                    {
                        std::cout << "Warning: Zeropoint invalid ranges. Guessing zero point." << std::endl;
                        std::cout << "s1 " << s1 << " s2 " << s2-1 << " power1 " << vertices_.r()[node1_id] <<
                                " power2 " << vertices_.r()[node2_id] << std::endl;
                        // assume both powers are larger than 0, from expierience
                        if (fabs(vertices_.r()[node1_id]) < fabs(vertices_.r()[node2_id]))
                        {
                            zero_points_.push_back(ZeroPoint((x1 + s1 * d) * inv_new_rad, node1, node2));
                        }
                        else
                        {
                            zero_points_.push_back(ZeroPoint((x1 + s2 * d) * inv_new_rad, node1, node2));
                        }
                    }
                }
            }
        }
    }
    void excludeSingleCircles()
    {
        if(single_circle_excluded_.size() < all_neighbours_.size())
        {
            single_circle_excluded_.setConstant(all_neighbours_.size(), 0);
            has_single_circle_.setConstant(all_neighbours_.size(), 0);
        }
        unsigned int stop = net_.size();
        for (unsigned int j = 0; j < stop; ++j)
        {
            const Vertex* node1 = net_.node(j);
            const UInt4& gen1 = node1->data().generators;
            for(unsigned int i=1; i<4; ++i)
            {
                has_single_circle_[gen1[i]] = 1;
            }
            if (vertices_.r()[j] > 0.0f) continue;
            for (int k = 0; k < 3; ++k)
            {
                const Vertex* node2 = node1->connects(k);
                if (node2 > node1) continue;
                // both powers have to be smaller than 0 to exclude single circle
                unsigned int node2_id = node2->id();
                if (vertices_.r()[node2_id] > 0.0f) continue;

                UInt2 tmp = node1->data().commonGenerators(node2->data());
                single_circle_excluded_[tmp[0]] = 1;
                single_circle_excluded_[tmp[1]] = 1;
            }
        }
        unsigned int zp_size = zero_points_.size();
        for(unsigned int i=0; i<zp_size; ++i)
        {
            UInt2 gen = zero_points_[i].generators();
            single_circle_excluded_[gen[0]] = 1;
            single_circle_excluded_[gen[1]] = 1;
        }
        stop = all_neighbours_.size();
        for(unsigned int i=0; i<stop; ++i)
        {
            single_circle_excluded_[i] |= !has_single_circle_[i];
        }
    }
    void setPowerCellData(PowerCell& pc) __attribute__((noinline))
    {
        pc.setNeighbours(all_neighbours_);
        pc.setPowers(vertices_.r());
        pc.setZeroPoints(zero_points_);
        pc.setSingleCircleExcluded(single_circle_excluded_, all_neighbours_.size());
        AlignedVector<unsigned int>::Type& real_nb = pc.realNeighbours();
        real_nb.resize(all_neighbours_.size()-6);
        for(unsigned int i=6; i<all_neighbours_.size(); ++i) // skip the first six dummy neighbours
        {
            real_nb[i-6] = plane_values_[i-6].second;
        }
    }
    std::string getCgo()
    {
        std::stringstream out;
        out << "10.0 2.0" << std::endl;
        out << "2.0 1.0" << std::endl;
        out << "6.0 1.0 0.0 0.0" << std::endl;
        for(unsigned int i=0; i<net_.size(); ++i)
        {
            Vertex* node = net_.node(i);
            for(int i=0; i<3; ++i)
            {
                if(node < node->connects(i))
                {
                    Coord3 p1 = vertices_.getPos(node->id());
                    Coord3 p2 = vertices_.getPos(node->connects(i)->id());
                    out << " 4.0 " << p1[0]  << " " << p1[1] << " "
                            << p1[2] << " 4.0 " << p2[0] << " " << p2[1]
                            << " " << p2[2] << std::endl;
                }
            }
        }
        out << "3.0" << std::endl;
        return out.str();
    }
    void computeAllPlaneValues() __attribute__((noinline))
    {
        unsigned int stop = all_neighbours_.size();
        if(plane_values_.size() < all_neighbours_.size()) plane_values_.resize(all_neighbours_.size());
        for(unsigned int i=6; i<stop; ++i)
        {
            Coord3 apos = all_neighbours_.getPos(i);
            Coord3 p0 = this->getPointOnPlane(apos, all_neighbours_.r()[i]);
            plane_values_[i-6].first = p0.squaredNorm();
            plane_values_[i-6].second = i;
        }
        std::sort(&plane_values_[0], &plane_values_[stop - 6], PlaneValueSorter());
    }
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    VertexBuilder()
    {
        init_neighbours_.setSize(6);
        init_neighbours_.x() << 1.01f, -1.01f, 0.0f, 0.0f, 0.0f, 0.0f;
        init_neighbours_.y() << 0.0f, 0.0f, 1.01f, -1.01f, 0.0f, 0.0f;
        init_neighbours_.z() << 0.0f, 0.0f, 0.0f, 0.0f, 1.01f, -1.01f;
        init_neighbours_.r().setZero();
    }
    bool build(const NeighbourGrid& nb, const BaseCoordArray& atoms, PowerCell& pc) __attribute__((noinline))
    {
        this->init(atoms);
        try
        {
            this->initNeighbours(nb, atoms, pc);
        }
        catch(const CoveredAtomException& e)
        {
            std::cout << "atom recieved covered atom exception " << pc.id() << std::endl;
            return 1;
        }
        this->computeAllPlaneValues();
        for(unsigned int i=6; i<all_neighbours_.size(); ++i) // skip the first six dummy neighbours
        {
            unsigned int nbid = plane_values_[i-6].second;
            Coord3 nb_pos = all_neighbours_.getPos(nbid);
            bool is_covered = this->addVerticesForAtom(nb_pos, all_neighbours_.r()[nbid], nbid);
            if(is_covered)
            {
                //this->check();
                //std::cout << pc.id() << " is covered" << std::endl;
                return 1;
            }
        }
        this->computeZeroPoints(pc);
        this->excludeSingleCircles();
        this->setPowerCellData(pc);
        return 0;
    }
    void visualize(const std::string& fname="vertices.cgo")
    {
        std::ofstream outfile;
        outfile.open(fname.c_str());
        std::cout << "writing CGO file " << fname << std::endl;
        if (outfile.is_open())
        {
            outfile << this->getCgo() << std::endl;
            outfile.close();
        }
    }
    void writeNeighbours(const std::string& filename)
    {
        all_neighbours_.writePqr(filename, 0.0f);
    }
    bool check()
    {
        bool b1 = net_.check();
        bool b2 = this->checkGenerators();
        return b1 && b2;
    }
};

} // end of namespace

#endif /* VERTEXBUILDER_H_ */
