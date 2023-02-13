/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef POWERDIAGRAM_H_
#define POWERDIAGRAM_H_

#include "Typedefs.h"
#include "NeighbourGrid.h"
#include "VertexBuilder.h"
#include "PowerCell.h"
#include "PowerSasa2.h"
#include "SasaPoints2.h"

namespace powerborn
{

class Diagram: public SasaPoints2
{
    Array sasa_;
    NeighbourGrid nb_;
    unsigned int has_exception_;

    openmp::ThreadPrivateObject<PowerCell> pc_;
    openmp::ThreadPrivateObject<VertexBuilder> vb_;
    openmp::ThreadPrivateObject<PowerSasa2> ps2_;

    unsigned int do_debug_, debug_aid_;
#ifdef DEBUG_POWERBORN_PARALLEL
    AlignedVector<PowerCell>::Type power_cells_;

    void test_serial(const NeighbourGrid& nb, const BaseCoordArray& atoms) __attribute__((noinline))
    {
        std::cout << "TEST SERIAL POWERDIAGRAM" << std::endl;
        unsigned int stop = atoms.size();
        //PowerCell& pc = pc_.get();
        VertexBuilder& vb = vb_.get();
        PowerSasa2& ps2 = ps2_.get();
        for(unsigned int i=0; i<stop; ++i)
        {
            try
            {
                PowerCell pc;
                pc.init(i, atoms);
                pc.isCovered() = vb.build(nb, atoms, pc);
                if(!(pc == power_cells_[i]))
                {
                    std::cout << "Power Cell " << i << " does not match" << std::endl;
                }
                if(pc.isCovered())
                {
                    if(sasa_[i] != 0.0f)
                    {
                        std::cout << "Sasa for " << i << " does not match, should be 0 since it is covered" << std::endl;
                    }
                    continue;
                }
                float sasa = ps2.calc_sasa_single(pc);
                sasa = pc.duplicateAtoms() == 0 ? sasa : sasa / float(pc.duplicateAtoms());
                if(sasa != sasa_[i])
                {
                    std::cout << "Sasa for " << i << " does not match" << std::endl;
                }
            } catch (PowerSasa2Exception& e)
            {
            }
        }
    }
#endif
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Diagram(): SasaPoints2(),
    has_exception_(0), do_debug_(0), debug_aid_(0)
    {
    }
    Diagram(unsigned int aid): SasaPoints2(),
            has_exception_(0), ps2_(PowerSasa2(aid)), do_debug_(1), debug_aid_(aid)
    {
        std::cout << "PowerDiagram debuggin atom " << aid << std::endl;
    }
    Diagram(float proberad): SasaPoints2(proberad),
            has_exception_(0), do_debug_(0), debug_aid_(0)
    {
    }
    Diagram(float proberad, float cut): SasaPoints2(proberad, cut),
            has_exception_(0), do_debug_(0), debug_aid_(0)
    {
    }
    void resetException()
    {
        has_exception_ = 0;
    }
    unsigned int hasException()
    {
        return has_exception_;
    }
    void build(const BaseCoordArray& atoms)
    {
        nb_.createNeighbourLists(atoms);
        this->build(nb_, atoms);
    }
    PowerCell createPowerCell(const BaseCoordArray& atoms, unsigned int aid)
    {
        nb_.createNeighbourLists(atoms);
        PowerCell pc;
        pc.init(aid, atoms);
        pc.isCovered() = vb_.get().build(nb_, atoms, pc);
        return pc;
    }
    void build(const NeighbourGrid& nb, const BaseCoordArray& atoms) __attribute__((noinline))
    {
        // global init stuff
        unsigned int stop = atoms.size();
#pragma omp single nowait
        {
            if(sasa_.size() != stop)
            {
                sasa_.resize(stop);
#ifdef DEBUG_POWERBORN_PARALLEL
                power_cells_.resize(stop);
#endif
            }
            this->clearWaters();
        }

        //thread init stuff
        PowerCell& pc = pc_.get();
        VertexBuilder& vb = vb_.get();
        PowerSasa2& ps2 = ps2_.get();
        CoordArray& ap = actual_points_.get();
        CoordArray& ac = thread_ac_.get();
        ac.clear();
#pragma omp barrier // wait for init stuff and neighbourGrid to finish
#pragma omp for schedule(dynamic) nowait
        for(unsigned int i=0; i<stop; ++i)
        {
            try
            {
                sasa_[i] = 0.0f;
                if(do_debug_ && i==debug_aid_)
                {
                    std::cout << "NEW debugging current power cell" << std::endl;
                }
                pc.init(i, atoms);
                pc.isCovered() = vb.build(nb, atoms, pc);
#ifdef DEBUG_POWERBORN_PARALLEL
                power_cells_[i] = pc;
#endif
                if(do_debug_ && i==debug_aid_)
                {
                    vb.check();
                    vb.writeNeighbours("test_neighbours.pqr");
                    vb.visualize("test_vertices.cgo");
                    atoms.dump("test_atoms.dump");
                    atoms.writePqr("test_atoms.pqr", 0.0);
                    pc.writeZeroPoints("test_zeropoints.pqr");
                    if(pc.isCovered())
                    {
                        std::cout << "debug atom is covered" << std::endl;
                    }
                }
                if(pc.isCovered())
                {
                    continue;
                }
                float sasa = ps2.calc_sasa_single(pc);
                sasa_[i] = pc.duplicateAtoms() == 0 ? sasa : sasa / float(pc.duplicateAtoms() + 1);
                this->watersForCell(pc, ap);
                if(ap.size()>0)
                {
                    Coord center = pc.position();
                    ap.translate(center[0],center[1],center[2]);
                    ac.append(ap);
                }
            } catch (PowerSasa2Exception& e)
            {
                has_exception_ = 1;
                std::cout << e.what() << "for atom " << i << std::endl;
                std::cerr << e.what() << "for atom " << i << std::endl;
                vb.check();
                pc.writeZeroPoints("ps2_exception_zeropoints.pqr");
                vb.writeNeighbours("ps2_exception_neighbours.pqr");
                vb.visualize("ps2_exception_vertices.cgo");
                atoms.dump("atoms_ps2_exception.dump");
                atoms.writePqr("atoms_ps2_exception.pqr");
                //abort();
            }
        }
#ifdef DEBUG_POWERBORN_PARALLEL
#pragma omp single
        {
            this->test_serial(nb, atoms);
        }
#endif
        // combine thread ac:
        this->append_water(ac);
    }
    unsigned int size() const {return sasa_.size();}
    const Array& getSasa() const {return sasa_;}
};


}

#endif /* POWERDIAGRAM_H_ */
