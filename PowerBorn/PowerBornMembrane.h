/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef POWERBORNMEMBRANE_H_
#define POWERBORNMEMBRANE_H_

#include "PowerDiagram.h"
#include "MembraneSurface.h"
#include "Integrator.h"
#include "SlabIntegrator.h"

namespace powerborn
{

class PowerBornSlab
{
	Diagram diagram_;
	MembraneSurface above_, below_;
	Integrator integrator_above_, integrator_below_;
    SlabIntegrator slab_integrator_;
	Array born_radii_;

	void init(const BaseCoordArray& atoms)
	{
		born_radii_.setZero(atoms.size());
		diagram_.resetException();
		above_.init(atoms);
		below_.init(atoms);
		integrator_above_.init(atoms);
		integrator_below_.init(atoms);
        slab_integrator_.init(atoms);
	}
	void buildPowerDiagram(const BaseCoordArray& atoms)
	{
		diagram_.build(atoms);
#pragma omp barrier
	}
	void buildOctrees(const BaseCoordArray& atoms)
	{
		if (above_.isNeeded()) above_.build(atoms, diagram_.getWaters());
		if (below_.isNeeded()) below_.build(atoms, diagram_.getWaters());
#pragma omp barrier
        if (above_.isNeeded()) above_.finalize();
        if (below_.isNeeded()) below_.finalize();
#pragma omp barrier
	}
	void integrateOctrees(const BaseCoordArray& atoms)
	{
		if (above_.isNeeded()) integrator_above_.addOctree(atoms, above_.getRoot(), above_.getSizeCenter());
		if (below_.isNeeded()) integrator_below_.addOctree(atoms, below_.getRoot(), below_.getSizeCenter());
#pragma omp barrier
	}
	void integrateSlabs(const BaseCoordArray& atoms)
	{
        slab_integrator_.integrate(atoms, above_, below_);
#pragma omp barrier
	}
    void computeBornRadii(const BaseCoordArray& atoms)
    {
        unsigned int stop = born_radii_.size();
        const Array& br1 = integrator_above_.getBornradii();
        const Array& br2 = integrator_below_.getBornradii();
        const Array& br3 = slab_integrator_.getBornradii();
        float pr = diagram_.getProbeRadius();
#pragma omp for nowait
        for(unsigned int i=0; i<stop; ++i)
        {
            float br = integrator_above_.integralToBornRadius(br1[i] + br2[i] + br3[i]);
            born_radii_[i] = atoms.r()[i] - pr > br ? atoms.r()[i] - pr : br;
            if(born_radii_[i] < atoms.r()[i] - diagram_.getProbeRadius())
            {
            	std::cout << "Warning: too smal radius for atom " << i << ": " << born_radii_[i] << " vs " << atoms.r()[i] - diagram_.getProbeRadius() << std::endl;
            }
        }
    }
#ifdef DEBUG_POWERBORN_PARALLEL
    bool testSerial(const BaseCoordArray& atoms)
    {
        Array backup = born_radii_;
        this->init(atoms);
        this->buildPowerDiagram(atoms);
        this->buildOctrees(atoms);
        this->integrateOctrees(atoms);
        this->integrateSlabs(atoms);
        this->computeBornRadii(atoms);
        bool ok = (backup == born_radii_).all();
        if(!ok)
        {
            std::cerr << "PowerBornMembrane parallel invariance test failed!" << std::endl;
        }
        return ok;
    }
#endif
public:
    const Array& update(const BaseCoordArray& atoms)
    {
        this->init(atoms);
        unsigned int has_exception = 0;
#pragma omp parallel
        {
            try
            {
                this->buildPowerDiagram(atoms);
                this->buildOctrees(atoms);
                this->integrateOctrees(atoms);
                this->integrateSlabs(atoms);
                this->computeBornRadii(atoms);

            } catch(...)
            {
                std::cerr << "Unhandeled exception in PowerBornSlab::update!" << std::endl;
#pragma omp atomic
                has_exception++;

            }
        }
        if(has_exception)
        {
            throw PowerSasa2Exception();
        }
#ifdef DEBUG_POWERBORN_PARALLEL
        this->testSerial(atoms);
#endif
        return born_radii_;
    }
    const Array& getSasa() const
    {
        return diagram_.getSasa();
    }
    void writeCgo(const std::string basename)
    {
        if(above_.isNeeded())
        {
            above_.writeCgo(basename + "_above.cgo");
            above_.writeBoundingBoxCgo(basename + "_above_bb.cgo");
        }
        if(below_.isNeeded())
        {
            below_.writeCgo(basename + "_below.cgo");
            below_.writeBoundingBoxCgo(basename + "_below_bb.cgo");
        }
    }
    const Array& getBornRadii()
    {
        return born_radii_;
    }
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    float gbEnergy(const CoordArray& atoms, const Array& q, float e_in = 2.0f, float e_out = 80.0f, float ks = 4.0f)
    {
        float kcal_factor=1389.3545660404325 / 4.1868; // used to convert apbs results to kcal per mol
        float epsilon_hat = (float) (1.0f / e_in  - 1.0f / e_out);
        float prefactor = -kcal_factor * epsilon_hat;
        float efac = -1.0f/ ks;
        double sum = 0.0;
        const Array& born = born_radii_;
        const Array& x = atoms.x();
        const Array& y = atoms.y();
        const Array& z = atoms.z();
        unsigned int stop = atoms.size();
#pragma omp parallel for reduction(+:sum)
        for(unsigned int id1 = 0; id1 < stop; id1++)
        {
            float xp = x[id1], yp = y[id1], zp = z[id1], qp = q[id1], brp = born[id1];
            sum += (double) ((q * qp) * ( ((x-xp).square() + (y-yp).square() + (z-zp).square()) +
                    (born * brp) * ( ((x-xp).square() + (y-yp).square() + (z-zp).square()) *
                            (born * brp).inverse() * efac).exp() ).inverse().sqrt()).sum();
        }
        sum = 0.5 * sum * prefactor;
        return (float) sum;
    }
	PowerBornSlab(float thickness=15.0f, float probe_radius=1.4f):
        // octree resolution and integration factor in agreement with SLIM paper values
	    diagram_(probe_radius),
        above_(0.5f, thickness, true, probe_radius),
        below_(0.5f, -thickness, false, probe_radius),
        integrator_above_(10.0f, 1.13015326386f, 0.0f, 0.0f, 1.0f),
        integrator_below_(10.0f, 1.13015326386f, 0.0f, 0.0f, 1.0f)
    {
        if(thickness < 0.0f)
        {
            std::cerr << "PowerBornSlab thickness must be larger than 0!" << std::endl;
            abort();
        }
        if(probe_radius <= 0.0f)
        {
            std::cerr << "PowerBornSlab probe radius must be larger than 0!" << std::endl;
            abort();
        }
    }
};


class PowerBornMembrane
{
    Diagram diagram_;
    MembraneSurface above_head_, below_head_, above_core_, below_core_;
    Integrator integrator_above_head_, integrator_below_head_, integrator_above_core_, integrator_below_core_;
    SlabIntegrator slab_integrator_head_, slab_integrator_core_;
    Array born_radii_head_, born_radii_core_;
    bool all_above_needed_, all_below_needed_, none_touch_;

    void init(const BaseCoordArray& atoms)
    {
        born_radii_head_.setZero(atoms.size());
        born_radii_core_.setZero(atoms.size());
        diagram_.resetException();
        above_head_.init(atoms);
        below_head_.init(atoms);
        above_core_.init(atoms);
        below_core_.init(atoms);
        integrator_above_head_.init(atoms);
        integrator_below_head_.init(atoms);
        integrator_above_core_.init(atoms);
        integrator_below_core_.init(atoms);
        slab_integrator_core_.init(atoms);
        slab_integrator_head_.init(atoms);

        all_above_needed_ = above_head_.isNeeded() && above_core_.isNeeded();
        all_below_needed_ = below_head_.isNeeded() && below_core_.isNeeded();
        none_touch_ = !above_head_.touchesMembrane() && !below_head_.touchesMembrane()
                         && !below_core_.touchesMembrane() && !above_core_.touchesMembrane();
    }
    void buildPowerDiagram(const BaseCoordArray& atoms)
    {
        diagram_.build(atoms);
    }
    void buildOctrees(const BaseCoordArray& atoms)
    {
        // special case if complete molecule is outside and bounding box does not touch both slabs
        // then only one octree is necessary
        if(all_above_needed_ && none_touch_)
        {
            above_head_.build(atoms, diagram_.getWaters());
        }
        else if(all_below_needed_ && none_touch_)
        {
            below_head_.build(atoms, diagram_.getWaters());
        }
        else
        {
            if (above_head_.isNeeded()) above_head_.build(atoms, diagram_.getWaters());
            if (below_head_.isNeeded()) below_head_.build(atoms, diagram_.getWaters());
            if (above_core_.isNeeded()) above_core_.build(atoms, diagram_.getWaters());
            if (below_core_.isNeeded()) below_core_.build(atoms, diagram_.getWaters());
        }
    }
    void finalizeOctrees()
    {
        if(all_above_needed_ && none_touch_)
        {
            above_head_.finalize();
        }
        else if(all_below_needed_ && none_touch_)
        {
            below_head_.finalize();
        }
        else
        {
            if (above_head_.isNeeded()) above_head_.finalize();
            if (below_head_.isNeeded()) below_head_.finalize();
            if (above_core_.isNeeded()) above_core_.finalize();
            if (below_core_.isNeeded()) below_core_.finalize();
        }
    }
    void integrateOctrees(const BaseCoordArray& atoms)
    {
        if(all_above_needed_ && none_touch_)
        {
            integrator_above_head_.addOctree(atoms, above_head_.getRoot(), above_head_.getSizeCenter());
        }
        else if(all_below_needed_ && none_touch_)
        {
            integrator_below_head_.addOctree(atoms, below_head_.getRoot(), below_head_.getSizeCenter());
        }
        else
        {
            if (above_head_.isNeeded()) integrator_above_head_.addOctree(atoms, above_head_.getRoot(), above_head_.getSizeCenter());
            if (below_head_.isNeeded()) integrator_below_head_.addOctree(atoms, below_head_.getRoot(), below_head_.getSizeCenter());
            if (above_core_.isNeeded()) integrator_above_core_.addOctree(atoms, above_core_.getRoot(), above_core_.getSizeCenter());
            if (below_core_.isNeeded()) integrator_below_core_.addOctree(atoms, below_core_.getRoot(), below_core_.getSizeCenter());
        }
    }
    void integrateSlabs(const BaseCoordArray& atoms)
    {
        slab_integrator_head_.integrate(atoms, above_head_, below_head_);
        slab_integrator_core_.integrate(atoms, above_core_, below_core_);
    }
    void computeBornRadii(const BaseCoordArray& atoms)
    {
        const float *br1, *br2, *br3, *br4, *br5, *br6;
        if(all_above_needed_ && none_touch_)
        {
            br1 = &integrator_above_head_.getBornradii()[0];
            br2 = &integrator_below_head_.getBornradii()[0];
            br3 = &slab_integrator_head_.getBornradii()[0];
            br4 = &integrator_above_head_.getBornradii()[0];
            br5 = &integrator_below_head_.getBornradii()[0];
            br6 = &slab_integrator_core_.getBornradii()[0];
        }
        else if(all_below_needed_ && none_touch_)
        {
            br1 = &integrator_above_head_.getBornradii()[0];
            br2 = &integrator_below_head_.getBornradii()[0];
            br3 = &slab_integrator_head_.getBornradii()[0];
            br4 = &integrator_above_head_.getBornradii()[0];
            br5 = &integrator_below_head_.getBornradii()[0];
            br6 = &slab_integrator_core_.getBornradii()[0];
        }
        else
        {
            br1 = &integrator_above_head_.getBornradii()[0];
            br2 = &integrator_below_head_.getBornradii()[0];
            br3 = &slab_integrator_head_.getBornradii()[0];
            br4 = &integrator_above_core_.getBornradii()[0];
            br5 = &integrator_below_core_.getBornradii()[0];
            br6 = &slab_integrator_core_.getBornradii()[0];
        }
        float pr = diagram_.getProbeRadius();
        unsigned int stop = born_radii_head_.size();
#pragma omp for nowait
        for(unsigned int i=0; i<stop; ++i)
        {
            float br_head = integrator_above_head_.integralToBornRadius(br1[i] + br2[i] + br3[i]);
            float br_core = integrator_above_head_.integralToBornRadius(br4[i] + br5[i] + br6[i]);
            born_radii_head_[i] = atoms.r()[i] - pr > br_head ? atoms.r()[i] - pr : br_head;
            born_radii_core_[i] = atoms.r()[i] - pr > br_core ? atoms.r()[i] - pr : br_core;
            /*if(born_radii_[i] < atoms.r()[i] - diagram_.getProbeRadius())
            {
                std::cout << "Warning: too smal radius for atom " << i << ": " << born_radii_[i] << " vs " << atoms.r()[i] - diagram_.getProbeRadius() << std::endl;
            }*/
        }
    }
#ifdef POWERBORN_DEBUG_PARALLEL
    bool testSerial(const BaseCoordArray& atoms)
    {
        Array backup1 = born_radii_head_;
        Array backup2 = born_radii_core_;
        this->init(atoms);
        this->buildPowerDiagram(atoms);
        this->buildOctrees(atoms);
        this->finalizeOctrees();
        this->integrateOctrees(atoms);
        this->integrateSlabs(atoms);
        this->computeBornRadii(atoms);
        bool test1 = (backup1 == born_radii_head_).all();
        bool test2 = (backup2 == born_radii_core_).all();
        if(!test1 || !test2)
        {
            std::cerr << "PowerBornMembrane parallel invariance test failed: " << test1 << " " << test2 << std::endl;
        }
        return test1 && test2;
    }
#endif // POWERBORN_DEBUG_PARALLEL
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    void update(const BaseCoordArray& atoms)
    {
        this->init(atoms);
        unsigned int has_exception = 0;
#pragma omp parallel
        {
            try
            {
                this->buildPowerDiagram(atoms);
    #pragma omp barrier
                this->buildOctrees(atoms);
                this->integrateSlabs(atoms);
    #pragma omp barrier
                this->finalizeOctrees();
    #pragma omp barrier
                this->integrateOctrees(atoms);
    #pragma omp barrier
                this->computeBornRadii(atoms);
            }
            catch(...)
            {
                std::cerr << "Unhandled exception in PowerBornMembrane!" << std::endl;
#pragma omp atomic
                has_exception++;
            }
        }
        if(has_exception)
        {
            throw PowerSasa2Exception();
        }
#ifdef POWERBORN_DEBUG_PARALLEL
        this->testSerial(atoms);
#endif
    }
    const Array& getSasa() const
    {
        return diagram_.getSasa();
    }
    const Array& getHeadRadii() const
    {
        return born_radii_head_;
    }
    const Array& getCoreRadii() const
    {
        return born_radii_core_;
    }
    float gbEnergy(const CoordArray& atoms, const Array& q, const Array& radii, float e_in = 2.0f, float e_out = 80.0f, float ks = 4.0f)
    {
        float kcal_factor=1389.3545660404325 / 4.1868; // used to convert apbs results to kcal per mol
        float epsilon_hat = (float) (1.0f / e_in  - 1.0f / e_out);
        float prefactor = -kcal_factor * epsilon_hat;
        float efac = -1.0f/ ks;
        double sum = 0.0;
        const Array& born = radii;
        const Array& x = atoms.x();
        const Array& y = atoms.y();
        const Array& z = atoms.z();
        unsigned int stop = atoms.size();
#pragma omp parallel for reduction(+:sum)
        for(unsigned int id1 = 0; id1 < stop; id1++)
        {
            float xp = x[id1], yp = y[id1], zp = z[id1], qp = q[id1], brp = born[id1];
            sum += (double) ((q * qp) * ( ((x-xp).square() + (y-yp).square() + (z-zp).square()) +
                    (born * brp) * ( ((x-xp).square() + (y-yp).square() + (z-zp).square()) *
                            (born * brp).inverse() * efac).exp() ).inverse().sqrt()).sum();
        }
        sum = 0.5 * sum * prefactor;
        return (float) sum;
    }
    float npEnergy(const CoordArray& atoms, float tension=0.03f)
    {
        float c = 0.32f;
        float za = 0.5f;
        float zb = 9.2f;
        float zc = 25.0f;
        float length_scale = 15.0f / above_head_.zBorder();
        float f_ab = (zb-za);
        f_ab *= f_ab * f_ab;
        f_ab = 1.0f / f_ab;
        float f_bc = zc*zc-zb*zb;
        f_bc *= f_bc * f_bc;
        f_bc = 1.0f / f_bc;
        double np_energy = 0.0f;
        for (unsigned int i=0; i<atoms.size(); ++i)
        {
            float z = length_scale * fabs(atoms.z()[i]);
            float rescale = float(1.0);
            if (z < za)
            {
                rescale = float(0.0);
            }
            else if(z < zb)
            {
                float t1 = (z-za);
                float t2 = (3.0*zb-2.0*z-za);
                t1 *= t1;
                rescale = c * t1 * t2 * f_ab;
            }
            else if(z < zc)
            {
                float t1 = (1.0-c);
                float t2 = (z*z-zb*zb);
                float t3 = (3.0*zc*zc-2.0*z*z-zb*zb);
                t2 *= t2;
                rescale = c + t1 * t2 *t3 * f_bc;
            }
            //std::cout << "XXX" << rescale << " " << tension << " " << diagram_.getSasa()[i] << std::endl;
            np_energy += rescale * tension * diagram_.getSasa()[i];
        }
        return np_energy;
    }
    PowerBornMembrane(float thickness_core=11.0f, float thickness_head=15.0f, float probe_radius=1.4f):
        // octree resolution and integration factor in agreement with SLIM paper values
        diagram_(probe_radius),
        above_head_(0.5f, thickness_head, true, probe_radius),
        below_head_(0.5f, -thickness_head, false, probe_radius),
        above_core_(0.5f, thickness_core, true, probe_radius),
        below_core_(0.5f, -thickness_core, false, probe_radius),
        integrator_above_head_(10.0f, 1.13015326386f, 0.0f, 0.0f, 1.0f),
        integrator_below_head_(10.0f, 1.13015326386f, 0.0f, 0.0f, 1.0f),
        integrator_above_core_(10.0f, 1.13015326386f, 0.0f, 0.0f, 1.0f),
        integrator_below_core_(10.0f, 1.13015326386f, 0.0f, 0.0f, 1.0f)
    {
        if(thickness_head < thickness_core)
        {
            std::cerr << "PowerBornMembrane core thickness must be smaller than head thickness!" << std::endl;
            std::cerr << "    current values: " << thickness_core << " " << thickness_head << std::endl;
            abort();
        }
        if(thickness_core < 0.0f || thickness_head < 0.0f)
        {
            std::cerr << "PowerBornMembrane core thickness and head thickness must be larger than 0!" << std::endl;
            abort();
        }
        if(probe_radius <= 0.0f)
        {
            std::cerr << "PowerBornMembrane probe radius must be larger than 0!" << std::endl;
            abort();
        }
    }
};

} // end of namepsace

#endif /* POWERBORNMEMBRANE_H_ */
