/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef POWERBORNMS2_H_
#define POWERBORNMS2_H_


#include <time.h>

#include "Atoms.h"
#include "Surface.h"
#include "PowerDiagram.h"
#include "SasaPoints2.h"
#include "Integrator.h"
#include "OpenMpUtil.h"

namespace powerborn
{

class PowerBornMS2
{
#ifdef DEBUG_POWERBORN_PARALLEL
	Array test1_, test2_;
	void test_serial(const CoordArray& atoms)
	{
	    std::cout << "PowerBornMS2 test parallel serial" << std::endl;
		test1_ = integrator_.getBornradii();
		test2_ = diagram_.getSasa();
        diagram_.resetException();
        surface_.init(atoms);
        integrator_.init(atoms);
		diagram_.build(atoms);
		surface_.build(atoms, diagram_.getWaters());
		surface_.finalize();
		integrator_.update(atoms, surface_.getRoot(), surface_.getSizeCenter(), diagram_.getProbeRadius());
		if(diagram_.hasException())
		{
		    throw PowerSasa2Exception();
		}
		bool passed1 = (test1_ == integrator_.getBornradii()).all();
		bool passed2 = (test2_ == diagram_.getSasa()).all();
		if(! (passed1 && passed2))
		{
			atoms.dump("parallel_failure.dump");
			std::cerr << "Test parallel failure: "  << passed1 << " " << passed2 << std::endl;
			abort();
		}

	}
#endif // DEBUG_POWERBORN_PARALLLEL
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	PowerBornMS2(float resolution, float probe_rad, float spherescale, float ifac, float fit1, float fit2, float alpha, float beta):
		diagram_(probe_rad, spherescale), surface_(resolution), integrator_(ifac, fit1, fit2, alpha, beta),
    t1_("PowerDiagram"), t2_("Surface"), t3_("Integrator")
	{
	}
    void setProbeRadius(float r, float spherescale = 0.8f)
    {
        diagram_.setParameters(r, spherescale);
    }
	const Array& update(const CoordArray& atoms)
	{
	    unsigned int has_exception = 0;
        diagram_.resetException();
        surface_.init(atoms);
        integrator_.init(atoms);
        //t1_.start();
#pragma omp parallel
		{
            try
            {
			diagram_.build(atoms);
#pragma omp barrier
/*#pragma omp single nowait
            {
                t1_.stop();
                t2_.start();
            }*/
			surface_.build(atoms, diagram_.getWaters());
#pragma omp barrier
			surface_.finalize();
#pragma omp barrier
			/*#pragma omp single nowait
            {
                t2_.stop();
                t3_.start();
            }*/
			integrator_.update(atoms, surface_.getRoot(), surface_.getSizeCenter(), diagram_.getProbeRadius());
            } catch (std::exception& e)
            {
                std::cerr << "PowerBornMS2: Unhandeled exception in parallel block! " << e.what() << std::endl;
#pragma omp atomic
                has_exception++;
            }
		}
        //t3_.stop();
		if(has_exception || diagram_.hasException())
		{
		    throw PowerSasa2Exception();
		}
#ifdef DEBUG_POWERBORN_PARALLEL
		this->test_serial(atoms);
#endif
		return integrator_.getBornradii();
	}
	const Array& getBornRadii()
	{
		return integrator_.getBornradii();
	}
	const Array& getSasa()
	{
	    return diagram_.getSasa();
	}
	void writeOctree(const std::string& outfile)
	{
		surface_.writeCgo(outfile);
	}
	void writeWaters(const std::string& outfile)
	{
		diagram_.getWaters().writePqr(outfile, 0.0f);
	}
	void updatePqr(const std::string& pqrname)
	{
		atoms_.parsePqr(pqrname, &charges_);
		this->update(atoms_);
	}
    void timings()
    {
        t1_.print();
        t2_.print();
        t3_.print();
    }
    void resetTimings()
    {
        t1_.reset();
        t2_.reset();
        t3_.reset();
    }
	float gbEnergy(float e_in = 1.0f, float e_out = 80.0f, float ks = 4.0f)
	{
		return this->gbEnergy(atoms_,charges_,e_in, e_out, ks);
	}
	float gbEnergy(const CoordArray& atoms, const Array& q, float e_in = 1.0f, float e_out = 80.0f, float ks = 4.0f)
	{
	    float kcal_factor=1389.3545660404325 / 4.1868; // used to convert apbs results to kcal per mol
	    float epsilon_hat = (float) (1.0f / e_in  - 1.0f / e_out);
	    float prefactor=-kcal_factor * epsilon_hat;
	    float efac = -1.0f/ ks;
	    double sum = 0.0;
	    const Array& born = integrator_.getBornradii();
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
	static PowerBornMS2 getAcc()
	{
		return PowerBornMS2(0.25f, 1.4f, 0.8f, 10.0f, 1.0667862f, 0.03516313f, 0.0f, 1.0f);
	}
	static PowerBornMS2 getFast()
	{
		return PowerBornMS2(0.4f, 1.4f, 0.8f, 10.0f, 1.04697536f, 0.03994111f, 0.0f, 1.0f);
	}
	static PowerBornMS2 getAccNeutral()
	{
		return PowerBornMS2(0.25f, 1.4f, 0.8f, 10.0f, 1.0f, 0.0f, 0.0f, 1.0f);
	}
	static PowerBornMS2 getFastNeutral()
	{
		return PowerBornMS2(0.4f, 1.4f, 0.8f, 10.0f, 1.0f, 0.0f, 0.0f, 1.0f);
	}
    static PowerBornMS2 getAccMembrane()
    {
        return PowerBornMS2(0.25f, 1.4f, 0.8f, 10.0f, 1.13015326386f, 0.0f, 0.0f, 1.0f);
    }
    static PowerBornMS2 getFastMembrane()
    {
        return PowerBornMS2(0.4f, 1.4f, 0.8f, 10.0f, 1.14499658222f, 0.0f, 0.0f, 1.0f);
    }
#ifdef WITH_CORRECTION
	static PowerBornMS2 getAccCorrection()
	{
		return PowerBornMS2(0.25f, 1.4f, 0.8f, 10.0f, 1.04838728e+00f, 5.40412131e-05f, -1.25872712f, 1.00361293f);
	}
	static PowerBornMS2 getFastCorrection()
	{
		return PowerBornMS2(0.4f, 1.4f, 0.8f, 10.0f, 1.05279530e+00f, 2.79398869e-04f, -1.22712151f, 1.00507344f);
	}
#endif
private:
	CoordArray atoms_;
	Array charges_;
	Diagram diagram_;
	Surface surface_;
	Integrator integrator_;
    openmp::Timer t1_, t2_, t3_;
};

}

#endif /* POWERBORNMS2_H_ */
