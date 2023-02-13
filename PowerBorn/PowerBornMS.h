/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef POWERBORNMS_H_
#define POWERBORNMS_H_

#include <time.h>

#include "Atoms.h"
#include "Surface.h"
#include "SasaPoints.h"
#include "Integrator.h"
#include "OpenMpUtil.h"

namespace powerborn
{

class PowerBornMSException: public std::exception
{
};

class PowerBornMS: public openmp::OpenMpInterface
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	PowerBornMS(float resolution, float probe_rad, float spherescale, float ifac, float fit1, float fit2, float alpha, float beta):
		sasapoints_(probe_rad, spherescale), surface_(resolution), integrator_(ifac, fit1, fit2, alpha, beta)
	{
	}
    void setProbeRadius(float r, float spherescale = 0.8f)
    {
        sasapoints_ = SasaPoints(r, spherescale);
    }
	const Array& update(const CoordArray& atoms, const PowerDiagram* pd)
	{
		waters_.clear();
        surface_.init(atoms);
        integrator_.init(atoms);
        unsigned int has_exception = 0;
#pragma omp parallel
		{
            try
            {
                sasapoints_.createWaters(waters_, pd);
    #pragma omp barrier
                surface_.build(atoms, waters_);
    #pragma omp barrier
                surface_.finalize();
    #pragma omp barrier
                integrator_.update(atoms, surface_.getRoot(), surface_.getSizeCenter(), sasapoints_.getProbeRadius());

            } catch (...)
            {
                std::cerr << "Unhandled exception in PowerBornMS!" << std::endl;
#pragma omp atomic
                has_exception++;
            }
		}
		if(has_exception)
		{
		    throw PowerBornMSException();
		}
		return integrator_.getBornradii();
	}
	const Array& getBornRadii()
	{
		return integrator_.getBornradii();
	}
	void writeOctree(const std::string& outfile)
	{
		surface_.writeCgo(outfile);
	}
	void writeWaters(const std::string& outfile)
	{
		waters_.writePqr(outfile, 0.0f);
	}
	void timings()
	{
	}
	void resetTimings()
	{
	}
	inline PowerDiagram getPowerDiagram(const CoordArray& atoms)
	{
		std::vector<unsigned int> bond_to;
		std::vector<Coord3> coords;
		std::vector<float> weights;
		bond_to.resize(atoms.size());
		coords.resize(atoms.size());
		weights.resize(atoms.size());
		bond_to[0] = 0;
		coords[0] =  Coord3(atoms.x()[0], atoms.y()[0], atoms.z()[0]);
		weights[0] = atoms.r()[0];
		for (unsigned int i = 1; i < atoms.size(); ++i)
		{
			bond_to[i] = i-1;
			coords[i][0] = atoms.x()[i];
			coords[i][1] = atoms.y()[i];
			coords[i][2] = atoms.z()[i];
			weights[i] = atoms.r()[i];
		}
		return PowerDiagram::create(coords.size(),coords.begin(),weights.begin(),bond_to.begin())
	                        		.with_radiiGiven(1).with_calculate(1).with_cells(1).with_Warnings(0);
	}
	void updatePqr(const std::string& pqrname)
	{
		atoms_.parsePqr(pqrname, &charges_);
		PowerDiagram pd = this->getPowerDiagram(atoms_);
		this->update(atoms_, &pd);
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
	static PowerBornMS getAcc()
	{
		return PowerBornMS(0.25f, 1.4f, 0.8f, 10.0f, 1.0667862f, 0.03516313f, 0.0f, 1.0f);
	}
	static PowerBornMS getFast()
	{
		return PowerBornMS(0.4f, 1.4f, 0.8f, 10.0f, 1.04697536f, 0.03994111f, 0.0f, 1.0f);
	}
	static PowerBornMS getAccNeutral()
	{
		return PowerBornMS(0.25f, 1.4f, 0.8f, 10.0f, 1.0f, 0.0f, 0.0f, 1.0f);
	}
	static PowerBornMS getFastNeutral()
	{
		return PowerBornMS(0.4f, 1.4f, 0.8f, 10.0f, 1.0f, 0.0f, 0.0f, 1.0f);
	}
    static PowerBornMS getAccMembrane()
    {
        return PowerBornMS(0.25f, 1.4f, 0.8f, 10.0f, 1.13015326386f, 0.0f, 0.0f, 1.0f);
    }
    static PowerBornMS getFastMembrane()
    {
        return PowerBornMS(0.4f, 1.4f, 0.8f, 10.0f, 1.14499658222f, 0.0f, 0.0f, 1.0f);
    }
#ifdef WITH_CORRECTION
	static PowerBornMS getAccCorrection()
	{
		return PowerBornMS(0.25f, 1.4f, 0.8f, 10.0f, 1.04838728e+00f, 5.40412131e-05f, -1.25872712f, 1.00361293f);
	}
	static PowerBornMS getFastCorrection()
	{
		return PowerBornMS(0.4f, 1.4f, 0.8f, 10.0f, 1.05279530e+00f, 2.79398869e-04f, -1.22712151f, 1.00507344f);
	}
#endif
private:
	CoordArray waters_, atoms_;
	Array charges_;
	SasaPoints sasapoints_;
	Surface surface_;
	Integrator integrator_;
};

}

#endif /* POWERBORNMS_H_ */
