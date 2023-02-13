/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef SASA_POINTS2_H_
#define SASA_POINTS2_H_

#include "Typedefs.h"
#include "Atoms.h"
#include "PowerCell.h"
#include "OpenMpUtil.h"

namespace powerborn
{

class SasaPoints2
{
	float probe_rad_, scale_, cutoff_;
	CoordArray sasa_points_[6];
    CoordArray waters_;

	//OpenMP and thread specific data
	openmp::Mutex water_append_mutex_;

	void initSasaPoints(float maxrad, CoordArray& points)
	{
		float theta_0=acos(1.0-cutoff_*cutoff_/(2.0*maxrad*maxrad));
		float pi(M_PI);
		unsigned int lines= pi/theta_0;
		points.clear();
		points.insert(0.0,0.0,1.0);
		for(unsigned int a=1; a < lines; ++a)
		{
			float theta=pi*float(a)/float(lines);
			float phi_0=acos(1.0-cutoff_*cutoff_/(2.0*sin(theta)*sin(theta)*maxrad*maxrad));
			unsigned int pointnum= 2.0 * pi / phi_0;
			for(unsigned int b=0; b < pointnum; ++b)
			{
				float phi=2.0 * pi / float(pointnum)*float(b);
				points.insert(sin(theta)*cos(phi), sin(theta)*sin(phi),cos(theta));
			}
		}
		points.insert(0.0,0.0,-1.0);
		points.fixPadding();
		points.setRad(probe_rad_);
		//printf("created %d atoms on sphere of %f Angstrom\n",(int) points.size(), maxrad);
	}
	void checkSasaPoints(const PowerCell& pc, CoordArray& ap) __attribute__((noinline))
	{
	    const AlignedVector<unsigned int>::Type& real_nb = pc.realNeighbours();
		unsigned int nbsize = real_nb.size();
		float r = pc.radius();
		for(unsigned int b=0; b < nbsize; ++b) // skip first six dummy neighbours
		{
		    unsigned int nbid = real_nb[b];
		    Coord3 nb = pc.position(nbid) * r;
		    float nb_rad = pc.radius(nbid) * r;
			ap.checkSphere(nb[0], nb[1], nb[2], nb_rad);
			if(ap.size()==0) break;
		}
		ap.r().setConstant(probe_rad_);
	}
	void resetSasaPoints(float r, unsigned int sasa_id, CoordArray& ap)
	{
		ap.getScaledCoords(sasa_points_[sasa_id],r);
	}
	void addZeroPoints(const PowerCell& pc, CoordArray& ap) __attribute__((noinline))
	{
	    unsigned int stop = pc.zeroPointSize();
	    for(unsigned int i=0; i<stop; ++i)
	    {
	        Coord3 p = pc.zeroPoint(i).position() * pc.radius();
	        ap.insert(p[0], p[1], p[2], probe_rad_);
	    }
	}
	void init()
	{
		this->initSasaPoints(float(1.5),sasa_points_[0]);
		this->initSasaPoints(float(2.0),sasa_points_[1]);
		this->initSasaPoints(float(3.0),sasa_points_[2]);
		this->initSasaPoints(float(4.0),sasa_points_[3]);
        this->initSasaPoints(float(5.0),sasa_points_[4]);
        this->initSasaPoints(float(10.0),sasa_points_[5]);
		actual_points_.setAll(sasa_points_[5]);
	}
protected:
    openmp::ThreadPrivateObject<CoordArray> actual_points_;
    openmp::ThreadPrivateObject<CoordArray> thread_ac_;
    void watersForCell(const PowerCell& pc, CoordArray& ap) __attribute__((noinline))
    {
        float r = pc.radius();
        unsigned int sasa_id = 0;
        if (r > float(5.0)) sasa_id = 5;
        else if(r > float(4.0)) sasa_id = 4;
        else if(r > float(3.0)) sasa_id = 3;
        else if (r > float(2.0)) sasa_id = 2;
        else if (r > float(1.0)) sasa_id = 1;
        this->resetSasaPoints(r, sasa_id, ap);
        this->checkSasaPoints(pc, ap);
        this->addZeroPoints(pc, ap);
    }

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    float getProbeRadius() const {return probe_rad_;}
	void clearWaters() {waters_.clear();}
	const CoordArray& getWaters() const {return waters_;}
	void append_water(const BaseCoordArray& ac) __attribute__((noinline))
	{
		openmp::ScopedLock lock(water_append_mutex_);
		waters_.append(ac);
	}
	void setParameters(float probe, float cut)
	{
	    probe_rad_ = probe;
	    cutoff_ = cut * probe;
	    this->init();
	}
	SasaPoints2(float pr = 1.4f, float s = 0.8f): probe_rad_(pr), scale_(s), cutoff_(s * pr)
	{
		this->init();
	}
	~SasaPoints2(){}
};

} // end of namespace

#endif /* SASA_POINTS2_H_ */
