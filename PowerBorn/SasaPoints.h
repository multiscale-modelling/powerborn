/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef SASA_POINTS_H_
#define SASA_POINTS_H_

#include <Eigen/StdVector>
#include "Atoms.h"
#include "Typedefs.h"
#include "PDTypedefs.h"
#include "OpenMpUtil.h"

namespace powerborn
{

class SasaPoints: public openmp::OpenMpInterface
{
	typedef Coord3 Coord;
	typedef PowerDiagram::vertex Vertex;
	typedef Vertex* VertexPtr;
	typedef PowerDiagram::cell Cell;
	typedef CoordArray AtomCollectionType;
	typedef AlignedVector<Float4>::Type Neighbours;

	float probe_rad_, scale_, cutoff_;
	CoordArray sasa_points_[6];
	ArrayU is_covered_;
	unsigned int index_;

	//OpenMP and thread specific data
	openmp::ThreadPrivateObject<CoordArray> actual_points_;
	openmp::ThreadPrivateObject<Neighbours> neighbours_;
	openmp::ThreadPrivateObject<CoordArray> thread_ac_;
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
	void checkSasaPoints(CoordArray& ap, const Neighbours& nb)
	{
		unsigned int nbsize = nb.size();
		for(unsigned int b=0; b < nbsize; ++b)
		{
			ap.checkSphere(nb[b][0],nb[b][1],nb[b][2],nb[b][3]);
			if(ap.size()==0) break;
		}
	}
	void resetSasaPoints(float r, unsigned int sasa_id, CoordArray& ap)
	{
		ap.getScaledCoords(sasa_points_[sasa_id],r);
	}
	void surfaceVertices(const PowerDiagram* pd, CoordArray& ac)
	{
		const std::vector<PowerDiagram::zeroPoint>& zp = pd->get_zeroPoints();
		const unsigned int stop = zp.size();
#pragma omp for schedule(dynamic, 16) nowait
		for(unsigned int a=0; a<stop; ++a)
		{
			if (zp[a].pos >= 0.0 and zp[a].pos <= 1.0)
			{
				Coord p = zp[a].getPos() + pd->center;
				ac.insert(p[0], p[1], p[2], probe_rad_);
			}
		}
	}
	void watersForCell(const Cell& c, CoordArray& ap, Neighbours& my_nb)
	{
		my_nb.clear();
		float r = c.r;
		unsigned int sasa_id = 0;
        if (r > float(5.0)) sasa_id = 5;
        else if(r > float(4.0)) sasa_id = 4;
		else if(r > float(3.0)) sasa_id = 3;
		else if (r > float(2.0)) sasa_id = 2;
		else if (r > float(1.0)) sasa_id = 1;
		this->resetSasaPoints(r, sasa_id, ap);
		const std::vector<Cell*>& nb = c.neighbours;
		for(std::vector<Cell*>::const_iterator it=nb.begin(); it!=nb.end(); ++it)
		{
			Coord tmp = (*it)->position - c.position;// pd->center-center;
			float d = r - tmp.norm();
			if( -(*it)->r < d and d < (*it)->r )
			{
				my_nb.push_back(Float4(tmp[0],tmp[1],tmp[2],(*it)->r));
			}
		}
		this->checkSasaPoints(ap, my_nb);
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
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    float getProbeRadius() const
    {
        return probe_rad_;
    }
	void createWaters(AtomCollectionType& waters, const PowerDiagram* pd)
	{
		std::vector<Cell> const & cells = pd->get_points();
		unsigned int stop = cells.size();
#pragma omp single
		{
			is_covered_.resize(stop);
		}
		const std::vector<Vertex>& vertices = pd->get_vertices();
		Cell* allcovered = vertices[0].generators[0];
		bool is_allcovered = true; // check for one atom including all others
		for(int a=0; a<8; ++a)
		{
			is_allcovered &= (allcovered == vertices[a].generators[0]);
		}
		CoordArray& ac = thread_ac_.get();
		ac.clear();
		CoordArray& ap = actual_points_.get();
		Neighbours& nb = neighbours_.get();
#pragma omp for
		for(unsigned int a = 0; a < stop; ++a)
		{
			const Cell& pd_cell = cells[a];
			bool is_covered=true;
			const std::vector<VertexPtr>& verts = pd_cell.myVertices;
			unsigned int vert_size=verts.size();

			for(unsigned int b=0; b < vert_size; ++b)
			{
				is_covered &= verts[b]->powerValue <= float(0.0);
			}
			is_covered = (is_allcovered and &pd_cell == allcovered)? false : is_covered;
			is_covered_[a] = !is_covered;
		}
#pragma omp single
		{
			index_ = 0;
			for(unsigned int i = 0; i < stop; ++i)
			{
				unsigned int tmp = is_covered_[i];
				is_covered_[index_] = i;
				index_ += tmp;
			}
		}
#pragma omp for schedule(dynamic) nowait
		for(unsigned int i = 0; i < index_; ++i)
		{
			this->watersForCell(cells[is_covered_[i]], ap, nb);
			if(ap.size()>0)
			{
				Coord center = cells[is_covered_[i]].position + pd->center;
				ap.translate(center[0],center[1],center[2]);
				ac.append(ap);
			}
		}
		this->surfaceVertices(pd, ac);
		// combine thread ac:
		this->append_water(ac, waters);
	}
	void append_water(const CoordArray& ac, CoordArray& waters)
	{
		openmp::ScopedLock lock(water_append_mutex_);
		waters.append(ac);
	}
	SasaPoints(float pr = 1.4f, float s = 0.8f): probe_rad_(pr), scale_(s), cutoff_(s * pr)
	{
		this->init();
	}
};

} // end of namespace

#endif /* POWERBON_SASA_POINTS_H_ */
