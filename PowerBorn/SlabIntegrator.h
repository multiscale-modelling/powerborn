/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef SLABINTEGRATOR_H_
#define SLABINTEGRATOR_H_

#include <math.h>
#include "Typedefs.h"
#include "MembraneSurface.h"

namespace powerborn
{

class SlabIntegrator: protected Integrator
{
	double squareIntegralOutside(double x, double y, double z)
	{
		const double cutoff_x=0.001;
		bool x_large=fabs(x) > cutoff_x;
		if (x_large) //none small
		{
			double invr = 1.0 / sqrt(x*x+y*y+z*z);
			double angley = asin(fabs(y)*invr);
			double anglez = asin(fabs(z)*invr);
			double xinv=1.0/x;
			double yt = y*xinv;
			double zt = z*xinv;
			double ytsq = yt*yt;
			double ztsq = zt*zt;
			double f1=1.0+ytsq;
			double f2=1.0+ztsq;
			double g1a=f1+f2;
			double g1=yt*zt;
			g1*=g1a;
			double g2a=f2+ytsq;
			double g2=f1*f2*g2a;
			double tmp1=sqrt(f1);
			double test3=anglez;
			if (zt<0) test3=-test3;
			double gb_1a=2.0*yt*test3;
			double gb_1b=0.5+f1;
			double gb_1=gb_1a*gb_1b;
			double gb_2=f1*tmp1;
			double tmp2=sqrt(f2);
			double test4=angley;
			if(yt<0) test4=-test4;
			double gc_1a=2.0*zt*test4;
			double gc_1b=0.5+f2;
			double gc_1=gc_1a*gc_1b;
			double gc_2=f2*tmp2;
			double gcb2 =gb_2*gc_2;
			double g_nom2=gb_1*g2*gc_2;
			double g_nom3=gc_1*g2*gb_2;
			g_nom2+=g_nom3;
			double g_nom1=g1*gcb2;
			double g_denom=g2*gcb2;
			g_nom2+=g_nom1;
			double g=-M_PI;
			if (z*y < 0) g=-g;
			g+=(g_nom2)/g_denom;
			double res2=xinv*xinv;
			const double res3=(-1.0/24.0)*xinv;
			return -g*res2*res3;
		}
		else
		{
			double ysq = y * y;
			double zsq = z * z;
			double test3 = -((y*z*(3.0*ysq*ysq + 2.0*ysq*zsq + 3.0*zsq*zsq) +
					3.0* (ysq + zsq) * (zsq*zsq*atan(z/y) +
							ysq*ysq*atan(y/z)))*x)/(96.0 * (ysq*ysq*zsq*zsq * (ysq + zsq)));
			return test3;
		}
	}
	float planeIntegral(const Coord3 p,
			const MembraneSurface& above,
			const MembraneSurface& below)
	{

		double xp[3] = {above.getSizeCenter()[1] - p[0] + above.getSizeCenter()[0], above.getSizeCenter()[2] - p[1] + above.getSizeCenter()[0], above.zBorder() - p[2]};
		double xm[3] = {above.getSizeCenter()[1] - p[0] - above.getSizeCenter()[0], above.getSizeCenter()[2] - p[1] - above.getSizeCenter()[0], below.zBorder() - p[2]};
		double above_contrib = 0.0;
        Coord3 center_above(above.getSizeCenter()[1], above.getSizeCenter()[2], above.getSizeCenter()[3]);
        Coord3 center_below(below.getSizeCenter()[1], below.getSizeCenter()[2], below.getSizeCenter()[3]);
		if(above.touchesMembrane())
		{
			double test1 = -(this->squareIntegralOutside(xp[2],xp[0],xp[1])
					+ this->squareIntegralOutside(xp[2],xm[0],xm[1])
					- this->squareIntegralOutside(xp[2],xm[0],xp[1])
					- this->squareIntegralOutside(xp[2],xp[0],xm[1]));
			double test2 = this->cubeIntegral(p, center_above, above.getSizeCenter()[0], false, true);
			above_contrib += test1;
			above_contrib += test2;
		}
		else if(above.zBorder() != 0.0f)
		{
			double d = above.zBorder() - p[2]; // distance to membrane boundary
			above_contrib += (M_PI / 6.0) / (d * d * d);
			if(above.isNeeded())
			{
				above_contrib += this->cubeIntegral(p, center_above, above.getSizeCenter()[0], true, true);
			}
		}
		else if(above.zBorder() == 0.0f)
		{
            above_contrib += this->cubeIntegral(p, center_above, above.getSizeCenter()[0], true, true);
		}

		double below_contrib = 0.0;
		if(below.touchesMembrane())
		{
			below_contrib += this->squareIntegralOutside(xm[2],xp[0],xp[1])
					+ this->squareIntegralOutside(xm[2],xm[0],xm[1])
					- this->squareIntegralOutside(xm[2],xm[0],xp[1])
					- this->squareIntegralOutside(xm[2],xp[0],xm[1]);
			below_contrib += this->cubeIntegral(p, center_below, below.getSizeCenter()[0], true, false);
		}
		else if(below.zBorder() != 0.0f)
		{
			double d = below.zBorder() - p[2]; // distance to membrane boundary
			below_contrib -= (M_PI / 6.0) / (d * d * d);
			if(below.isNeeded())
			{
				below_contrib += this->cubeIntegral(p, center_below, below.getSizeCenter()[0], true, true);
			}
		}
		else if(below.zBorder() == 0.0f)
		{
            below_contrib += this->cubeIntegral(p, center_below, below.getSizeCenter()[0], true, true);
		}
		return (float) (above_contrib + below_contrib);
	}

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    SlabIntegrator(): Integrator(10.0f, 1.0f, 0.0f, 1.0f, 1.0f) // values are not important, since only analytic solution is used
    {
    }
    void init(const BaseCoordArray& atoms)
    {
        Integrator::init(atoms);
    }
    const Array& getBornradii()
    {
        return Integrator::getBornradii();
    }
	void integrate(const BaseCoordArray& atoms,
			const MembraneSurface& surface_above,
			const MembraneSurface& surface_below)
	{
        unsigned int stop = atoms.size();
#pragma omp for schedule(dynamic, 16) nowait
        for(unsigned int i=0; i<stop; ++i)
        {
            born_radii_[i] += this->planeIntegral(atoms.getPos(i), surface_above, surface_below);
        }
	}
};

} // end of namespace


#endif /* SLABINTEGRATOR_H_ */
