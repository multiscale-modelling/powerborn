/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include "Atoms.h"
#include "OctreeCell.h"
#include "Typedefs.h"

namespace powerborn {

class Integrator
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Integrator(float ifac, float fit1, float fit2, float alpha, float beta):
			integration_factor_(ifac), fit1_(fit1), fit2_(fit2), alpha_(alpha), beta_(beta)
	{
#if defined(__SSE__) && !defined(POWERBORN_NO_SSE)
		const uint32_t one = 1;
		ione = _mm_load1_ps((float*) &one);
		const uint32_t sixteen = 16;
		isixteen = _mm_load1_ps((float*) &sixteen);
#endif
		num_grid_ << 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0;
		//printf("iface: %f  fit1: %f  fit2 %f  alpha %f  beta %f", ifac, fit1, fit2, alpha, beta);
	}
	const Array& update(const BaseCoordArray& atoms, const OctreeCell& root, const Float4 size_center, float probe_radius)
	{
		unsigned int stop = atoms.size();
		const float size = size_center[0] * 0.5f; // size of first 8 cells is half that of the bounding box
		Coord center(size_center[1], size_center[2], size_center[3]);

#pragma omp for schedule(dynamic) nowait
		for(unsigned int i=0; i<stop; ++i)
		{
			float integral = this->integrateOctree(size, i, atoms, root);
			Coord point(atoms.x()[i], atoms.y()[i], atoms.z()[i]);
			float box_integral = this->cubeIntegral(point, center, size_center[0]);
			float br = this->integralToBornRadius(integral + box_integral);
			born_radii_[i] = atoms.r()[i] - probe_radius > br ? atoms.r()[i] - probe_radius : br;
		}
		return born_radii_;
	}
    void addOctree(const BaseCoordArray& atoms, const OctreeCell& root, const Float4 size_center)
    {
        unsigned int stop = atoms.size();
        const float size = size_center[0] * 0.5f; // size of first 8 cells is half that of the bounding box
        Coord center(size_center[1], size_center[2], size_center[3]);
#pragma omp for schedule(dynamic) nowait
        for(unsigned int i=0; i<stop; ++i)
        {
            born_radii_[i] = this->integrateOctree(size, i, atoms, root);
        }
    }
	const Array& getBornradii()
	{
		return born_radii_;
	}
	void init(const BaseCoordArray& atoms)
	{
		born_radii_.setZero(atoms.size());
	}
	float cubeIntegral(const Coord& p, const Coord& center, float size, bool bottom=true, bool top=true)
	{
		Coord tmp = p - center;
		Coord xpf(size, size, size);
		Coord xmf(-size, -size, -size);
		xpf -= tmp;
		xmf -= tmp;

		double angley[24], anglez[24];
		double xp[3]={xpf[0],xpf[1],xpf[2]};
		double xm[3]={xmf[0],xmf[1],xmf[2]};
		this->computeAngles(xm,xp,angley,anglez);
		double contribs[24];
		for(unsigned int a=0; a<24; ++a)
		{
			contribs[a]=0.0;
		}
		contribs[0] = +this->squareIntegral(xm[0], xp[1], xp[2], angley[0], anglez[0]); // back
		contribs[1] = -this->squareIntegral(xm[0], xp[1], xm[2], angley[1], anglez[1]);
		contribs[2] = -this->squareIntegral(xm[0], xm[1], xp[2], angley[2], anglez[2]);
		contribs[3] = +this->squareIntegral(xm[0], xm[1], xm[2], angley[3], anglez[3]);
		contribs[4] = -this->squareIntegral(xp[0], xp[1], xp[2], angley[4], anglez[4]); //front
		contribs[5] = +this->squareIntegral(xp[0], xp[1], xm[2], angley[5], anglez[5]);
		contribs[6] = +this->squareIntegral(xp[0], xm[1], xp[2], angley[6], anglez[6]);
		contribs[7] = -this->squareIntegral(xp[0], xm[1], xm[2], angley[7], anglez[7]);
		contribs[8] = +this->squareIntegral(xm[1], xp[0], xp[2], angley[8], anglez[8]); //left
		contribs[9] = -this->squareIntegral(xm[1], xp[0], xm[2], angley[9], anglez[9]);
		contribs[10] = -this->squareIntegral(xm[1], xm[0], xp[2], angley[10], anglez[10]);
		contribs[11] = +this->squareIntegral(xm[1], xm[0], xm[2], angley[11], anglez[11]);
		contribs[12] = -this->squareIntegral(xp[1], xp[0], xp[2], angley[12], anglez[12]); //right
		contribs[13] = +this->squareIntegral(xp[1], xp[0], xm[2], angley[13], anglez[13]);
		contribs[14] = +this->squareIntegral(xp[1], xm[0], xp[2], angley[14], anglez[14]);
		contribs[15] = -this->squareIntegral(xp[1], xm[0], xm[2], angley[15], anglez[15]);
		if (bottom)
		{
			contribs[16] = +this->squareIntegral(xm[2], xp[1], xp[0], angley[16], anglez[16]); //bottom
			contribs[17] = -this->squareIntegral(xm[2], xp[1], xm[0], angley[17], anglez[17]);
			contribs[18] = -this->squareIntegral(xm[2], xm[1], xp[0], angley[18], anglez[18]);
			contribs[19] = +this->squareIntegral(xm[2], xm[1], xm[0], angley[19], anglez[19]);
		}
		if(top)
		{
			contribs[20] = -this->squareIntegral(xp[2], xp[1], xp[0], angley[20], anglez[20]); //top
			contribs[21] = +this->squareIntegral(xp[2], xp[1], xm[0], angley[21], anglez[21]);
			contribs[22] = +this->squareIntegral(xp[2], xm[1], xp[0], angley[22], anglez[22]);
			contribs[23] = -this->squareIntegral(xp[2], xm[1], xm[0], angley[23], anglez[23]);
		}
		double integral[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		for (unsigned int a=0; a<4; ++a)
		{
			integral[0]+=contribs[a];
			integral[1]+=contribs[a+4];
			integral[2]+=contribs[a+8];
			integral[3]+=contribs[a+12];
			integral[4]+=contribs[a+16];
			integral[5]+=contribs[a+20];
		}
		return (float) ((integral[0]+integral[1])+(integral[2]+integral[3])+(integral[4]+integral[5]));
	}
    float integralToBornRadius(float integral)
    {
        const float pi(3.0f / 4.0f / M_PI);
        float tmp = pow(pi * integral, 0.33333333333333333333f);
        return float(1.0f) / (fit1_ * tmp + fit2_);
    }
protected:
	typedef Eigen::Vector3f Coord;
	Float8 num_grid_;
	Array born_radii_;
	float integration_factor_, fit1_, fit2_, alpha_, beta_;

	float integrateOctree(float size, unsigned int i, const BaseCoordArray& atoms, const OctreeCell& root)
	{
		float integral;
#if defined(__SSE__) && !defined(POWERBORN_NO_SSE)
		__m128 x,y,z;
		x = _mm_set1_ps(atoms.x()[i]);
		y = _mm_set1_ps(atoms.y()[i]);
		z = _mm_set1_ps(atoms.z()[i]);
		__m128 mintegral = this->integrate(x, y, z, size, root);
		float *p = (float*) &mintegral;
		integral = p[0] + p[2] + p[1] + p[3];
#else
		integral = this->integrate(atoms.x()[i], atoms.y()[i], atoms.z()[i], size, root);
#endif
		return integral;
	}
#if defined(__SSE__) && !defined(POWERBORN_NO_SSE)
	__m128 ione, isixteen;
	inline __m128 distanceSquare(const __m128& x, const __m128& y, const __m128& z,
			const float *cx, const float* cy, const float * cz)
	{
		__m128 cx1 = _mm_load_ps(cx);
		__m128 cy1 = _mm_load_ps(cy);
		__m128 cz1 = _mm_load_ps(cz);
		cx1 = _mm_sub_ps(cx1,x);
		cy1 = _mm_sub_ps(cy1,y);
		cz1 = _mm_sub_ps(cz1,z);
		cx1 = _mm_mul_ps(cx1,cx1);
		cy1 = _mm_mul_ps(cy1,cy1);
		cz1 = _mm_mul_ps(cz1,cz1);
		cx1 = _mm_add_ps(cx1,cy1);
		return _mm_add_ps(cx1,cz1);
	}
#ifndef WITH_CORRECTION
	inline __m128 fastIntegrand(const __m128& d, const float *cs, const __m128& /*tmp*/)
	{
		__m128 res1 = _mm_load_ps(cs);
		__m128 invd = _mm_rcp_ps(d);
		res1 = _mm_mul_ps(res1, invd);
		invd = _mm_mul_ps(invd, invd);
		return _mm_mul_ps(res1, invd);
#else
	inline __m128 fastIntegrand(const __m128& d, const float *cs, const __m128& tmp)
	{
		const __m128 alpha = _mm_set1_ps(alpha_);
		const __m128 beta = _mm_set1_ps(beta_);
		__m128 res1 = _mm_mul_ps(beta, _mm_load_ps(cs));
		__m128 d1 = _mm_add_ps(d,_mm_mul_ps(alpha, tmp));
		__m128 invd = _mm_rcp_ps(d1);
		res1 = _mm_mul_ps(res1, invd);
		invd = _mm_mul_ps(invd, invd);
		return _mm_mul_ps(res1, invd);
#endif
	}
	inline __m128 integrate(const __m128& x, const __m128& y, const __m128& z, const float s, const OctreeCell& c)
	{
		__m128 cx1 = this->distanceSquare(x, y, z, &c.x[0], &c.y[0], &c.z[0]);
		__m128 cx2 = this->distanceSquare(x, y, z, &c.x[4], &c.y[4], &c.z[4]);
		__m128 ts = _mm_set1_ps(s * s);
		const __m128 ifac = _mm_set1_ps(integration_factor_);
		__m128 threshold = _mm_mul_ps(ts, ifac);
		__m128 mask1 = _mm_cmpgt_ps(cx1, threshold); // criteria to integrate cell or not
		__m128 mask2 = _mm_cmpgt_ps(cx2, threshold);

		__m128 dmask1 = _mm_load_ps((float*) &c.dtype[0]);
		__m128 dmask2 = _mm_load_ps((float*) &c.dtype[4]);
		dmask1 = _mm_andnot_ps(mask1, dmask1);
		dmask2 = _mm_andnot_ps(mask2, dmask2);
		UInt8 mask;
		_mm_store_ps((float*) &mask[0], dmask1);
		_mm_store_ps((float*) &mask[4], dmask2);

		__m128 res = this->fastIntegrand(cx1, &c.s[0], ts);
		res = _mm_and_ps(mask1, res);
		__m128 res2 = this->fastIntegrand(cx2, &c.s[4], ts);
		res2 = _mm_and_ps(mask2, res2);
		res = _mm_add_ps(res,res2);

		Coord tmp( ((float *) &x)[0], ((float *) &y)[0] ,((float *) &z)[0]);
		__m128 nres = _mm_setzero_ps();
		float *p = (float*) &nres;
		float newsize = 0.5f * s;
		for(unsigned int i=0; i<8; ++i)
		{
			switch(mask[i])
			{
			case 0: {
				break;
			}
			case 1: {
				Coord center(c.x[i], c.y[i], c.z[i]);
				p[i%4] += this->numericalVolumeIntegral(tmp,center,s);
				break;
			}
			case 16: {
				nres = _mm_add_ps(nres, this->integrate(x, y, z, newsize, c.getChild(i)));
				break;
			}
			}
		}
		res = _mm_add_ps(res, nres);
		return res;
	}
	__m128 numericalIntegrand(__m128 xysq, const __m128 zsq)
	{
	    xysq = _mm_add_ps(xysq, zsq);
	    __m128 res = _mm_mul_ps(xysq, xysq);
	    res = _mm_mul_ps(res, xysq);
	    return _mm_rcp_ps(res); //  1/r^6
	    //return _mm_div_ps(_mm_set1_ps(1.0f),res);
	}
	float numericalVolumeIntegral(const Coord& point, const Coord& center, float size)
	{
	    float step = 0.25f * size;
	    float offset = 0.5f * step - size;
	    Coord xm = center - point;

	    xm[0] += offset;
	    xm[1] += offset;
	    xm[2] += offset;

	    __m128 step4 = _mm_set1_ps(step);
	    __m128  x = _mm_set1_ps(xm[0]);
	    __m128 z1 = _mm_mul_ps(_mm_load_ps((float*) &num_grid_[0]),step4);
	    __m128 z2 = _mm_mul_ps(_mm_load_ps((float*) &num_grid_[4]),step4);
	    z1 = _mm_add_ps(z1,_mm_set1_ps(xm[2]));
	    z2 = _mm_add_ps(z2,_mm_set1_ps(xm[2]));
	    __m128 res1 = _mm_setzero_ps();
	    __m128 res2 = _mm_setzero_ps();

	    z1 = _mm_mul_ps(z1,z1);
	    z2 = _mm_mul_ps(z2,z2);
	    for(unsigned int a=0; a<8; a++)
	    {
	        __m128 xsq = _mm_mul_ps(x,x);
#ifdef WITH_CORRECTION
	        // numerical correction to better fit analytical formula (ACC used only)
	        const __m128 f1 = _mm_set1_ps(-0.47379885f);
	        xsq = _mm_add_ps(xsq, _mm_mul_ps(step4, f1));
#endif
	        __m128  y = _mm_set1_ps(xm[1]);
	        for(unsigned int b=0; b<8; b++)
	        {
	            __m128 ysq = _mm_mul_ps(y,y);
	            ysq = _mm_add_ps(xsq,ysq);
	            res1 = _mm_add_ps(res1, this->numericalIntegrand(ysq,z1));
	            res2 = _mm_add_ps(res2, this->numericalIntegrand(ysq,z2));
	            y = _mm_add_ps(y,step4);
	        }
	        x = _mm_add_ps(x,step4);
	    }
	    res1 = _mm_add_ps(res1, res2);
	    float* p = (float*) &res1;
#ifndef WITH_CORRECTION
	    float vol = step * step * step;
#else
	    float vol = step * step * (step * 0.99629983f); // numerical correction to better fit analytic formula (ACC used only)
#endif
	    return vol*((p[0]+p[2]) + (p[1]+p[3]));
	}
#else
	float integrate(const float x, const float y, const float z, const float s, const OctreeCell& c)
	{
		float tmp = s * s;
		float threshold = tmp * integration_factor_;
		const Float8 d = (c.x-x).square() + (c.y-y).square() + (c.z-z).square();
		const UInt8 mask = (d > threshold).cast<unsigned int>();
#ifndef WITH_CORRECTION
		Float8 res = (mask != 0 ).select(c.s / d.cube(), 0.0f);
#else
		Float8 res = (mask != 0 ).select(c.s * beta_ / (d + alpha_ * tmp).cube(),0.0f);
#endif
		for(unsigned int i=0; i<8; ++i)
		{
			if (c.dtype[i] == 1 && mask[i]==0)
			{
				Coord tmp(x, y ,z);
				Coord center(c.x[i], c.y[i], c.z[i]);
				res[i] = this->numericalVolumeIntegral(tmp,center,s);
			}
			else if(c.dtype[i] == 16 && mask[i] == 0)
			{
				res[i] = this->integrate(x, y, z, 0.5f * s, c.getChild(i));
			}
		}
		return res.sum();
	}
	float numericalVolumeIntegral(const Coord& point, const Coord& center, float size)
	{
		float step = 0.25f * size;
	    float offset = 0.5f * step - size;
	    Coord xm = center - point + Coord(offset, offset, offset);

	    Float8 res = Float8::Zero();
	    Float8 z = (num_grid_ * step + xm[2]).square();
	    float x = xm[0];
	    for(unsigned int a=0; a<8; a++)
	    {
#ifndef WITH_CORRECTION
	        float xsq = x * x;
#else
	        float xsq = x * x + (-0.47379885f * step * step); // numerical correction to better fit analytic formula (ACC used only)
#endif
	        float  y = xm[1];
	        for(unsigned int b=0; b<8; b++)
	        {
	            float ysq = y * y;
	            ysq += xsq;
	            res += (z + ysq).cube().inverse();
	            y += step;
	        }
	        x += step;
	    }
#ifndef WITH_CORRECTION
	    float vol = step * step * step;
#else
	    float vol = step * step * (step * 0.99629983f); // numerical correction to better fit analytic formula (ACC used only)
#endif
	    return vol*res.sum();
	}
#endif
	void computeAngles(const double* xm, const double* xp,
			double* angley, double* anglez)
	{
		double xm0=fabs(xm[0]), xm1=fabs(xm[1]), xm2=fabs(xm[2]);
		double xp0=fabs(xp[0]), xp1=fabs(xp[1]), xp2=fabs(xp[2]);

		double invr1=this->invR(xm0,xp1,xp2);
		angley[0] = this->calcAngle(invr1, xp1);
		anglez[0] = this->calcAngle(invr1, xp2);
		angley[14] = this->calcAngle(invr1, xm0);

		double invr2=this->invR(xm0,xp1,xm2);
		angley[1] = this->calcAngle(invr2, xp1);
		anglez[1] = this->calcAngle(invr2, xm2);
		angley[15] = this->calcAngle(invr2, xm0);

		double invr3=this->invR(xm0,xm1,xp2);
		angley[2] = this->calcAngle(invr3, xm1);
		anglez[2] = this->calcAngle(invr3, xp2);
		angley[10] = this->calcAngle(invr3, xm0);

		double invr4=this->invR(xm0,xm1,xm2);
		angley[3] = this->calcAngle(invr4, xm1);
		anglez[3] = this->calcAngle(invr4, xm2);
		angley[11] = this->calcAngle(invr4, xm0);

		double invr5=this->invR(xp0,xp1,xp2);
		angley[4] = this->calcAngle(invr5, xp1);
		anglez[4] = this->calcAngle(invr5, xp2);
		angley[12] = this->calcAngle(invr5, xp0);

		double invr6=this->invR(xp0,xp1,xm2);
		angley[5] = this->calcAngle(invr6, xp1);
		anglez[5] = this->calcAngle(invr6, xm2);
		angley[13] = this->calcAngle(invr6, xp0);

		double invr7=this->invR(xp0,xm1,xp2);
		angley[6] = this->calcAngle(invr7, xm1);
		anglez[6] = this->calcAngle(invr7, xp2);
		angley[8] = this->calcAngle(invr7, xp0);

		double invr8=this->invR(xp0,xm1,xm2);
		angley[7] = this->calcAngle(invr8, xm1);
		anglez[7] = this->calcAngle(invr8, xm2);
		angley[9] = this->calcAngle(invr8, xp0);

		angley[21] = angley[0];// = this->calcAngle(xp[1], r);
		anglez[14] = anglez[0];// = this->calcAngle(xp[2], r);
		angley[17] = angley[1];// = this->calcAngle(xp[1], r);
		anglez[15] = anglez[1];// = this->calcAngle(xm[2], r);

		angley[23] = angley[2];// = this->calcAngle(xm[1], r);
		anglez[10] = anglez[2];// = this->calcAngle(xp[2], r);
		angley[19] = angley[3];// = this->calcAngle(xm[1], r);
		anglez[11] = anglez[3];// = this->calcAngle(xm[2], r);

		angley[20] = angley[4];// = this->calcAngle(xp[1], r);
		anglez[12] = anglez[4];// = this->calcAngle(xp[2], r);
		angley[16] = angley[5];// = this->calcAngle(xp[1], r);
		anglez[13] = anglez[5];// = this->calcAngle(xm[2], r);

		angley[22] = angley[6];// = this->calcAngle(xm[1], r);
		anglez[ 8] = anglez[6];// = this->calcAngle(xp[2], r);
		angley[18] = angley[7];// = this->calcAngle(xm[1], r);
		anglez[ 9] = anglez[7];// = this->calcAngle(xm[2], r);

		anglez[22] =angley[8];// = this->calcAngle(xp[0], r);
		anglez[18] =angley[9];// = this->calcAngle(xp[0], r);
		anglez[23] =angley[10];// = this->calcAngle(xm[0], r);
		anglez[19] =angley[11];// = this->calcAngle(xm[0], r);
		anglez[20] =angley[12];// = this->calcAngle(xp[0], r);
		anglez[16] =angley[13];// = this->calcAngle(xp[0], r);
		anglez[21] =angley[14];// = this->calcAngle(xm[0], r);
		anglez[17] =angley[15];// = this->calcAngle(xm[0], r);
	}
	double squareIntegral(double x, double y, double z, double angley, double anglez)
    {
    	// used only for bounding box
    	double xinv = 1.0/x;
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
    	g_nom2 /= g_denom;
    	double res2=xinv*xinv;
    	const double res3=(-1.0/24.0)*xinv;
    	return g_nom2*res2*res3;
    }
	double invR(double x, double y, double z)
	{
		x*=x;
		y*=y;
		z*=z;
		return 1.0/sqrt(x+y+z);
	}
    double calcAngle(double r, double y)
    {
    	return asin(y*r);
    }
};

} // end of namespace


#endif /* INTEGRATOR_H_ */
