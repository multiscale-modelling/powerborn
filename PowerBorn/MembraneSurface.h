/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef MEMBRANESURFACE_H_
#define MEMBRANESURFACE_H_

#include "Surface.h"

namespace powerborn
{

class MembraneSurface: public Surface
{
	float z_border_, probe_radius_;
	bool above_, is_needed_, touches_membrane_; // wether the octree should be above or below the z_border
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	MembraneSurface(float resolution, float z_border, bool above, float probe_radius=1.4f): Surface(resolution),
		z_border_(z_border), probe_radius_(probe_radius),
		above_(above), is_needed_(false), touches_membrane_(false)
	{
	}
	float zBorder() const
	{
		return z_border_;
	}
	bool isNeeded() const
	{
		return is_needed_;
	}
	bool touchesMembrane() const
	{
		return touches_membrane_;
	}
	bool isAbove() const
	{
		return above_;
	}
	void oldBoundingBox(const BaseCoordArray& atoms, float probe_radius)
	{
	    // old bounding box computation to be compatible with SLIM paper version
	    atoms.boundingBox();
	    float xmax = (atoms.x() + atoms.r() - probe_radius).maxCoeff();
	    float xmin = (atoms.x() - atoms.r() + probe_radius).minCoeff();
	    float ymax = (atoms.y() + atoms.r() - probe_radius).maxCoeff();
	    float ymin = (atoms.y() - atoms.r() + probe_radius).minCoeff();
	    float zmax = (atoms.z() + atoms.r() - probe_radius).maxCoeff();
	    float zmin = (atoms.z() - atoms.r() + probe_radius).minCoeff();

        float size = 0.5f * std::max(std::max(xmax - xmin, ymax - ymin), zmax - zmin);
        size_center_[0] = size;
        size_center_[1] = xmin + size;
        size_center_[2] = ymin + size;
        size_center_[3] = zmin + size;
	}
protected:
	virtual void initSizeCenter(const BaseCoordArray& atoms)
	{
	    this->oldBoundingBox(atoms, probe_radius_);
        float s = minsize_;
        end_level_ = 0;
        while (s < size_center_[0])
        {
            end_level_++;
            s *= 2.0f;
        }
        end_level_--;
        if(z_border_ == 0.0f && above_)
        {
            is_needed_ = true;
            touches_membrane_ = false;
            return;
        }
        else if(z_border_ == 0.0f && !above_)
        {
            is_needed_ = false;
            touches_membrane_ = false;
            return;
        }
		if(above_)
		{
			is_needed_ = size_center_[3] + size_center_[0] > z_border_; // bb reaches out of slab at top
			touches_membrane_ = is_needed_ && (size_center_[3] - size_center_[0] <= z_border_);
			if( is_needed_ && touches_membrane_)
			{
				size_center_[3] = z_border_ + size_center_[0];
			}
		}
		else
		{
			is_needed_ = size_center_[3] - size_center_[0] < z_border_; // bb reaches out of slab at bottom
			touches_membrane_ = is_needed_ && (size_center_[3] + size_center_[0] >= z_border_);
			if( is_needed_ && touches_membrane_)
			{
				size_center_[3] = z_border_ - size_center_[0];
			}
		}
	}
};


} // end of namespace

#endif /* MEMBRANESURFACE_H_ */
