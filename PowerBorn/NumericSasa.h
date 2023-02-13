/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef NUMERICSASA_H_
#define NUMERICSASA_H_

#include "PowerCell.h"
#include "Atoms.h"

namespace powerborn
{

class NumericSasa: private SasaPoints2
{
    void initGrid(unsigned int grid_points)
    {
        grid_.clear();
        float phi_step = 2.0f * M_PI / float(grid_points);
        float theta_step = M_PI / float(grid_points);
        float phi = 0.0f;
        float theta = 0.5f * theta_step;
        for(unsigned int i=0; i<grid_points; ++i)
        {
            for(unsigned int j=0; j<grid_points; ++j)
            {
                Coord pos;
                pos[0] = cos(phi) * sin(theta);
                pos[1] = sin(phi) * sin(theta);
                pos[2] = cos(theta);
                float weight = phi_step * theta_step * fabs(sin(theta));
                grid_.insert(pos[0], pos[1], pos[2], weight);
                phi += phi_step;
            }
            phi = i % 2 == 1 ? 0.0f : 0.5f * phi_step;
            theta += theta_step;
        }
        float sum = grid_.r().sum() / (4.0f * M_PI);
        std::cout << "theta final " << theta + theta_step << std::endl;
        std::cout << "normalization " << sum << std::endl;
        grid_.r() /= sum;
        //grid_.r().setConstant(0.01f);
        //grid_.writePqr("NumericSasaGrid.pqr", 0.0f);
    }

    CoordArray grid_, actual_grid_;
public:
    NumericSasa(unsigned int grid_points=200)
    {
        this->initGrid(grid_points);
    }
    float getSasa(const PowerCell& pc)
    {

        // TODO test this
        actual_grid_.getScaledCoords(grid_, pc.normedRadius());
        for(unsigned int i=0; i<pc.neighbourSize(); ++i)
        {
            Coord pnb = pc.position(i);
            float nbrad = pc.radius(i);
            actual_grid_.checkSphere(pnb[0], pnb[1], pnb[2], nbrad);
        }
        float sasa = actual_grid_.r().sum() * pc.radius() * pc.radius() * pc.normedRadius() * pc.normedRadius();
        //actual_grid_.r().setConstant(0.01f);
        //actual_grid_.writePqr("test_num_sasa.pqr", 0.0f);
        return sasa;
    }
};


}

#endif // NUMERICSASA_H_
