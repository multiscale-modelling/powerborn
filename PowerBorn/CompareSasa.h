/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef COMPARESASA_H_
#define COMPARESASA_H_

#include "Typedefs.h"
#include "power_sasa.h"
#include "NumericSasa.h"

#include <math.h>


namespace powerborn
{

class CompareSasa
{
    POWERSASA::PowerSasa<float,Coord>* ps_;
    Diagram pd_;
    NumericSasa numeric_sasa_;
    std::vector<float> weights_, rmsd_, max_, total_sasa_new_, total_sasa_org_;
    std::vector<Coord> coords_;
    unsigned int count_;

    void dumpError(const BaseCoordArray& atoms, float max_err)
    {
        std::stringstream s;
        std::cout << "sasa comparison max error: " << max_err << " run: " << count_ << std::endl;
        s << "sasa_comp_err_" << count_ << ".dump";
        atoms.dump(s.str());
    }
public:
    CompareSasa(): ps_(NULL), count_(0)
    {
    }
    ~CompareSasa()
    {
        this->printResults();
        if(ps_) delete ps_;
    }
    void compare(const Array& new_sasa, const BaseCoordArray& atoms)
    {
        coords_.clear();
        weights_.clear();
        for(unsigned int i=0; i < atoms.size(); ++i)
        {
            coords_.push_back(powerborn::Coord(atoms.x()[i], atoms.y()[i], atoms.z()[i]));
            weights_.push_back(atoms.r()[i]);
        }
        if(! ps_)
        {
            ps_ = new POWERSASA::PowerSasa<float,powerborn::Coord>(coords_, weights_, 1, 0, 0, 0);
            ps_->update_coords(coords_, weights_);
            ps_->calc_sasa_all();
        }
        else
        {
            delete ps_;
            ps_ = new POWERSASA::PowerSasa<float,powerborn::Coord>(coords_, weights_, 1, 0, 0, 0);
            ps_->update_coords(coords_, weights_);
            ps_->calc_sasa_all();
        }
        const std::vector<float> & org_sasa = ps_->getSasa();
        float max = 0.0f;
        float sum=0.0f;
        float sumsq = 0.0f;
        float sasa_org = 0.0;
        float sasa_new = 0.0;
        for(unsigned int i=0; i<atoms.size(); ++i)
        {
            float os = org_sasa[i];
            float ns = new_sasa[i];
            float err = os - ns;
            max = fabs(err) > max ? fabs(err) : max;
            sum += err;
            sumsq += err * err;
            sasa_org += os;
            sasa_new += ns;
            if(fabs(err) > 1.0f)
            {
                PowerCell pc = pd_.createPowerCell(atoms, i);
                float num_sasa = numeric_sasa_.getSasa(pc);
                std::cout << "Sasa Error "<< i << " " << err << " "
                    << os << " " << ns << " "
                    <<  4.0f * M_PI * atoms.r()[i] * atoms.r()[i] <<
                    " num sasa " << num_sasa << std::endl;
                //ps_->calc_sasa_single(i, true);
                if(fabs(num_sasa - ns) > fabs(num_sasa - os))
                {
                    this->dumpError(atoms, (num_sasa - ns));
                }
            }
        }
        float rmsd = sqrt(sumsq / float(atoms.size()));
        rmsd_.push_back(rmsd);
        max_.push_back(max);
        total_sasa_new_.push_back(sasa_new);
        total_sasa_org_.push_back(sasa_org);
        count_++;
    }
    unsigned int count()
    {
        return count_;
    }
    void printResults()
    {
        for(unsigned int i=0; i<rmsd_.size(); ++i)
        {
            std::cout << i << " " << rmsd_[i] << " " << max_[i] << " " << total_sasa_org_[i] << " " << total_sasa_new_[i] << std::endl;
        }
    }
};

}

#endif /* COMPARESASA_H_ */
