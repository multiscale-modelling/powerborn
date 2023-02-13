/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#include "PowerBornMS2.h"

int main(int argc, char* argv[])
{
    if(argc < 8)
    {
        std::cout << "Usage: " << std::endl << std::endl
                << "gb_energy <fit1> <fit2> <alpha> <beta> <e_in> <e_out> [pqr1.pqr] [energy1] ..."
                << std::endl;
    }
    float fit1 = atof(argv[1]);
    float fit2 = atof(argv[2]);
    float alpha = atof(argv[3]);
    float beta = atof(argv[4]);
    float e_in = atof(argv[5]);
    float e_out = atof(argv[6]);
    powerborn::PowerBornMS2 pbr(0.25f, 1.4f, 0.8f, 10.0f, fit1, fit2, alpha, beta);
    double sumsq = 0;
    unsigned int count = 0;
    for(int i=7; i<argc; i+=2)
    {
        std::string pqrfile = argv[i];
        double ref = atof(argv[i+1]);
        pbr.updatePqr(pqrfile);
        double energy = pbr.gbEnergy(e_in, e_out);
        sumsq += (energy - ref) * (energy - ref) / (ref * ref);
        count++;
        std::cout << pqrfile << " " << ref << " " << energy << std::endl;
    }
    std::cout << sqrt(sumsq / double(count)) << std::endl;
	return 0;
}

