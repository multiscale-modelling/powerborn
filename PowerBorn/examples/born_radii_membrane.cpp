/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */
#include "PowerBornMembrane.h"

void born_radii(int argc, char* argv[])
{
    if (argc < 2 || argc > 3)
    {
        printf(\
"Usage:\n\
 ./born_radii <input.pqr> [tension]\n\
\n\
Computes GB and NP energies for SLIM model.\n\
Output: gb_core  gb_head  np  total\n");
        abort();
    }
    std::string infile = argv[1];
    float tension = 0.030f;
    if (argc == 3)
        tension = atof(argv[2]);

    powerborn::PowerBornMembrane pbr3(11.0f, 15.0f);
    powerborn::BaseCoordArray atoms;

    powerborn::Array charges;
    atoms.parsePqr(infile, &charges);
    pbr3.update(atoms);
    float gb_core = pbr3.gbEnergy(atoms, charges, pbr3.getCoreRadii(), 2.0f,
            6.0f);
    float gb_head = pbr3.gbEnergy(atoms, charges, pbr3.getHeadRadii(), 6.0f,
            80.0f);
    float np = pbr3.npEnergy(atoms, tension);
    std::cout << gb_core << " " << gb_head << " " << np << " " << gb_core + gb_head + np << std::endl;
    /*const powerborn::Array& br_core = pbr1.getBornRadii();
    const powerborn::Array& br_head = pbr2.getBornRadii();
    for (unsigned int i = 0; i < br_core.size(); ++i)
    {
        std::cout << br_core[i] << " " << br_head[i] << std::endl;
    }
    openmp::Timer t1;
    t1.start();
    for (unsigned int i = 0; i < 1000; ++i)
        pbr3.update(atoms);
    t1.stop();
    t1.print();*/
}

int main(int argc, char* argv[])
{
    born_radii(argc, argv);
    return 0;
}

