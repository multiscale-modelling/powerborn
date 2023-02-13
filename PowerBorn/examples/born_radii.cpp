/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */
#include <time.h>
#include <boost/algorithm/string/predicate.hpp>

#include "PowerBornMS.h"
#include "PowerBornMS2.h"


void born_radii(int argc, char* argv[])
{
    if(argc != 3) {
        printf("Usage:\n ./born_radii <input1.pqr> <input2.pqr>\n\nComputes and prints Born radii and GB energies for ACC and FAST.\n");
        abort();
    }
    std::string infile1 = argv[1];
    std::string infile2 = argv[2];

    powerborn::PowerBornMS pbr1 = powerborn::PowerBornMS::getAcc();
    powerborn::PowerBornMS2 pbr2 = powerborn::PowerBornMS2::getAcc();
    powerborn::PowerBornMS2 pbr3 = powerborn::PowerBornMS2::getAcc();
    pbr3.setProbeRadius(1.4f);
    powerborn::CoordArray atoms1;
    powerborn::CoordArray atoms2;
    if (boost::ends_with(infile1, ".dump"))
    {
        atoms1.parseDump(infile1);
    }
    else
    {
        atoms1.parsePqr(infile1);
    }
    if (boost::ends_with(infile2, ".dump"))
    {
        atoms2.parseDump(infile2);
    }
    else
    {
        atoms2.parsePqr(infile2);
    }
    powerborn::Array charges;
    charges.setConstant(atoms1.size(), 0.25f);

    {
        powerborn::PowerDiagram pd1 = pbr1.getPowerDiagram(atoms1);
        pbr1.update(atoms1, &pd1);
        pbr2.update(atoms1);
        pbr3.update(atoms1);
        const powerborn::Array& br1 = pbr1.getBornRadii();
        const powerborn::Array& br2 = pbr2.getBornRadii();
        const powerborn::Array& br3 = pbr3.getBornRadii();
        printf("# MS    MS2\n");
        for(unsigned int i=0; i< br1.size(); ++i)
        {
            printf("%10.6f  %10.6f  %10.6f  %10.6f\n",br1[i], br2[i], br3[i], br2[i] - br1[i]);
        }
        printf("# GB Energies: %f  %f %f\n", pbr1.gbEnergy(atoms1, charges), pbr2.gbEnergy(atoms1, charges), pbr3.gbEnergy(atoms1, charges));
    }
    {
        powerborn::PowerDiagram pd2 = pbr1.getPowerDiagram(atoms2);
        pbr1.update(atoms2, &pd2);
        pbr2.update(atoms2);
        pbr3.update(atoms2);
        const powerborn::Array& br1 = pbr1.getBornRadii();
        const powerborn::Array& br2 = pbr2.getBornRadii();
        const powerborn::Array& br3 = pbr3.getBornRadii();
        printf("# MS    MS2\n");
        for(unsigned int i=0; i< br1.size(); ++i)
        {
            printf("%10.6f  %10.6f  %10.6f  %10.6f\n",br1[i], br2[i], br3[i], br2[i] - br1[i]);
        }
        printf("# GB Energies: %f  %f %f\n", pbr1.gbEnergy(atoms2, charges), pbr2.gbEnergy(atoms2, charges), pbr3.gbEnergy(atoms2, charges));
    }
}

int main(int argc, char* argv[])
{
    born_radii(argc, argv);
    return 0;
}
