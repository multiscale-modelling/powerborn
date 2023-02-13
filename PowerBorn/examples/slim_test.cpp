/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */
#include <time.h>

#include "PowerBornMembrane.h"

void born_radii(int argc, char* argv[])
{
	if(argc != 2) {
		printf("Usage:\n ./born_radii <input.pqr> \n\nComputes and prints Born radii and GB energies for ACC and FAST.\n");
		abort();
	}
	std::string infile = argv[1];
    powerborn::BaseCoordArray atoms;
    powerborn::Array charges;
    atoms.parsePqr(infile, &charges);

    // move structure to com
    float com_x = atoms.x().sum() / float(atoms.size());
    float com_y = atoms.y().sum() / float(atoms.size());
    float com_z = atoms.z().sum() / float(atoms.size());
    atoms.translate(-com_x, -com_y, -com_z);

    powerborn::PowerBornMembrane pbrm(0.0f, 15.0f);

    for(unsigned int i=0; i<1000; ++i)
    {
        pbrm.update(atoms);
        float core_energy = pbrm.gbEnergy(atoms, charges, pbrm.getCoreRadii(), 1.0f, 13.0f);
        float head_energy = pbrm.gbEnergy(atoms, charges, pbrm.getHeadRadii(), 13.0f, 80.0f);

        float z = atoms.z().sum() / float(atoms.size());
        std::cout << z << " " << core_energy + head_energy << std::endl;
        atoms.translate(0.0f, 0.0f, 0.1f);
    }
}

int main(int argc, char* argv[])
{
    born_radii(argc, argv);
	return 0;
}
