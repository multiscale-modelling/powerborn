/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */


#include <deque>
#include <Eigen/StdVector>
#include <Eigen/StdDeque>

#include "PowerDiagram.h"
#include "Atoms.h"
#include "PowerBornMS.h"
#include "PowerBornMS2.h"
#include "power_sasa.h"
#include "CompareSasa.h"

void compare_powersasa(powerborn::CoordArray atoms, unsigned int);

void debug_pd(const std::string dump, unsigned int aid)
{
    powerborn::CoordArray atoms;
    atoms.parseDump(dump);
    //atoms.writePqr("debug_atoms.pqr", 0.0);
    /*powerborn::Diagram d(aid);
    d.build(atoms);
    float s1 = d.getSasa().sum();
    d.build(atoms);
    float s2 = d.getSasa().sum();
    if(s1 != s2)
    {
        std::cout << "repetition Test failed!!! " << s1 << " " << s2 << std::endl;
    }*/
    compare_powersasa(atoms, aid);
}

void test_pd(const std::string& pqrfile)
{
	powerborn::CoordArray atoms;
	atoms.parsePqr(pqrfile);

	{
	    openmp::Timer t("new PD & PS");
	    t.start();
	    powerborn::Diagram d;
#pragma omp parallel
	    for(unsigned int i=0; i<1000; ++i)
	    {
            d.build(atoms);
#pragma omp barrier
	    }
	    t.stop();
	    t.print();
	}
}

void test_old_pd(const std::string& pqrfile)
{
    powerborn::CoordArray atoms;
    atoms.parsePqr(pqrfile);
    std::vector<powerborn::Coord> coords;
    std::vector<float> weights;
    for(unsigned int i=0; i < atoms.size(); ++i)
    {
        coords.push_back(powerborn::Coord(atoms.x()[i], atoms.y()[i], atoms.z()[i]));
        weights.push_back(atoms.r()[i]);
    }
    POWERSASA::PowerSasa<float,powerborn::Coord> ps(coords, weights, 1, 0, 0, 0);
    {
        openmp::Timer t("old PD & PS");
        t.start();
        for(unsigned int i=0; i<1000; ++i)
        {
            ps.update_coords(coords, weights);
            ps.calc_sasa_all();
        }
        t.stop();
        t.print();
    }
}

void test_pbr(const std::string& pqrfile)
{
    powerborn::CoordArray atoms;
    atoms.parsePqr(pqrfile);
    openmp::Timer t1("new pbr");
    t1.start();
    powerborn::PowerBornMS2 pbr2 = powerborn::PowerBornMS2::getAcc();
    for(unsigned int i=0; i<1000; ++i)
    {
        pbr2.update(atoms);
    }
    pbr2.timings();
    t1.stop();
    t1.print();
    
    openmp::Timer t2("old pbr");
    t2.start();
    powerborn::PowerBornMS pbr = powerborn::PowerBornMS::getAcc();
    powerborn::PowerDiagram pd =  pbr.getPowerDiagram(atoms);
    for(unsigned int i=0; i<1000; ++i)
    {
        pbr.update(atoms, &pd);
    }
    pbr.timings();
    t2.stop();
    t2.print();
}

void compare_powersasa(const std::string pqrfile)
{
    powerborn::CoordArray atoms;
    atoms.parsePqr(pqrfile);
    compare_powersasa(atoms, 1000000);
}

void compare_powersasa(powerborn::CoordArray atoms, unsigned int aid)
{
    powerborn::Diagram d(aid);
    powerborn::CompareSasa comparison;
    for(unsigned int j=0; j<1; ++j)
    {
        d.build(atoms);
        comparison.compare(d.getSasa(), atoms);
    }
}

void usage() {
	std::cout << "Usage:" << std::endl << "./test_powerdiagram input1.pqr [input2.pqr] ..." << std::endl;
    std::cout << "Usage for debugging:" << std::endl << "./test_powerdiagram input1.dump atom_id" << std::endl;

	abort();
}

int main(int argc, char* argv[])
{
	if(argc == 1) usage();
	std::string fname(argv[1]);
	if (fname.substr(fname.size()-3,3) == std::string("pqr"))
	{
		for(int i=1; i<argc; ++i)
		{
            compare_powersasa(argv[i]);
			test_pd(argv[i]);
			test_old_pd(argv[i]);
			test_pbr(argv[i]);
		}
	}
	else if (fname.substr(fname.size()-4,4) == std::string("dump"))
	{
	    unsigned int aid = atoi(argv[2]);
        debug_pd(fname, aid);
	}
	else
	{
		usage();
	}
	return 0;
}

