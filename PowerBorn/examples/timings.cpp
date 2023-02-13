/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#include "PowerBornMS.h"

void timing(const std::string& pqrfile)
{
	powerborn::PowerBornMS acc = powerborn::PowerBornMS::getAcc();
	powerborn::PowerBornMS fast = powerborn::PowerBornMS::getFast();
    powerborn::CoordArray c;
    c.parsePqr(pqrfile);
    powerborn::PowerDiagram pd = acc.getPowerDiagram(c);
    openmp::Timer t("");
    t.start();
    for(unsigned int i=0; i<1000; ++i) acc.update(c, &pd);
    t.stop();
    printf("ACC: ");
    t.print();
    t.reset();
    t.start();
    for(unsigned int i=0; i<1000; ++i) fast.update(c, &pd);
    t.stop();
    printf("FAST: ");
    t.print();
}

void usage() {
	printf("Usage:\n./timings input1.pqr\n\nMeasure time of 1000 PowerBorn computations for a given pqr\n");
	abort();
}

int main(int argc, char* argv[])
{
	if(argc != 2) usage();
	std::string fname(argv[1]);
    timing(fname);
	return 0;
}

