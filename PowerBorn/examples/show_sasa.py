#!/usr/bin/python
#    This file is part of the PowerBorn Library
#    Please refer to the file LICENSE for License
#    and Copyright information.
import time
import sys
import pymol



def show_sasa(inputfile, probe_radius):
    radii=[]
    with open(inputfile) as fin:
        for line in fin:
            if not line.startswith("ATOM"): continue
            try:
                radii.append(float(line.split()[-1]))
            except:
                radii.append(float(line.split()[-2])) # if the last entry is the atom element
    #load the pdb in pymol
    pymol.cmd.load(inputfile,inputfile,format='pqr')
    # alter the vdw spheres in pymol
    for a in range(len(radii)):
        pymol.cmd.alter("%s&id %s" % (inputfile,a+1), "vdw=%s" % (radii[a]+probe_radius))
    pymol.cmd.show("spheres")

if __name__=='__main__':
    pymol.finish_launching()
    time.sleep(1.0)
    show_sasa(sys.argv[1],float(sys.argv[2]))


