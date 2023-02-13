#!/usr/bin/python

import sys
import getopt
from show_sasa import show_sasa
from show_cgo import show_cgo
import time

# start pymol
import pymol
pymol.finish_launching()
time.sleep(1.0)

def show_pdb(file):
    pymol.cmd.load(inputfile,format='pdb')

if __name__=="__main__":
    opts,ar = getopt.getopt(sys.argv[1:], "s:c:p:v:")
    print opts
    for a in opts:
        if a[0]=='-s':
            print "sasa: ", a[1]
            show_sasa(a[1], 1.4)
        elif a[0]=='-c':
            print "cgo: ", a[1]
            show_cgo(a[1])
        elif a[0]=='-p':
            print "pdb: ", a[1]
            show_pdb(a[1])
        elif a[0]=='-v':
            print "vdw: ", a[1]
            show_sasa(a[1], 0.0)

