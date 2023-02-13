#!/usr/bin/python
import time
import sys
import pymol


def show_cgo(inputfile):
    #load cgo from file
    data=[]
    with open(inputfile) as fin:
        for line in fin:
            data+=map(lambda x: float(x), line.split())
    pymol.cmd.load_cgo(data, inputfile, 1)

if __name__=='__main__':
# start pymol
    pymol.finish_launching()
    time.sleep(1.0)
    try:
        for file in sys.argv[2:]:
            show_cgo(file)
    except:
        print "usage: ./show_cgo.py  pdbfile.pdb  cgo1.cgo cgo2.cgo ..."
        raise
    try:
        pymol.cmd.load(sys.argv[1],sys.argv[1])
        pymol.cmd.show('surface')
    except:
        pass

