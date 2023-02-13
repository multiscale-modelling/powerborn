from scipy.optimize import minimize
import os
from subprocess import check_output
import numpy as np

class ParamOptimizer:

    def __init__(self, inputfile):
        self.files = []
        self.ref_energies = []

        with open(inputfile) as fin:
            for line in fin:
                self.files.append(line.split()[0])
                self.ref_energies.append(float(line.split()[1]))
        self.ref_energies = np.asanyarray(self.ref_energies)
        
        
    def func(self, p):
        command = ["gb_energy_gpu"]
        command.extend([str(x) for x in p])
        command.extend(self.files)
        output = check_output(command)
        #print output
        result = np.asanyarray([float(x.split()[-1]) for x in output.splitlines() if x.startswith("GB ENERGY")])
        result -= self.ref_energies
        result /= self.ref_energies
        tmp = np.sqrt((result * result).mean())
        print p, tmp
        return tmp 
    
    def optimize(self):
        start = [1.0, 0.0]
        res = minimize(self.func, start, method='Powell')
        print res.x, res.fun

if __name__ == "__main__":
    import sys
    infile = sys.argv[1]
    p = ParamOptimizer(infile)
    if len(sys.argv) == 4:
        # run with given params
        param1 = float(sys.argv[2])
        param2 = float(sys.argv[3])
        print param1, param2, p.func([param1, param2])
    else:
        p.optimize()
