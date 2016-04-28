#!/usr/bin/env python

import numpy
import sys
import matplotlib
matplotlib.use("agg") # Change to "agg" to run on FRI
from pylab import *
import tsase
import ase
from basin import *

##### Setting up the atoms object with ASE and TSASE #########
##############################################################



def basin_hop(p):
    displacement = 0.3
    opt = BasinHopping(  atoms = p, 
                         temperature = 3000 * ase.units.kB,
                         fmax = 0.01,
                         dr = displacement,
                         adjust_cm = True,
                         minenergy = -173.927,
                         distribution = 'uniform',
                         optimizer_logfile = None # uncomment this line if you do not want the information from each optimization to print 

                         )
    opt.run(5000)
    return opt.steps

tcalls=0
failCalls=0
sum_opts=0
for i in range(10): 
    p = tsase.io.read_con('../10_lj38/'+str(i)+'.con')

    lj = tsase.calculators.lj(cutoff=35.0)
    p.center(100.0)
    p.set_calculator(lj)

    opts = basin_hop(p)
    sum_opts+=opts

    print 'force calls: ', lj.force_calls
    print 'opts: ', opts
    print 'pe: ', p.get_potential_energy()

    if(lj.force_calls!=5000):
        tcalls+=lj.force_calls
    else:
    	failCalls+=1

print "avg: ", tcalls/(10-failCalls)
print "avg opts: ", sum_opts/10

sys.exit()













