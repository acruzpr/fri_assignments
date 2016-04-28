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
                         temperature = 8000 * ase.units.kB,
                         fmax = 0.01,
                         dr = displacement,
                         adjust_cm = True,
                         minenergy = -108.315,
                         #minenergy = -1,
                         distribution = 'uniform',
                         optimizer_logfile = None # uncomment this line if you do not want the information from each optimization to print 

                         )
    opt.run(5000)
    for i in range(len(opt.norm_dists)):
        print i, ' ', opt.norm_dists[i]

        figure()
        hist(opt.norm_dists[i], 50, normed=1)
        xlabel('Distance from center for ' + str(i))
        ylabel('Prob. Density')
        savefig('histo_38_' + str(i) + '.png')

        figure()
        plot(opt.norm_dists[i])
        xlabel('Iteration')
        ylabel('Distance from center for ' + str(i))
        savefig('line_38_' + str(i) + '.png')

tcalls=0
for i in range(1): 
    p = tsase.io.read_con('cluster_38.con')

    lj = tsase.calculators.lj(cutoff=35.0)
    p.center(100.0)
    p.set_calculator(lj)

    basin_hop(p)

    print 'force calls: ', lj.force_calls
    print 'pe: ', p.get_potential_energy()
    print 'distances: ', p.get_all_distances()

    tcalls+=lj.force_calls

print "avg: ", tcalls/10

sys.exit()

