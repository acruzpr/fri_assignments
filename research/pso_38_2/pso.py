#!/usr/bin/env python

import numpy
import sys
import matplotlib
matplotlib.use("agg") # Change to "agg" to run on FRI
from pylab import *
from math import *
import tsase
import ase
from ase import Atoms
from basin import BasinHopping
import scipy
import time
from bcm import BCM
import os

#Constants
POPULATION = 30
NUM_BESTS = 5
MIN_ENERGY = -173.928
#MIN_ENERGY = -108.315
BEST_BCM_CUTOFF = -1
REMOVE_FOR_DIVERSITY = max(1,int(POPULATION / 5))
STRUCTURE_SIZE = 38


def mag(cls,vec):
    return sqrt(dot(vec,vec))

def get_basin_hop(p,energy_goal):
    displacement = 0.3
    opt = BasinHopping(  atoms = p, 
                         temperature = 4000 * ase.units.kB,
                         fmax = 0.01,
                         dr = displacement,
                         minenergy = energy_goal,
                         optimizer_logfile = None # uncomment this line if you do not want the information from each optimization to print 

                         )
    return opt

lj = tsase.calculators.lj(cutoff=35.0)
def random_structure():
    #r = tsase.io.read_con(random_structure_arr[int(random()*100)])

    r = tsase.io.read_con(random_structure_arr[0])
    r = r[0:STRUCTURE_SIZE]
    r.set_positions(r.get_positions()+np.random.uniform(-1,1,(STRUCTURE_SIZE,3)))

    #r = tsase.io.read_con("structure.con")
    #r.set_positions(r.get_positions()[0:25])
    r.center(100.0)
    r.set_calculator(lj)
    ret = get_basin_hop(r,MIN_ENERGY)
    return ret

class placeholder_structure:
    def get_bcm(self):
        return [0] * 6
    def get_energy(self):
        return 0

#inserts into the list of bests if better
def insert_best(bests, structure):
    nenergy = structure.get_energy()
    #print 'nenergy: ', nenergy
    if (bests[0].get_energy()>nenergy):

        tsase.io.write_con("best.con",structure.atoms,w='a')
        i=0
        a=0
        within_cutoff=None
        while (a < len(bests)):
            if bests[a].get_energy()>nenergy:
                i = a
                #print "i: ", i
            temp_mag = structure.get_bcm()-bests[a].get_bcm()
            mag_mag = mag(temp_mag, temp_mag)
            print "mag_mag: ", mag_mag
            if (within_cutoff is None and mag_mag < BEST_BCM_CUTOFF):
                #print 'within_cutoff: ', temp_mag
                within_cutoff = a
            a+=1
            #if(a<len(bests)):
                #print 'energy: ', bests[a].get_energy()
        #i+=1
        #print i, within_cutoff
        if (within_cutoff is not None):
            if (bests[within_cutoff].get_energy()>nenergy):
                scopy = structure.copy()
                bests[within_cutoff] = scopy
        
                #while (within_cutoff < i and within_cutoff < len(bests)-1):
                #    within_cutoff += 1
                #    bests[within_cutoff-1] = bests[within_cutoff]
                #    bests[within_cutoff] = scopy
        else:
            temp2 = structure.copy()
            bests[i] = temp2
           # while i>0:
           #     i-=1
           #     temp = bests[i]
           #     bests[i] = temp2
           #     temp2 = temp
        return sorted(bests, key=lambda at: at.get_energy(), reverse=True)
    return bests

def insert_all(bests, structures):
    for i in range(len(structures)):
        bests = insert_best(bests, structures[i])
    print 'current bests:'
    for i in range(NUM_BESTS):
        print i, ': ', bests[i].get_energy(), " bcm: ", bests[i].get_bcm()
    return bests

def remove_worst(structures):
    #print 'remove:'
    struct = sorted(structures, key=lambda at: at.get_energy(), reverse=True)
    #for i in struct:
        #print i.get_energy()
    #print struct
    sum_opts = 0
    sum_bcms = 0
    for i in range(REMOVE_FOR_DIVERSITY):
        sum_opts += struct[i].local_optimizations
        sum_bcms += struct[i].bcms
        struct[i] = random_structure()
    #print struct
    return struct,sum_opts,bcms


print 'start'

if os.path.exists("best.con"):
	os.remove("best.con")
if os.path.exists("local_minima.con"):
	os.remove("local_minima.con")

random_structure_arr = [0]*100
for i in range(100):
    random_structure_arr[i] = '../100_lj38/'+str(i)+'.con'

structures = [None] * POPULATION
for i in range(POPULATION):
    structures[i] = random_structure()
    #print i, ': ',  structures[i].get_bcm()


bests = [placeholder_structure()] * NUM_BESTS

best_energy = 1.e32
local_optimizations = 0
bcms = 0
bests = insert_all(bests,structures)
while best_energy>MIN_ENERGY:
    for i in structures:
        i.run(1,bests)
    structures,remove_opts,remove_bcms = remove_worst(structures)
    bests = insert_all(bests,structures)
    local_optimizations += remove_opts
    bcms += remove_bcms
    best_energy = min(bests[-1].get_energy(),best_energy)
    print "Optimizations: ", local_optimizations, "Bcms: ", bcms, "Best Energy: ", best_energy

for i in structures:
    local_optimizations += i.local_optimizations
    bcms += i.bcms
print "Optimizations: ", local_optimizations, "Bcms: ", bcms, "Best Energy: ", best_energy

for i in range(NUM_BESTS):
    print i, ': ', bests[i].get_energy(), " bcm: ", bests[i].get_bcm()
    #tsase.io.write_con("bests.con",



#sys.exit()


