import numpy as np

from ase.optimize.optimize import Dynamics
from ase.optimize.fire import FIRE
from tsase.optimize.sdlbfgs import SDLBFGS
from ase.units import kB
from ase.parallel import world
from ase.io.trajectory import PickleTrajectory
import tsase
import sys
from bcm import BCM
import traceback
import copy
import random
import math


class BasinHopping(Dynamics):
    """Basin hopping algorithm.

    After Wales and Doye, J. Phys. Chem. A, vol 101 (1997) 5111-5116

    and 

    David J. Wales and Harold A. Scheraga, Science, Vol. 285, 1368 (1999)
    """

    def __init__(self, atoms,
                 temperature=100 * kB,
                 optimizer=SDLBFGS,
                 fmax=0.1,
                 dr=0.1,
                 logfile='-', 
                 trajectory=None,
                 optimizer_logfile='-',
                 local_minima_trajectory='local_minima.con',
                 mss=0.2,
                 minenergy=None,
                 adjust_step_size=None,
                 adjust_every = None,
                 target_ratio = 0.5,
                 adjust_fraction = 0.05,
                 pushapart = 0.4,
                 jumpmax=2,
                 inertia_weight=.4,
                 method='PSO'
                 ):
        Dynamics.__init__(self, atoms, logfile, trajectory)
        self.local_optimizations = 0
        self.bcms = 0
        self.kT = temperature
        self.optimizer = optimizer
        self.fmax = fmax
        self.dr = dr
        self.method = method

        self.cm = atoms.get_center_of_mass()

        self.optimizer_logfile = optimizer_logfile
        self.lm_trajectory = local_minima_trajectory
        self.minenergy = minenergy
        self.energy = 0
        self.adjust_step = adjust_step_size
        self.adjust_every = adjust_every
        self.target_ratio = target_ratio
        self.adjust_fraction = adjust_fraction 
        self.pushapart = pushapart
        self.jumpmax = jumpmax
        self.mss = mss
        self.inertia_weight = inertia_weight
        self.velocity = 0.0 * self.atoms.get_positions()
        self.initialize()


    def initialize(self):
        self.positions = 0.0 * self.atoms.get_positions()
        pos = self.push_apart(self.atoms.get_positions())
        self.Emin = self.get_energy(pos) or 1.e32
        self.rmin = pos
        #self.positions = self.atoms.get_positions()
        self.call_observers()
        self.bcm_calculator = BCM
        self.bcm_changed = True
        self.log(-1, self.Emin, self.Emin,self.dr)
        if isinstance(self.lm_trajectory, str):
            tsase.io.write_con(self.lm_trajectory,self.atoms,w='a')
    def run(self, steps, bests):
        """Hop the basins for defined number of steps."""
        self.steps = 0
        ro = self.positions
        Eo = self.get_energy(ro)
        acceptnum = 0
        rejectnum = 0
        #for step in range(steps):
        while acceptnum is 0:
            En = None
            self.steps += 1
            while En is None:
                En,rn,velocity = self.move(ro, bests)
            if En < self.Emin:
                self.Emin = En
                self.rmin = self.atoms.get_positions()
                self.call_observers()
            self.log(self.steps, En, self.Emin,self.dr)
            if Eo >= En:
                accept = True
            else:
                accept = np.exp((Eo - En) / self.kT) > np.random.uniform()
            if rejectnum > self.jumpmax:
                acceptnum = 1
            print 'Accept: ', accept
            if accept:
                acceptnum += 1.
                rejectnum = 0
                ro = rn.copy()
                self.velocity = velocity
                Eo = En
                self.bcm_changed = True
                if self.lm_trajectory is not None:
                    tsase.io.write_con(self.lm_trajectory,self.atoms,w='a')
            else:
                rejectnum += 1
            if self.minenergy != None:
                if Eo < self.minenergy:
                    break
            if self.adjust_step == True:
                if step % self.adjust_every == 0:
                    ratio = float(acceptnum)/float(self.adjust_every)
                    if ratio > self.target_ratio:
                       self.dr = self.dr * (1+self.adjust_fraction)
                    elif ratio < self.target_ratio:
                       self.dr = self.dr * (1-self.adjust_fraction)

    def log(self, step, En, Emin,dr):
        if self.logfile is None:
            return
        name = self.__class__.__name__
        self.logfile.write('%s: step %d, energy %15.6f, emin %15.6f, dr %15.6f\n'
                           % (name, step, En, Emin,dr))
        self.logfile.flush()

    def __str__(self):
        return 'energy %15.6f, emin %15.6f, dr %15.6f' % (self.get_energy(), self.Emin, self.dr)

    def __repr__(self):
        return self.__str__()

    def move(self, ro, bests):
        """Move atoms by a random step."""
        atoms = self.atoms
        velocity = copy.deepcopy(self.velocity)

        do_pso = True
        if self.method == 'firefly':
            do_pso = False
        elif self.method == 'PSO_split':
            best_to_use = bests[-1]
            best_dist = 1e32
            for i in bests:
                d = i.get_bcm() - self.get_bcm()
                di = np.sqrt(np.vdot(d,d))
                if(di<best_dist):
                    best_to_use = i
                    best_dist = di

            do_pso = True
        else:
            do_pso = True
            best_to_use = bests[-1]

        # PSO heuristic
        if (not do_pso):
            r_arr = [None] * len(bests)
            r_sum = 0
            pe_sum = 0
            for i in range(len(bests)):
                    d = self.get_bcm() - bests[i].get_bcm()
                    r_arr[i] = np.sqrt(np.vdot(d,d))
                    r_sum += r_arr[i]
                    pe_sum += bests[i].get_energy()
            r_mean = r_sum / len(bests)
            pe_mean = abs(pe_sum / len(bests))
            denom = math.sqrt(r_mean*pe_mean)
            attract_arr = [None] * len(bests)
            for i in range(len(bests)):
                    numer = bests[i].get_energy() * r_arr[i]
                    numer *= numer
                    attract_arr[i] = max(1e-25,math.e**(-1*numer/denom))
            sort_attract = sorted(attract_arr)
            print "sorted: ", sort_attract

            velocity = self.inertia_weight*velocity
            E = np.random.normal(1,max(0.01,.01*min(25,sort_attract[-1]/sort_attract[-2])), size=(len(atoms),3))
            for i in range(len(bests)):
                    velocity += attract_arr[i]*(bests[i].positions-ro)*E/sort_attract[-1]
        else:
            print "using best: ", best_to_use.get_bcm()
            for i in range(len(ro)):
                for j in range(3):
                    velocity[i][j] = self.inertia_weight*velocity[i][j] + random.uniform(0,2)*(self.rmin[i][j]-ro[i][j]) + random.uniform(0,2)*(best_to_use.positions[i][j]-ro[i][j])

        rn = ro + velocity
        rn = self.push_apart(rn)
        atoms.set_positions(rn)
        En = self.get_energy(rn)
        cm = atoms.get_center_of_mass()
        atoms.translate(self.cm - cm)

        rn = atoms.get_positions()
        world.broadcast(rn, 0)

        atoms.set_positions(rn)
        return En,atoms.get_positions(),velocity

    def get_minimum(self):
        """Return minimal energy and configuration."""
        atoms = self.atoms.copy()
        atoms.set_positions(self.rmin)
        return self.Emin, atoms

    def get_energy(self, positions=None):
        #print 'getting energy, pos: ', positions
        if  positions is None:
            #print 'energy is: ', self.energy
            return self.energy
        """Return the energy of the nearest local minimum."""
        if np.sometrue(self.positions != positions):
            self.positions = positions
            self.atoms.set_positions(positions)
 
            try:
                opt = self.optimizer(self.atoms, 
                                     logfile=self.optimizer_logfile,
                                     maxstep=self.mss)
                opt.run(fmax=self.fmax)
                self.energy = self.atoms.get_potential_energy()
                #print 'energy good ', self.energy
                self.local_optimizations += 1
            except:
                print 'Exception: ', sys.exc_info()[0]
                traceback.print_exc()
                # Something went wrong.
                # In GPAW the atoms are probably to near to each other.
                return None
            
        #print 'old energy: ', self.energy
        return self.energy

    def push_apart(self,positions):
        movea = np.zeros(np.shape(positions))
        alpha = 0.025
        for w in range(500):
            moved = 0
            movea = np.zeros(np.shape(positions))
            for i in range(len(positions)):
                for j in range(i+1,len(positions)):
                    d = positions[i] - positions[j]
                    magd = np.sqrt(np.vdot(d,d))
                    if magd < self.pushapart and magd > 0:
                        moved += 1
                        vec = d/magd
                        movea[i] += alpha *vec
                        movea[j] -= alpha *vec
            positions += movea
            if moved == 0:
                break
        return positions

    def get_dist_geo_center(self):
        position = self.atoms.get_positions()
        geocenter = np.sum(position,axis=0)/float(len(position))
        distance = np.zeros(len(position))
        for i in range(len(distance)):
            vec = position[i]-geocenter
            distance[i] = np.sqrt(np.vdot(vec,vec))
        distance /= np.max(distance)  #np.sqrt(np.vdot(distance,distance))
        return distance 

    def get_bcm(self):
        if self.bcm_changed:
            self.bcms += 1
            self.bcm_changed = False
            self.bcm = self.bcm_calculator.get_bcm(self.atoms.get_positions(), self.cm)
        return self.bcm
    def copy(self):
        return BasinCopy(self.positions, self.get_energy(), self.get_bcm())


class BasinCopy:
    def __init__(self,
            positions,
            energy,
            bcm
            ):
        self.positions = copy.deepcopy(positions)
        self.bcm = copy.copy(bcm)
        self.energy = energy
    def get_energy(self):
        return self.energy
    def get_bcm(self):
        return self.bcm

