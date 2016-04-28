import numpy as np
import sys
from pylab import *
from math import *
import tsase
import ase
from basin import *
import scipy
import time

class BCM():
    @classmethod
    def mag(cls,vec):
        return sqrt(dot(vec,vec))

    @classmethod
    def q_avg(cls,l,m,positions,cm):
        sum = complex(0,0)
        bonds = 0
        for a in range(len(positions)):
            b = a+1
            while b<len(positions):
            #for b in range(len(positions)):
                dist =  cls.mag(positions[a]-positions[b])
                if dist<1.35 and a!=b:
#            	print positions[a]
#            	print positions[b]
#            	print cm
                    rvec = (positions[a]+positions[b]-cm-cm)/2
                    #rvec = (positions[a]-positions[b])
                    #rvec = (positions[a]+positions[b])/2
                    #if rvec gives [x,y,z]
                    theta = atan2(rvec[1],rvec[0])
                    phi = atan2(rvec[2],hypot(rvec[0],rvec[1])) + pi/2
                    sh = scipy.special.sph_harm(m,l,theta,phi)
                    sum += sh
                    bonds += 1
                    #print "sh: ", sh
                b += 1
        #print "sum: ", sum
        #print "bonds: ",bonds
        return sum / bonds

    @classmethod
    def q_norm(cls,l,positions,cm):
        sum = 0
        m = -l
        while m <= l:
            qavg = abs(cls.q_avg(l,m,positions,cm))
            sum += qavg * qavg
            m+=1
        #print sum
        return sqrt(4*pi*sum/(2*l+1))

    @classmethod
    def get_bcm(cls, positions,cm):
        #print 'cm: ',cm
        ans = np.array([0.0]*6)
        for i in range(6):
            #print i
            ans[i] = cls.q_norm(i*2,positions,cm)
        return ans

    @classmethod
    def cm(cls, positions):
        ans = np.array([0]*3)
        for i in positions:
            ans += i
        ans /= len(positions)
        return ans




#print 'start'
#print BCM()
#arr = [0]*100
#for i in range(100):
#    arr[i] = tsase.io.read_con('../100_lj38/'+str(i)+'.con');
#    print i
#    #p = tsase.io.read_con('answer_38.con')
#    #bcm = [0]
#
#start = time.time()
#for i in range(100):
#    p = arr[i]
#    bcm = BCM().get_bcm(p.get_positions(),cm(p.get_positions()))
#    print bcm
#    print sqrt(dot(bcm,bcm))
#print 'time: ', time.time()-start
#
#
#sys.exit()


