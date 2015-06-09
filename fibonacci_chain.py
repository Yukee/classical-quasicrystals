# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 16:38:59 2014

@author: nicolas

At short times the Fibonacci spring chain behaves as a periodic spin chain,
whose spings have a stiffness equal to the harmonic mean of the Fibonacci sprins stiffness.
TODO: behaviour at large times? Is there energy transfer to the higher modes in the Fibonacci case?

Here we evolve the solution in time using the leapfrog integration method. 
This method is stable for an oscillatory motion, as long as dt < 2/nu, where nu is any eigenfrequency of the system.
"""

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import hmean

#compute Fibonacci numbers
def fib(n):
    a, b = 0, 1
    for i in range(n):
        a, b = b, a + b
    return a

# total time of the simulation
Ntime = 6000
n = 13
# length of the spring chain
L = fib(n)-1
# Fibonacci numbers (used to generate the sequence of couplings)
p = fib(n-2)
q = fib(n-1)
# quasiperiodic spring stiffness
ks = 6.
kw = 1.
# return kw or ks according to the Fibonacci sequence
def k(i):
    if np.remainder(p*i,p+q) >= p:
        ktemp = kw
    else: ktemp = ks
    return ktemp    
# compute the sequence of couplings on the Fib chain
kf = [k(i+1) for i in range(L+1)]
# harmonic mean stiffness
meank = hmean(kf) #should be equal to (p+q)/(p/ks+q/kw)

"""
if ks >= kw, then we have for the eigenmodes:
omega_max <= omega^p_max,
where omega^p_max is the max eigenmode of the periodic chain where ks = kw.
We have omega^p_max = 2 \sqrt{ts}.
Since the leapfrog scheme is stable if dt < 2/omega_max, it is enougth to take
dt = 1/\sqrt{ts}
"""
# timestep
dt = 0.9/math.sqrt(ks)

# force from harmonic potential
def HARMONIC(u,i):
    return meank*(u[i+1]+u[i-1]-2.*u[i])

# force from harmonic potential with Fibonacci stiffness
def FIBO(u,i):
    return k(i+1)*u[i+1]+k(i)*u[i-1]-(k(i)+k(i+1))*u[i]

# Fibonacci potential between sites n and n-1
def pot_fib(u,i):
    return 0.5*k(i)*(u[i]-u[i-1])**2
    
# compute the total energy in the harmonic chain
def en_harmo(u,v):
    return sum([0.5*v[0]**2+0.5*meank*u[0]**2]+[0.5*v[n]**2 + 0.5*meank*(u[n]-u[n-1])**2 for n in range(1,L)]+[0.5*v[L]**2+0.5*meank*(u[L]-u[L-1])**2+0.5*meank*u[L]**2])

def en_fibo(u,v):
    return sum([0]+[0.5*v[n]**2 + pot_fib(u,n) for n in range(1,L)]+[pot_fib(u,L)])

# Fourier transform of a function e, on the chain of size L
def FT(e,k):
    return sum([e[n]*math.sin(n*k*math.pi/L) for n in range(1,L)])/math.sqrt(L/2.)
    
#initial condition for the position u and velocity v
u = [0]+[math.sin(i*math.pi/L)for i in range(1,L)]+[0]
v = [0. for i in range(L+1)]
# store energy as a function of time
enf = [en_fibo(u,v)]
# initial velocity for the control system
vc = v[:]
# we advance velocities from 1/2 time step to initiate the leapfrog integration
v = [0]+[v[i]+0.5*dt*FIBO(u,i) for i in range(1,L)]+[0]

# looking at the first Nmodes Fourier energy modes
Nmodes = 5
modc = [[] for i in range(Nmodes)]
mod = [[] for i in range(Nmodes)]

# vector of displacement and velocity
# /!\ velocity is avdanced by 1/2 time step
ut = [u]
vt = [v]
# control system (harmonic potential)
uc = u[:]
vc = [0]+[vc[i]+0.5*dt*HARMONIC(u,i) for i in range(1,L)]+[0]
utc = [u]
vtc = [v]
# store energy as a function of time
enc = [en_harmo(uc,vc)]

#dynamics with fixed boundary conditions: u[0] = u[L-1] = 0. Using Euler integrator
for itime in range(Ntime):
    olduc = uc[:]
    oldvc = vc[:]
    # compute u at time t+1 for u at time t
    uc = [0]+[olduc[i]+dt*oldvc[i] for i in range(1,L)]+[0]
    # compute v at time t+1/2 from v at time t-1/2
    vc = [0]+[oldvc[i]+dt*HARMONIC(uc,i) for i in range(1,L)]+[0]
    # compute  v at time t+1 from v at time t+1/2
    vsynchro = [0]+[vc[i]+0.5*dt*HARMONIC(uc,i) for i in range(1,L)]+[0]
        
    utc.append(uc)
    vtc.append(vsynchro)
    enc.append(en_harmo(uc,vsynchro))
    
    # first Nmodes Fourier modes, harmonic
#    for k in range(Nmodes):
#        enc = FT(vc,k)**2/2.+2.*meank*(math.sin(k*math.pi/(2*L)))**2*FT(uc,k)**2
#        modc[k].append(enc)
    
    oldu = u[:]
    oldv = v[:]
    u = [0]+[oldu[i]+dt*oldv[i] for i in range(1,L)]+[0]
    v = [0]+[oldv[i]+dt*FIBO(u,i) for i in range(1,L)]+[0]
    # compute  v at time t+1 from v at time t+1/2
    vsynchro = [0]+[vc[i]+0.5*dt*HARMONIC(uc,i) for i in range(1,L)]+[0]
    
    ut.append(u)
    vt.append(vsynchro)
    enf.append(en_fibo(u,vsynchro))
    
    # first Nmodes Fourier modes, Fibonacci    
#    for k in range(Nmodes):
#        en = FT(v,k)**2/2.+2.*meankg*(math.sin(k*math.pi/(2*L)))**2*FT(u,k)**2
#        mod[k].append(en)
    

def en():
    plt.title('Power spectrum of the Fibonacci and harmonic chains')
    plt.xlabel('time')
    plt.ylabel('energy of modes ' + str(0) + ' to ' + str(Nmodes-1))
    time = [dt*itime for itime in range(Ntime)]
    #exact = [L*(L*math.sin(math.pi/(2*L)))**2 for itime in range(Ntime)]
    #plt.plot(time,exact,'b-')
    for k in range(Nmodes): plt.plot(time,modc[k],'r-')
    for k in range(Nmodes): plt.plot(time,mod[k],'b-')
    modsum = [sum(sublist) for sublist in zip(mod[0],mod[1],mod[2],mod[3],mod[4])]
    plt.plot(time,modsum,'g-')
    plt.axis([0,Ntime*dt,min(modc[1])-0.2,max(modc[1])+0.2])#,75.,82.])
    plt.show()

r = range(L+1)
def PLOT(timeit):
    plt.title('The Fibonacci chain at time ' + str(round(dt*timeit,2)))
    plt.xlabel('position')
    plt.ylabel('displacement')
    plt.axis([0,L+1,-2.,2.])
    plt.plot(r,utc[timeit])
    plt.plot(r,vtc[timeit])
    plt.show()
    
def plot_c(timeit):
    plt.title('The Fibonacci and harmonic chains at time ' + str(round(dt*timeit,2)))
    plt.xlabel('position')
    plt.ylabel('displacement')
    plt.axis([0,L+1,-2.,2.])
    plt.plot(r,ut[timeit])
    plt.plot(r,utc[timeit])
    plt.show()

def record_FIBO(timeit):
    plt.title('The Fibonacci chain at time ' + str(round(dt*timeit,2)))
    plt.xlabel('position')
    plt.ylabel('displacement')
    plt.axis([0,L+1,-2.,2.])
    plt.plot(r,ut[timeit])
    plt.savefig('data/fibo_realspace_'+str(timeit)+'.png')
    # clear plot
    plt.clf()
    
def record_control(timeit,skip):
    plt.title('The Fibonacci and harmonic chains at time ' + str(round(skip*dt*timeit,2)))
    plt.xlabel('position')
    plt.ylabel('displacement')
    plt.axis([0,L+1,-2.,2.])
    plt.plot(r,ut[skip*timeit])
    plt.plot(r,utc[skip*timeit])
    plt.savefig('data/fibo_harmo_'+str(timeit)+'.png')
    # clear plot
    plt.clf()

# For video output: mencoder "mf://FPU_realspace_%d.png" -mf fps=12:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell:vbitrate=7000 -vf scale=800:600 -o FPU_movie.avi


#r=range(L+1)
#for t in range(100,400,100):
#    plt.title('The Fibonacci chain at several times')
#    plt.xlabel('position')
#    plt.ylabel('displacement')
#    plt.plot(r,ut[int(t/dt)])

#time = [dt*itime for itime in range(Ntime)]
#plt.title('The FPU chain')
#plt.xlabel('time')
#plt.ylabel('energy of modes 1 to 4')
#plt.plot(time,modes[1],'r-')
#plt.plot(time,modes[2],'b-')
#plt.plot(time,modes[3],'b-')
#plt.plot(time,modes[4],'b-')
#plt.show()