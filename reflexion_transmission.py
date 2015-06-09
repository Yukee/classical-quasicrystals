# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 13:31:11 2014

@author: nicolas

At short times the Fibonacci spring chain behaves as a periodic spin chain,
whose spings have a stiffness equal to the harmonic mean of the Fibonacci sprins stiffness.

"""

import math
import matplotlib.pyplot as plt

#compute Fibonacci numbers
def fib(n):
    a, b = 0, 1
    for i in range(n):
        a, b = b, a + b
    return a

# total time of the simulation
Ntime = 9000
n = 11
# length of the Fibonacci spring chain
L = fib(n)-1
# length of the periodic sping chain
Lper = 120
# timestep
dt = 0.02
# strenght of the quasiperiodic addition to the stiffness
delta = 5.
# harmonic mean stiffness
meank = (1.+delta)/(1.+delta*fib(n-2)/fib(n))
# geometric mean stiffness
meankg = 1. + delta*fib(n-1)/fib(n)
# stiffness of the connected periodic chain. kper plays the role of a**2, where a is the spacing btw two consecutive atoms
kper = meank
# amplitude of the incident wave
amp = .2

# force from harmonic potential
def HARMONIC(u,n):
    return meank*(u[n+1]+u[n-1]-2.*u[n])

# force for the connected periodic (harmonic chain)
def PERIODIC(u,n):
    return kper*(u[n+1]+u[n-1]-2.*u[n])
    
# inverse golden ratio
om = 2./(1.+math.sqrt(5))

# force from harmonic potential with Fibonacci stiffness
def FIBO(u,n):
    fibn = int(om*(n+1)) - int(om*n) # an integer taking values 0 or 1 occordingly to Fibonacci sequence
    fibnm = int(om*n) - int(om*(n-1))
    cn = 1.+delta*fibn
    cnm = 1.+delta*fibnm
    return cn*u[n+1]+cnm*u[n-1]-(cn+cnm)*u[n]

# Fourier transform of a function e, on the chain of size L
def FT(e,k):
    return sum([e[n]*math.sin(n*k*math.pi/L) for n in range(1,L)])/math.sqrt(L/2.)
    
#initial condition for the position u and velocity v
puls = .1

u = [0.]+[0. for i in range(1,Lper)]+[0. for i in range(Lper,L+Lper)]+[0.]

v = [0. for i in range(L+Lper+1)]

# looking at the first Nmodes Fourier energy modes
Nmodes = 5
modc = [[] for i in range(Nmodes)]
mod = [[] for i in range(Nmodes)]

# vector of displacement and velocity
ut = [u]
vt = [v]
# control system (harmonic potential)
utc = [u]
vtc = [v]
uc = u[:]
vc = v[:]

#dynamics with fixed boundary conditions: u[0] = u[L-1] = 0. Using Euler integrator
for itime in range(Ntime):
    olduc = uc[:]
    oldvc = vc[:]
    
    if(itime*dt < 2*math.pi/puls):
        vc = [amp*puls*math.cos(puls*itime*dt)]+[oldvc[i]+dt*PERIODIC(olduc,i) for i in range(1,Lper)]+[oldvc[i]+dt*HARMONIC(olduc,i) for i in range(Lper,L+Lper)]+[0.]
        uc = [amp*math.sin(puls*itime*dt)]+[olduc[i]+dt*vc[i] for i in range(1,Lper)]+[olduc[i]+dt*vc[i] for i in range(Lper,L+Lper)]+[0]
    else:
        vc = [0]+[oldvc[i]+dt*PERIODIC(olduc,i) for i in range(1,Lper)]+[oldvc[i]+dt*HARMONIC(olduc,i) for i in range(Lper,L+Lper)]+[0]
        uc = [0]+[olduc[i]+dt*vc[i] for i in range(1,Lper)]+[olduc[i]+dt*vc[i] for i in range(Lper,L+Lper)]+[0]
        
    utc.append(uc)
    vtc.append(vc)
    
    oldu = u[:]
    oldv = v[:]
    
    if(itime*dt < 2*math.pi/puls):
        v = [amp*puls*math.cos(puls*itime*dt)]+[oldv[i]+dt*PERIODIC(oldu,i) for i in range(1,Lper)]+[oldv[i]+dt*FIBO(oldu,i) for i in range(Lper,L+Lper)]+[0.]
        u = [amp*math.sin(puls*itime*dt)]+[oldu[i]+dt*v[i] for i in range(1,Lper)]+[oldu[i]+dt*v[i] for i in range(Lper,L+Lper)]+[0]
    else:
        v = [0]+[oldv[i]+dt*PERIODIC(oldu,i) for i in range(1,Lper)]+[oldv[i]+dt*FIBO(oldu,i) for i in range(Lper,L+Lper)]+[0]
        u = [0]+[oldu[i]+dt*v[i] for i in range(1,Lper)]+[oldu[i]+dt*v[i] for i in range(Lper,L+Lper)]+[0]

    ut.append(u)
    vt.append(v)

    
aper = math.sqrt(kper)
ac = math.sqrt(meank)
r = [i*aper for i in range(Lper)]+[Lper*aper+i*ac for i in range(L+1)]
    
def record_control(timeit,skip):
    plt.title('The Fibonacci and harmonic chains at time ' + str(round(skip*dt*timeit,2)) + '.\n Incident plane wave at $\omega$ = ' +str(puls))
    plt.xlabel('position')
    plt.ylabel('displacement')
    plt.axis([0,ac*(L+1)+aper*Lper,-1.2,1.2])
    plt.plot(r,ut[skip*timeit])
    plt.plot(r,utc[skip*timeit])
    plt.savefig('data/propagation_test'+str(timeit)+'.png',format='png',dpi=120)
    # clear plot
    plt.clf()

# For video output: mencoder "mf://FPU_realspace_%d.png" -mf fps=12:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell:vbitrate=7000 -vf scale=800:600 -o FPU_movie.avi
# or rather mencoder "mf://propagation_test%d.png" -mf fps=12:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell:vbitrate=7000 -vf scale=960:720 -o reflection.avi