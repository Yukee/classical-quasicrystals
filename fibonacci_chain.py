# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 16:38:59 2014

@author: nicolas

At short times the Fibonacci spring chain behaves as a periodic spin chain,
whose spings have a stiffness equal to the harmonic mean of the Fibonacci sprins stiffness.
TODO: behaviour at large times? Is there energy transfer to the higher modes in the Fibonacci case?
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
Ntime = 6000
n = 11
# length of the spring chain
L = fib(n)-1
# stength of the anharmonicity
alpha = 0.8
# timestep
dt = 0.2
# strenght of the quasiperiodic addition to the stiffness
delta = 5.
# mean stiffness
#meank = 1 + delta*fib(n-1)/fib(n)
meank = (1.+delta)/(1.+delta*fib(n-2)/fib(n))
# force from harmonic potential
def HARMONIC(u,n):
    return meank*(u[n+1]+u[n-1]-2.*u[n])

# force from anharmonic (Fermi Pasta Ulam) potential 
def FPU(u,n):
    return u[n+1]+u[n-1]-2.*u[n]+alpha*((u[n+1]-u[n])**2-(u[n]-u[n-1])**2)

# inverse golden ratio
om = 2./(1.+math.sqrt(5))

# force from harmonic potential with Fibonacci stiffness
def FIBO(u,n):
    fibn = int(om*(n+1)) - int(om*n) # an integer to takes values 0 or 1 occordingly to Fibonacci sequence
    return (1+delta*fibn)*(u[n+1]+u[n-1]-2.*u[n])

# Fourier transform of a function e, on the chain of size L
def FT(e,k):
    return sum([e[n]*math.sin(n*k*math.pi/L) for n in range(1,L)])/math.sqrt(L/2.)
    
#initial condition for the position u and velocity v
u = [0]+[math.sin(i*math.pi/L)+0.5*math.sin(2*i*math.pi/L)+0.25*math.sin(3*i*math.pi/L) for i in range(1,L)]+[0]
#u=[0]+[1.-2./L*abs(L/2.-i) for i in range(1,L)]+[0]
# x0 = L/4 +4
# x1 = L+1 - x0
# v0 = 0.02
v = [0. for i in range(L+1)]

# looking at the first Nmodes fourier energy modes
Nmodes = 5
modes = [[] for i in range(Nmodes)]

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
    vc = [0]+[oldvc[i]+dt*HARMONIC(olduc,i) for i in range(1,L)]+[0]
    uc = [0]+[olduc[i]+dt*vc[i] for i in range(1,L)]+[0]
        
    utc.append(uc)
    vtc.append(vc)
    
    oldu = u[:]
    oldv = v[:]
    v = [0]+[oldv[i]+dt*FIBO(oldu,i) for i in range(1,L)]+[0]
    u = [0]+[oldu[i]+dt*v[i] for i in range(1,L)]+[0]
    
    ut.append(u)
    vt.append(v)
#    
#    for k in range(Nmodes):
#        en = FT(v,k)**2/2.+2.*(math.sin(k*math.pi/(2*L)))**2*FT(u,k)**2
#        modes[k].append(en)
    

def en(k):
    plt.title('The Fibonacci chain chain')
    plt.xlabel('time')
    plt.ylabel('energy of mode ' + str(k))
    time = [dt*itime for itime in range(Ntime)]
    #exact = [L*(L*math.sin(math.pi/(2*L)))**2 for itime in range(Ntime)]
    #plt.plot(time,exact,'b-')
    plt.plot(time,modes[k],'r-')
    plt.axis([0,Ntime*dt,min(modes[k])-0.2,max(modes[k])+0.2])#,75.,82.])
    plt.show()

r = range(L+1)
def PLOT(timeit):
    plt.title('The Fibonacci chain at time ' + str(round(dt*timeit,2)))
    plt.xlabel('position')
    plt.ylabel('displacement')
    plt.axis([0,L+1,-2.,2.])
    plt.plot(r,utc[timeit])
    plt.show()
    
def plot_c(timeit):
    plt.title('The Fibonacci and harmonic chains at time ' + str(round(dt*timeit,2)))
    plt.xlabel('position')
    plt.ylabel('displacement')
    plt.axis([0,L+1,-2.,2.])
    plt.plot(r,ut[timeit])
    plt.plot(r,utc[timeit])
    plt.show()

def record(timeit):
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

def record_FPU(timeit):
    plt.title('The FPU chain at time ' + str(round(dt*timeit,2)))
    plt.xlabel('position')
    plt.ylabel('displacement')
    plt.axis([0,L+1,-2.,2.])
    plt.plot(r,ut[timeit])
    plt.savefig('data/FPU_realspace_'+str(int(timeit/3))+'.png')
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