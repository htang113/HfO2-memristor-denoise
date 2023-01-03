# -*- coding: utf-8 -*-
"""
Created on Mon May 16 15:03:28 2022

@author: 17000
"""

import numpy as np;
from scipy import special;
import matplotlib.pyplot as plt;
from itertools import product;
import json;
import os;

recalc_pot = True;

a = 3;  #nm
#d = 1;  #nm
const = 0.03810;  #nm^2 eV
m = 0.11; #effective mass, in me
coef = const/m;
L = 6; #nm

epsilon = 16;
T = 300;
n = 5*10**18; #doping concentration, in cm^-3
lD = np.sqrt(epsilon*T/n)*0.0690*10**9;
E = 0.2; #eV; 
k0 = np.sqrt(E/coef)

nx = 2;
ny = 2;
nz = 100;
chan = [];
k1 = 1.440;
d1,d2 = 1, 0.8;
for kl in [(0,0),(0,1),(1,0),(1,1)]: #eV*nm;
#k1=0

    mnl = [[n+1,m+1] for m,n in product(range(nx),range(ny))];
    kmn = [np.sqrt(k0**2-(mn[0]*np.pi/a)**2-(mn[1]*np.pi/a)**2+2*(np.pi/a)**2+0j) for mn in mnl];
    zl = np.linspace(-L/2,L/2,nz+1);
    dz = zl[1]-zl[0];
    
    def V(x,y,z):
        
        r1 = np.sqrt((x+d1)**2+(y-a/2-0.5)**2+z**2);
        r2 = np.sqrt((x+d2)**2+(y-a/2+0.5)**2+z**2);
        
        return kl[0]*k1*np.exp(-r1/lD)/r1+kl[1]*k1*np.exp(-r2/lD)/r1;
    
    Nx = 20;
    xl = np.linspace(a/Nx/2,a-a/Nx/2,Nx);
    
    def Vk(mn,mnp,z):
        
        ls = [[V(x,y,z)*(2/a)**2*np.sin(mn[0]*np.pi*x/a)*np.sin(mn[1]*np.pi*y/a)*np.sin(mnp[0]*np.pi*x/a)*np.sin(mnp[1]*np.pi*y/a) for y in xl] for x in xl]
        
        return np.sum(ls)*a**2/Nx**2;
    
    if(('V.json' in os.listdir()) and recalc_pot==False):
        with open('V.json','r') as file:
            Vl = json.load(file);
    
    else:
        Vl = [[[Vk(mn,mnp,z) for mnp in mnl] for mn in mnl] for z in zl]
        with open('V.json','w') as file:
            json.dump(Vl,file);
            
    A = np.zeros([nx*ny*(nz+3),nx*ny*(nz+3)])*1j;
    b = np.zeros(nx*ny*(nz+3));
    
    for i in range(nx*ny):
        A[i,i],A[i,nx*ny+i]=-1,1;
        A[nx*ny+i,i] = kmn[i];
        A[nx*ny+i,nx*ny+i],A[nx*ny+i,2*nx*ny+i] = 1j/dz, -1j/dz;
        b[0] = 1;
        b[nx*ny] = k0;
    
    for zi in range(1,nz):
        for i,j in product(range(nx*ny),range(nx*ny)):
            A[(zi+1)*nx*ny+i][(zi+1)*nx*ny+j] += Vl[zi][i][j];
            if(i==j):
                A[(zi+1)*nx*ny+i][(zi+1)*nx*ny+i] += -E+coef*((mnl[i][0]**2+mnl[i][1]**2-2)*(np.pi/a)**2+2/dz**2);
                A[(zi+1)*nx*ny+i][(zi)*nx*ny+i] += -coef/dz**2;
                A[(zi+1)*nx*ny+i][(zi+2)*nx*ny+i] += -coef/dz**2;
                
    for i in range(nx*ny):
        A[(nz+1)*nx*ny+i,(nz+1)*nx*ny+i],A[(nz+1)*nx*ny+i,(nz+2)*nx*ny+i]=-1,1;
        A[(nz+2)*nx*ny+i,(nz+2)*nx*ny+i] = kmn[i];
        A[(nz+2)*nx*ny+i,(nz+1)*nx*ny+i],A[(nz+2)*nx*ny+i,(nz)*nx*ny+i] = 1j/dz, -1j/dz;
    
    x = np.linalg.solve(A,b);
    channel = [np.linalg.norm(x[i*nx*ny:(i+1)*nx*ny])**2 for i in range(1,nz+2)];
    thick = int(nz/5)
    source = [np.abs(np.exp(-2j*k0*dz*(thick-t))+x[0])**2+np.linalg.norm([x[i]*np.exp(1j*(thick-t)*dz*kmn[i]) for i in range(1,nx*ny)])**2 for t in range(thick)];
    drain = [np.linalg.norm([x[(nz+2)*nx*ny+i]*np.exp(1j*t*dz*kmn[i]) for i in range(nx*ny)])**2 for t in range(thick)];
    plt.plot(channel);
    chan.append(channel);
    plt.savefig('wave_function.pdf');
    
with open('wave_function.txt','w') as file:
    
    file.write('z (nm)\t No trap\t Trap (0.8 nm)\t Trap (1 nm)\t Trap (both)\n')
    
    for i in range(len(chan[0])):
        file.write(str(zl[i])+'\t')
        for ch in chan:
            file.write(str(ch[i])+'\t');
        file.write('\n');
