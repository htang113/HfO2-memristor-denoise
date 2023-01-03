# -*- coding: utf-8 -*-
"""
Created on Tue May 10 13:35:56 2022

@author: 17000
"""

import numpy as np;
import matplotlib.pyplot as plt;
from matplotlib import cm
from itertools import product;
import os;
import json;

dt = 0.0005;
K = 0.05;
size = (640,480);
steps = 50;

cmap1 = 'jet'

with open('energy_profile.txt','r') as file:
    raw = file.readlines();
data = np.array([[float(v) for v in u.split()] for u in raw]);

with open('initial_state.json','r') as file:
    data1 = json.load(file);
ma = np.max(data1);
mi = np.min(data1);
sp = np.shape(data1);
data1 = data1/(ma-mi)-mi/(ma-mi)*np.ones(sp);
sp = [int(sp[0]/10),int(sp[1]/10)];

init = [[data1[-10*i][10*j] for j in range(sp[1])] for i in range(sp[0])]
evolution = [np.array(init)];

def df(x):
    n1 = int(round(20*x));
    if(n1==20):
        n1=19;
    return (data[n1+1][1]-data[n1][1])*20;
def renorm(y1):
    for u,v in product(range(sp[0]),range(sp[1])):
        if(y1[u,v]>1):
            y1[u,v]=1;
        elif(y1[u,v]<0):
            y1[u,v]=0;
    return y1;

def evol(steps,bias):
    for i in range(steps):
        var = np.zeros([sp[0],sp[1]]);
        for u,v in product(range(sp[0]),range(sp[1])):
            if(u==sp[0]-1 or v==sp[1]-1):
                var[u,v]=0;
            else:
                data0 = evolution[-1];
                var[u,v] = -df(data0[u,v]) + K*(data0[u+1,v]+data0[u,v+1]+data0[u-1,v]+data0[u,v-1]-4*data0[u,v])*sp[0]*sp[1]-bias;
        res = renorm(evolution[-1]+var*dt);
        evolution.append(res);

fig, axs = plt.subplots(1, 4, figsize=(14,5))

axs[0].contourf((ma-mi)*evolution[0]+mi*np.ones(sp),np.linspace(0,1,50),cmap=cmap1);

evol(60,0)
axs[1].contourf((ma-mi)*evolution[-1]+mi*np.ones(sp),np.linspace(0,1,50),cmap=cmap1);

evol(10,30);
evol(10,0);
axs[2].contourf((ma-mi)*evolution[-1]+mi*np.ones(sp),np.linspace(0,1,50),cmap=cmap1);

evol(10,30);
evol(10,0);
axs[3].contourf((ma-mi)*evolution[-1]+mi*np.ones(sp),np.linspace(0,1,50),cmap=cmap1);
plt.savefig('Fig_4.pdf')
