# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 22:49:47 2021

@author: 17000
"""

import os;
import matplotlib.pyplot as plt;
import numpy as np;
import json;

route = os.getcwd();
flist = os.listdir(route+'\\experimental_data');
rg1 = [0.4,1.4];
rg2 = [1.5,2.5];
rg3 = [2.6,3.6];

rg1 = [2000,7000];
rg2 = [7500,12500];
rg3 = [13000,18000];

def level(rg,cret):
    res = mat[2][rg[0]:rg[1]];
    sigma = np.std(res);
    k1,k2 = 0,0;
    k1n = np.average(res)+sigma;
    k2n = k1n-2*sigma;
    while(abs(k1n-k1)>cret or abs(k2n-k2)>cret):
        ave = (k1n+k2n)/2;
        k1l = [];
        k2l = [];
        for x in res:
            if(x<ave):
                k1l.append(x);
            else:
                k2l.append(x);
        k1,k2 = k1n,k2n;
        k1n = np.average(k1l);
        k2n = np.average(k2l);
    return {'ave':ave,'dI':abs(k1n-k2n)};

def tau(rg,ave):
    res = mat[2][rg[0]:rg[1]];
    bl = [int(u1>ave) for u1 in res]
    jmp = [bl[i+1]-bl[i] for i in range(len(bl)-1)];
    tau0 = [];
    tau1 = [];
    init=0;
    current = 0;
    while(current<len(jmp)):
        val = jmp[current];
        if(val==1):
            if(init!=0):
                tau0.append(current-init);
            init=current;
        elif(val==-1):
            if(init!=0):
                tau1.append(current-init);
            init=current;
        current += 1;
    return {'tau0':np.array(tau0)*2*10**-4,'tau1':np.array(tau1)*2*10**-4};

t0 = [[] for i in range(3)];
t1 = [[] for i in range(3)];
dI = [[] for i in range(3)];
rglist = [rg1,rg2,rg3];

for file_name in flist:
    file = 'experimental_data\\'+file_name;
    if(file[-4:]=='.csv'):
        with open(file,'r') as file:
            data = [[float(v) for v in u[:-1].split(',')] for u in file.readlines()[1:]];
        mat = np.transpose(data);
        for i1 in range(3):
            res = level(rglist[i1],10**-8);
            ave = res['ave'];
            tres = tau(rglist[i1],ave);
            t0[i1] += list(tres['tau0']);
            t1[i1] += list(tres['tau1']);
            dI[i1].append(res['dI'])

def smear(tlist, sigma, t):
    return np.average([np.exp(-(t-t0)**2/(2*sigma**2))/(sigma*np.sqrt(2*np.pi)) for t0 in tlist]);

out = {};
for i1 in range(3):
    tlist = [i*0.0008 for i in range(50)];
    res = [smear(t0[i1],0.002,t) for t in tlist];
    plt.plot(tlist,res);
    plt.savefig('t0_'+str(i1)+'.png');
    plt.close();
    out['t0_'+str(i1)]=res[:];
    res = [smear(t1[i1],0.001,t) for t in tlist];
    plt.plot(tlist,res);
    plt.savefig('t1_'+str(i1)+'.png');
    plt.close();
    out['t1_'+str(i1)]=res[:];
out['tlist'] = tlist;
with open('data.json','w') as file:
    json.dump(out,file);
    

import json;
import numpy as np;
from itertools import product;

with open('data.json','r') as file:
    data = json.load(file);

file = open('fig_output.txt','w');
file.write('time (s)\t tau_0 (-0.15 V, exp) \t tau_0 (-0.20 V, exp)\t tau_0 (-0.25 V, exp)\t');
file.write('tau_1 (-0.15 V, exp)\t tau_1 (-0.20 V, exp)\t tau_1 (-0.25 V, exp)\t');
file.write('tau_0 (simu)\t')
file.write('tau_1 (-0.15 V, simu)\t tau_1 (-0.20 V, simu)\t tau_1 (-0.25 V, simu)\n');

for i in range(1,len(data['tlist'])):
    file.write(str(data['tlist'][i])+'\t');
    for u,v in product(range(2),range(3)):
        file.write(str(data['t'+str(u)+'_'+str(v)][i]));
        file.write('\t');
            
    A = data['tlist'][i];
    def fitting(C,t,td):
        return C/(t-td)*np.exp(-A/t)*(1-np.exp(-(t-td)/t/td*A));

    fitting1 = fitting(0.92719,0.01152,0.001);    
    fitting2 = fitting(0.95,0.01316,0.0005);
    fitting3 = fitting(0.95,0.00754,0.0005);
    fitting4 = fitting(0.95,0.00435,0.0005);
    file.write(str(fitting1)+'\t'+str(fitting2)+'\t'+str(fitting3)+'\t'+str(fitting4)+'\n');
file.close();
