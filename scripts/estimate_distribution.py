#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import sys 

filename = sys.argv[1]
kmap = {}

with open(filename, 'r') as f:
    for line in f:
        x, y = line.strip().split('\t')
        kmap[int(x)] = float(y)

arrx = np.arange(1,max(kmap.keys()) + 1)
arry = [kmap[x] if x in kmap else 0.0 for x in np.arange(1,max(kmap.keys()) + 1)]
df = pd.Series(arry, index=arrx)
rolling5 = df.rolling(window=5).mean()
rolling20 = df.rolling(window=20).mean()
rolling50 = df.rolling(window=50).mean()
ewm5 = df.ewm(span=5,adjust=False).mean()
ewm20 = df.ewm(span=20,adjust=False).mean()
ewm50 = df.ewm(span=50,adjust=False).mean()

def get_eq_pos(data):
    try:
        return [len(set(dlist))==1 and dlist[0] > 0.0 for dlist in zip(data[0:],data[1:],data[2:],data[3:],data[4:])].index(True)
    except ValueError:  
        return 'NA'

def get_smaller_pos(data1, data2, length):
    smaller = data1 < data2
    tozip = [smaller[x:] for x in np.arange(0,length)]
    try:
        return [len(set(dlist))==1 and dlist[0] > 0.0 for dlist in zip(*tozip)].index(True)
    except ValueError:  
        return 1

ewm5min = get_smaller_pos(ewm5.cummin(), ewm5, 6)
rolling5min = get_smaller_pos(rolling5.cummin(), rolling5, 6)
dfmin = get_smaller_pos(df.cummin(), df, 6)
first = get_smaller_pos(df.cummin(), df, 1)
inflexion = np.diff(df, 2).tolist()
inf_point = inflexion.index(list(filter(lambda x: x<0, inflexion))[0])

#plt.semilogy(df, color='blue', label='original data')
#plt.semilogy((df.cummin()), '--', color='blue', label='cum min original data')
#plt.semilogy(rolling5, color='red', label='rolling window (length 5)')
#plt.semilogy((rolling5.cummin()), '--', color='red', label='cum min rolling')
#plt.semilogy(ewm5, color='green', label='exponential window (length5)')
#plt.semilogy((ewm5.cummin()), '--', color='green', label='cum min exponential')
#plt.axis([0,200,0.1,max(kmap.values())])
#plt.xlabel("Depth of kmer")
#plt.ylabel("#kmers with depth x")
#plt.title("Smoothing the kmer histogram")

#plt.loglog(rolling5, color='red', label='rolling window (length 5)')
#plt.loglog((rolling5.cummin()), '--', color='red', label='cum min rolling')
plt.loglog(df, color='blue', label='original data')
plt.loglog((df.cummin()), '--', color='blue', label='cum min original data')
#plt.loglog(ewm5, color='green', label='exponential window (length 5)')
#plt.loglog((ewm5.cummin()), '--', color='green', label='cum min exponential')
plt.axis([1,max(kmap.keys()),0.1,max(kmap.values())])
plt.xlabel("Depth of kmer")
plt.ylabel("#kmers with depth x")
plt.title("k-mer histogram of cluster")

#plt.plot(df, color='blue')
#plt.plot((df.cummin()), '--', color='blue')
#plt.plot(rolling5, color='red')
#plt.plot((rolling5.cummin()), '--', color='red')
#plt.plot(ewm5, color='green')
#plt.plot((ewm5.cummin()), '--', color='green')
#plt.axis([0,max(ewm5min, rollinblueg5min, dfmin) + 50,0,5])

plt.axvline(dfmin - 1, color='red')
#plt.axvline(first, color='black')
#plt.axvline(inf_point - 1, color='green')
#plt.axvline(turn_point, color='red')
#plt.axvline(ewm5min, color='green')
#plt.axvline(rolling5min, color='red')
#plt.axvline(2, color='black')
plt.legend(loc='best', prop={'size' : 8})
plt.savefig(sys.argv[2], dpi=300)

#print get_eq_pos(ewm5.cummin())
#print get_eq_pos(rolling5.cummin())
#print get_eq_pos(df.cummin())

#a = np.array([kmap[x] if x in kmap else 0.0 for x in np.arange(1,max(kmap.keys()))])
