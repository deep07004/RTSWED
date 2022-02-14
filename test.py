#!/usr/bin/env python
import sys
from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt
from geographiclib.geodesic import Geodesic
from obspy import read_inventory

gl = Geodesic.WGS84

sta_list=[]
inv = read_inventory(sys.argv[1])
for nt in inv.networks:
    for sta in nt.stations:
        sta_list.append([nt.code+'.'+sta.code,sta.latitude,sta.longitude,sta.elevation])
n=len(sta_list)
pp = np.zeros((n,2))
for i in range(n):
    ff = sta_list[i]
    pp[i,0] = np.float64(ff[2])
    pp[i,1] = np.float64(ff[1])
tri = Delaunay(pp)
ss=(np.zeros(np.shape(tri.simplices)[0])==0)
### for l, xx in enumerate(tri.simplices):
###     for i in range(3):
###         j = (i+1)%3
###         a = xx[i]
###         b = xx[j]
###     # Check the side length
###         side_length = gl.Inverse(pp[a,1],pp[a,0],pp[b,1],pp[b,0])
###         if side_length['s12']/1000 < 10 or side_length['s12']/1000 > 600:
###             ss[l] = False
###     # Check for interior angle
for l, xx in enumerate(tri.simplices):
    dist = np.zeros(3)
    ang = np.zeros(3)
    dist[0] = gl.Inverse(pp[xx[0],1],pp[xx[0],0],pp[xx[1],1],pp[xx[1],0])['s12']/1000.0
    dist[1] = gl.Inverse(pp[xx[0],1],pp[xx[0],0],pp[xx[2],1],pp[xx[2],0])['s12']/1000.0
    dist[2] = gl.Inverse(pp[xx[1],1],pp[xx[1],0],pp[xx[2],1],pp[xx[2],0])['s12']/1000.0
    dist.sort()
    ang[0] = np.arccos((dist[1]**2 + dist[2] **2 - dist[0] **2) / (2 * dist[1] * dist[2]))
    ang[1] = np.arccos((dist[0]**2 + dist[2] **2 - dist[1] **2) / (2 * dist[0] * dist[2]))
    ang[2] = np.arccos((dist[0]**2 + dist[1] **2 - dist[2] **2) / (2 * dist[0] * dist[1]))
    ang *= 180.0/np.pi
    ang.sort()
    if dist[0] < 50 or dist[2] > 600 :
        ss[l] = False
    if ang[0] < 30 or ang[2] > 120:
        ss[l] = False

plt.triplot(pp[:,0],pp[:,1],tri.simplices[ss])
plt.show()