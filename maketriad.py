#!/usr/bin/env python

from scipy.spatial import Delaunay
import numpy as np
import pandas as pd
from geographiclib.geodesic import Geodesic

def maketriad(inv, minlen=10, maxlen=600, minang=30, maxang=120):
    """
    This program form triad subarrays from the input station inventory (fdsnxml)
    using delaunay trianglation
    Inputs:
      inv    : station inventory in fdsnxml format
      minlen : minimum side length of a triad in km
      maxlen : maximum side length of a triad in km
      minang : minimum interior angle of a triad
      maxang : maximum interior angle of a triad
    Outoput:
      traids : panda data frame of triads containing nt.code of
               three stations in each row.
      stations : panda data frame of stations indexed by nt.sta
                column 1 is latitude, column 2 is longitude and
                column 3 is elevation 
    """
    gl = Geodesic.WGS84
    sta_list = []
    for nt in inv.networks:
        for sta in nt.stations:
            sta_list.append([nt.code+'.'+sta.code,sta.latitude,sta.longitude,sta.elevation])

    stations = pd.DataFrame(sta_list)
    stations.set_index(0,inplace=True)
    pp = np.array(stations[[2,1]])
    tri = Delaunay(pp)
    ss=(np.zeros(np.shape(tri.simplices)[0])==0)
    # Check the side length
    # Check for interior angle
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
        if dist[0] < minlen or dist[2] > maxlen :
            ss[l] = False
        if ang[0] < minang or ang[2] > maxang:
            ss[l] = False
    traids_list = []
    for t in tri.simplices[ss]:
        _tmp =[]
        for i in t:
            _a  = stations[stations[1] == pp[i,1]]
            _tmp.append(_a[_a[2] == pp[i,0]].index[0])
        traids_list.append(_tmp)
    triads = pd.DataFrame(traids_list)
    return triads, stations