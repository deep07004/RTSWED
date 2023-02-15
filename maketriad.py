#!/usr/bin/env python

from scipy.spatial import Delaunay
import numpy as np
import pandas as pd
from geographiclib.geodesic import Geodesic
from triad import Triad, Triads

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
      traids : list of triad
    """
    gl = Geodesic.WGS84
    sta_list = [[nt.code+'.'+sta.code,sta.latitude,sta.longitude,sta.elevation] for nt in inv for sta in nt.stations]

    stations = pd.DataFrame(sta_list)
    #stations.set_index(0,inplace=True)
    pp = np.array(stations[[2,1]])
    tri = Delaunay(pp)
    ss=(np.zeros(np.shape(tri.simplices)[0])==0)
    td = Triads()
    for l, sipmlices in enumerate(tri.simplices):
        sta = []
        for i in sipmlices:
            a = [stations.iloc[i][0],stations.iloc[i][1],stations.iloc[i][2]]
            sta.append(a)
        _tmp_td = Triad(sta)
        # Check the side length
        # Check for interior angle
        dist = _tmp_td.dist
        ang = _tmp_td.ang
        dist.sort()
        ang.sort()
        if dist[0] < minlen or dist[2] > maxlen :
            continue
        if ang[0] < minang or ang[2] > maxang:
            continue
        td.append(_tmp_td)
    return td