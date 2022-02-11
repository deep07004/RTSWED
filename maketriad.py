#!/usr/bin/env python

# This program form triad subarrays from the input station list
# using delaunay trianglation
# Inputs:
#   station : array with station information ['ID', lat, lon, elevation]
#   minlen : minimum side length of a triad in km
#   maxlen : maximum side length of a triad in km
#   minang : minimum interior angle of a triad
#   maxang : maximum interior angle of a triad
# Outoput:
#   traids :  array ['ID', lat, lon, elevation, triad index]
#   ntriad : no of triads
from scipy.spatial import Delaunay
import numpy as np

def maketriads(stations, minlen=10, maxlen=600, minang=30, maxang=120):
    nsta = len(stations)
    pp = np.zeros((nsta,2))
    for i in range(nsta):
        pp[i,0] = stations[i,1]
        pp[i,0] = stations[i,2]
    return triads, ntriad
