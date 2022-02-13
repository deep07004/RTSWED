#!/usr/bin/env python
from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt

fp=open('stations.txt').readlines()
n=len(fp)
pp = np.zeros((n,2))
for i in range(n):
    ff = fp[i].split(' ')
    pp[i,0] = np.float64(ff[2])
    pp[i,1] = np.float64(ff[1])
tri = Delaunay(pp)
for xx in tri.simplices:
    for i in range(xx):
        j=(i+1)%3
        a = pp[i,0]
        b = pp[j,0]


plt.triplot(pp[:,0],pp[:,1],tri.simplices.copy())
a = tri.simplices[0]
plt.plot(pp[a,0],pp[a,1],'r*')
plt.show()

###  thresh = 1.0  # user defined threshold
###  small_edges = set()
###  large_edges = set()
###  for tr in tri.vertices:
###      for i in xrange(3):
###          edge_idx0 = tr[i]
###          edge_idx1 = tr[(i+1)%3]
###          if (edge_idx1, edge_idx0) in small_edges:
###              continue  # already visited this edge from other side
###          if (edge_idx1, edge_idx0) in large_edges:
###              continue
###          p0 = points[edge_idx0]
###          p1 = points[edge_idx1]
###          if np.linalg.norm(p1 - p0) <  thresh:
###              small_edges.add((edge_idx0, edge_idx1))
###          else:
###              large_edges.add((edge_idx0, edge_idx1))
###  
###  # Plotting the output
###  figure()
###  plot(points[:, 0], points[:, 1], '.')
###  for i, j in small_edges:
###      plot(points[[i, j], 0], points[[i, j], 1], 'b')
###  for i, j in large_edges:
###      plot(points[[i, j], 0], points[[i, j], 1], 'c')
###  