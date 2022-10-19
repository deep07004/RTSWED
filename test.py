#!/usr/bin/env python
from obspy import read_inventory
import pandas as pd
import pygmt
import maketriad

inv = read_inventory('Inventory.xml')
triads, stations = maketriad.maketriad(inv, minlen=10, maxlen=1200, minang=10, maxang=120)
reg = [ 105, 182, -48, -2 ]
grid = pygmt.datasets.load_earth_relief(resolution="02m", region=reg)
fig = pygmt.Figure()
pygmt.makecpt(cmap="terra", series=[-8000, 8000])
fig.basemap(region=reg, projection="M15c",frame=True)
fig.grdimage(grid=grid,shading=True)
for i in triads.index:
    lat = []
    lon = []
    for a in triads.loc[i]:
        _tmp = stations.loc[a]
        lat.append(_tmp[1])
        lon.append(_tmp[2])
    lat.append(lat[0])
    lon.append(lon[0])
#    for j in range(3):
#        fig.plot(x=lon[j:j+2], y=lat[j:j+2], pen='0.35p')
fig.plot(stations[[2,1]],style="t0.25c", color="red")
fig.savefig('triad.jpg')
fig.show()
print(triads)
