#!/usr/bin/env python
import sys, os, argparse, yaml, logging, time
from obspy.core import UTCDateTime
from obspy import read_inventory, Inventory
from obspy import read_inventory
import pandas as pd
import pygmt
import maketriad

# To read yaml configuration file to stream
def read(filepath):
    with open(filepath, 'r') as _f:
        return yaml.load(_f, Loader=yaml.FullLoader)

# Parse the commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config', metavar=('FILEPATH'), help='configuration file',\
    default='./config.yaml', type=str)

args = parser.parse_args()
# Read the configuration file and check that all informations are there 
# to proceed further 
cfg = read(args.config)
inv = read_inventory('Inventory.xml')
sta_list = [[nt.code+'.'+sta.code,sta.latitude,sta.longitude,sta.elevation] for nt in inv for sta in nt.stations]
stations = pd.DataFrame(sta_list)
triads = maketriad.maketriad(inv, minlen=10, maxlen=1200, minang=10, maxang=120)
reg = [ 105, 182, -48, -2 ]
grid = pygmt.datasets.load_earth_relief(resolution="02m", region=reg)
fig = pygmt.Figure()
pygmt.makecpt(cmap="terra", series=[-8000, 8000])
fig.basemap(region=reg, projection="M15c",frame=True)
fig.grdimage(grid=grid,shading=True)
for td in triads.triads:
    lat = []
    lon = []
    for a in td.stations:
        lat.append(a[1])
        lon.append(a[2])
    lat.append(lat[0])
    lon.append(lon[0])
    for j in range(3):
        fig.plot(x=lon[j:j+2], y=lat[j:j+2], pen='0.35p')
fig.plot(stations[[2,1]],style="t0.25c", color="red")
fig.savefig('triad.jpg')
fig.show()
print(triads)
