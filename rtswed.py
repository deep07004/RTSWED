#!/usr/bin/env python
import sys, os, argparse, yaml, logging, time
from obspy.core import UTCDateTime
from obspy import read_inventory, Inventory
from obspy import read_inventory
import pandas as pd
import numpy as np
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
inv = read_inventory(os.path.join(cfg['Work_dir'],cfg['Inventory']))
# Cleanup the inventory and select only three channels per station
_tmp_HH = inv.select(channel='HH?')
_tmp_BH = inv.select(channel='BH?')
_tmp_BH_stations = [ sta.code for nt in _tmp_BH for sta in nt]
_tmp_HH_stations = [ sta.code for nt in _tmp_HH for sta in nt]
for a in _tmp_BH_stations:
    if a not in _tmp_HH_stations:
        _tmp_HH.extend(_tmp_BH.select(station=a))
inv = _tmp_HH

#inv =  read_inventory('../Turkey.xml')
triads = maketriad.maketriad(inv, minlen=cfg['Triad']['min_side_length'], \
    maxlen=cfg['Triad']['max_side_length'], minang=cfg['Triad']['min_angle'],\
    maxang=cfg['Triad']['max_angle'])

#ot = UTCDateTime("2022-1-15T04:14:45")
ot = UTCDateTime("2023-02-06T01:17:35") 
fp = open("Tout.txt","w")
for i in range(20):
    t1= ot + 180*i
    t2 = t1 + 600
    triads.get_waveform(start=t1,end=t2,cfg=cfg)
    for td in triads:
        td.correlate(shift=3000)
        td.beamform()
    triads.polarization(method=cfg['Source'])
    a = triads.select_active()
    for i, td in enumerate(triads):
        if td.detection:
            s = "%s,%0.4f,%0.4f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%s,%0.2f\n" %(t1,td.clat, \
               td.clon, np.mean(td.cc),np.mean(td.dt), td.azm, td.polazm, td.apvel, td.beam[0],td.beam[1])
            print(s)
            fp.write(s)
fp.close()
    #triads.map_plot(center=[140,-35])