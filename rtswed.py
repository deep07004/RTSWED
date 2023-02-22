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
#inv =  read_inventory('../Turkey.xml')
triads = maketriad.maketriad(inv, minlen=10, maxlen=1200, minang=10, maxang=120)

t1 = UTCDateTime("2022-1-15T04:14:45")
#t1 = UTCDateTime("2023-02-06T01:17:35") + 400
t2= t1+3600

triads.get_waveform(start=t1,end=t2,source = ['sds','/home/deep/CWP/Landslide/Tonga/buffer'])
#triads.get_waveform(start=t1,end=t2,source = ['sds','/home/deep/CWP/Landslide/Turkey'])

for td in triads:
    td.correlate(shift=3000)
    td.beamform()
a = triads.select_active()
for i, td in enumerate(triads):
    if td.detection:
        s = "%d %0.4f %0.4f %0.2f %0.2f %s %0.2f\n" %(i,td.clat, td.clon, td.azm, td.apvel, td.beam[0],td.beam[1])
        print(s)
triads.map_plot(center=[140,-35])