#!/usr/bin/env python3
import sys, os, argparse
from tqdm import tqdm
import requests
from obspy.core import UTCDateTime

td = UTCDateTime()
basic_url="http://service.iris.edu/fdsnws/station/1/query?"
## Parse the arguments ##
parser = argparse.ArgumentParser()
parser.add_argument("-n",type=str,default="",help="Netwotk code")
parser.add_argument("-s",type=str,default="*",help="Station code. \
Defalut is *, means all stations of the network")
parser.add_argument("-c",type=str,default="*",help="Channel code. \
Default is all chennels. Wild cards are also accepted. e.g HH?")
parser.add_argument("-l",type=str,default="station",help="level of \
information requested. Default is station. Options are station, network, channel\
 and response.")
parser.add_argument("-o",type=str,default="tmp.xml",help="Output xml")
parser.add_argument("-r",action='store_true',help="Recently active only")

args = parser.parse_args()
if args.n:
    basic_url += 'net=' + args.n +'&'
if args.s:
    basic_url += 'sta=' + args.s +'&'
if args.c:
    basic_url += 'cha=' + args.c +'&'
if args.l:
    basic_url += 'level=' + args.l +'&'
if args.r:
    basic_url += 'startbefore=' + td.isoformat()[0:22] +'&endafter='+ td.isoformat()[0:22] + '&'
basic_url += "format=xml&includecomments=true&includeavailability=true&nodata=404"
print(basic_url)
response = requests.get(basic_url, stream=True)
total_size_in_bytes= int(response.headers.get('content-length', 0))
block_size = 1024 #1 Kibibyte
progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
with open(args.o, 'wb') as file:
    for data in response.iter_content(block_size):
        progress_bar.update(len(data))
        file.write(data)
progress_bar.close()