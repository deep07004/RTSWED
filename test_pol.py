import numpy as np
from scipy.signal import hilbert
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from geographiclib.geodesic import Geodesic
import matplotlib.pylab as plt
import sys

gl = Geodesic.WGS84
def corr_corff (st):
    z = st.select(channel="??Z")[0].data
    n = st.select(channel="??N")[0].data
    e = st.select(channel="??E")[0].data
    zh = np.imag(hilbert(z))
    zz = np.correlate(zh,zh)
    czr = []
    aczr = []
    for theta in range(0,360):
        tmp = theta * np.pi/180
        m11 = np.cos(tmp)
        m12 = np.sin(tmp)
        r = n*m11 + e*m12
        #t = n*-m12 + e*m11
        zr = np.correlate(zh,r)
        rr = np.correlate(r,r)
        czr.append(zr/np.sqrt(rr*zz))
        aczr.append(zr/zz)
    return np.array(czr), np.array(aczr)


OT = UTCDateTime("2023-02-06T01:17:35") 
evla = 37.2251
evlo = 37.0209
net =  sys.argv[1]
sta =  sys.argv[2]
loc =  sys.argv[3]
channel = sys.argv[4]
stla = float(sys.argv[5])
stlo = float(sys.argv[6])
print(evla,evlo,stla,stlo)
distaz = gl.Inverse(stla,stlo,evla,evlo)
or_baz = distaz['azi1']
distkm = distaz['s12']/1000.0

cl = Client("IRIS")
t1 = OT + distkm/4 - 20.0
t2 = t1 + 300
st = cl.get_waveforms(net,sta, loc, channel, t1,t2)
st.select(channel="BH2")[0].stats.channel ="BHE"
st.select(channel="BH1")[0].stats.channel ="BHN"
st.filter('bandpass',freqmin=0.02, freqmax=0.04)
st.plot()
x, y = corr_corff(st)
print(np.argmax(x),np.argmax(y),or_baz)
plt.plot(x)
plt.plot(y)
plt.show()
