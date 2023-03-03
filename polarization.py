import numpy as np
from scipy.signal import hilbert


def corr_corff(st, da=1):
    z = st.select(channel="??Z")[0].data
    n = st.select(channel="??N")[0].data
    e = st.select(channel="??E")[0].data
    zh = np.imag(hilbert(z))
    zz = np.correlate(zh,zh)
    czr = []
    aczr = []
    for theta in range(0,360,da):
        tmp = theta * np.pi/180
        m11 = np.cos(tmp)
        m12 = np.sin(tmp)
        r = n*m11 + e*m12
        #t = n*-m12 + e*m11
        zr = np.correlate(zh,r)
        rr = np.correlate(r,r)
        czr.append(zr/np.sqrt(rr*zz))
        aczr.append(zr/zz)
    x = np.array(czr)
    y = np.array(aczr)
    baz =  np.argmax(x) * da
    return baz