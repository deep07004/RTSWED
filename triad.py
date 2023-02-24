import numpy as np
import copy
from obspy import UTCDateTime, Stream
from obspy.clients.filesystem.sds import Client as SDS
from obspy.clients.fdsn.client import Client
from geographiclib.geodesic import Geodesic
from obspy.signal import cross_correlation
from obspy.signal.filter import envelope

def request_data(i, td, start, end, sps,fl,fh, method):
    st = Stream()
    streams = [(sta[0].split('.')[0], sta[0].split('.')[1],'*',"?HZ",start, end) for sta in td.stations]

    if method[0] == "SDS":
        sds = SDS(method[1])
        try:
            st += sds.get_waveforms_bulk(streams)
        except:
            pass
    if method[0] == "FDSNWS":
        fdsnws = Client(method[1])
        try:
            st += fdsnws.get_waveforms_bulk(streams)
        except:
            pass

    if st.count() > 0:
        st.merge(fill_value=0)
        st.detrend("demean")
        st.detrend("linear")
        st.filter('bandpass', freqmin=fl, freqmax=fh)
        st.resample(sps, window='cosine')
        sss = [tr.stats.starttime for tr in st]
        eee = [tr.stats.endtime for tr in st]
        sss.sort(reverse=True)
        eee.sort()
        t1 = sss[0] + (1000000-sss[0].microsecond)/1000000
        t2 = eee[0] - eee[0].microsecond/1000000
        if (t2-t1) < 100 :
            st.clear()
        else:
            st.trim(t1,t2)
    
    tr = Stream()
    for sta in td.stations:
        s = sta[0].split('.')[1]
        try:
            tr +=  st.select(station=s).sort(reverse=True)[0]
        except:
            print("No data for stations %s" %sta)
            return (-1,-1)
    npts = [x.stats.npts for x in tr]
    st = tr.copy()
    if st.count() == 3 and np.mean(npts) == npts[0]:
        st.normalize()
        for tr in st:
            tr.data = envelope(tr.data)
        return (i, st)
    else:
        return (-1,-1)


class Triad(object):
    def __init__(self, stations=[]):
        if len(stations) != 3:
            return None
        self.stations = stations
        self.dist, self.ang = self.calc_dis_ang()
        self.clat, self.clon, self.CX, self.CY, self.Amat = self.centroid()
        self.cc = np.zeros(3)
        self.dt = np.zeros(3)
        self.alpha = 0.0
        self.apvel = 0.0
        self.azm = 0.0
        self.detection = False
        self.have_data = False
        self.data = None
        self.beam = None
        self.BFS = None
    
    def calc_dis_ang(self):
        gl = Geodesic.WGS84
        sta = self.stations
        dist = np.zeros(3)
        ang = np.zeros(3)
        dist[0] = gl.Inverse(sta[0][1],sta[0][2],sta[1][1],sta[1][2])['s12']/1000.0
        dist[1] = gl.Inverse(sta[0][1],sta[0][2],sta[2][1],sta[2][2])['s12']/1000.0
        dist[2] = gl.Inverse(sta[1][1],sta[1][2],sta[2][1],sta[2][2])['s12']/1000.0
        ang[0] = np.arccos((dist[1]**2 + dist[2] **2 - dist[0] **2) / (2 * dist[1] * dist[2]))
        ang[1] = np.arccos((dist[0]**2 + dist[2] **2 - dist[1] **2) / (2 * dist[0] * dist[2]))
        ang[2] = np.arccos((dist[0]**2 + dist[1] **2 - dist[2] **2) / (2 * dist[0] * dist[1]))
        ang *= 180.0/np.pi
        return dist, ang
    def centroid(self):
        gl = Geodesic.WGS84
        R = gl.a
        sta = self.stations
        lat = np.array([90-sta[i][1] for i in range(3)])
        lon = np.array([sta[i][2] for i in range(3)])
        lat *= np.pi/180
        lon *= np.pi/180
        X = R*np.sin(lat)*np.cos(lon)
        Y = R*np.sin(lat)*np.sin(lon)
        Z = R*np.cos(lat)
        X3DC = X.mean()
        Y3DC = Y.mean()
        Z3DC = Z.mean()
        L = np.sqrt(X3DC**2 + Y3DC**2 + Z3DC**2)
        XC = X3DC/L
        YC = Y3DC/L
        ZC = Z3DC/L
        clat = 90 - (np.arctan2(np.sqrt(XC**2 + YC**2),ZC))*180/np.pi
        clon = np.arctan2(YC, XC) * 180/np.pi
        # Position of the stations in a cartesian coordinate system with centroid as origin.
        cosref = np.cos(clat)
        cx = np.array([111.2*(sta[i][2]-clon)*cosref for i in range(3)])
        cy = np.array([111.2*(sta[i][1]-clat) for i in range(3)])
        Amat = np.zeros([3,2])
        Amat[0,0] = cx[1] - cx[0]
        Amat[1,0] = cx[2] - cx[0]
        Amat[2,0] = cx[2] - cx[1]
        Amat[0,1] = cy[1] - cy[0]
        Amat[1,1] = cy[2] - cy[0]
        Amat[2,1] = cy[2] - cy[1]
        return clat, clon, cx, cy, Amat
    
    def copy(self):
        return copy.deepcopy(self)
    
    def correlate(self,shift=300):
        if self.have_data:
            sps = self.data[0].stats.sampling_rate
            shift_sps = int(shift * sps)
            a = cross_correlation.correlate(self.data[0].data,self.data[1].data,shift=shift_sps)
            b = cross_correlation.correlate(self.data[0].data,self.data[2].data,shift=shift_sps)
            c = cross_correlation.correlate(self.data[1].data,self.data[2].data,shift=shift_sps)
            self.dt[0], self.cc[0] = cross_correlation.xcorr_max(a)
            self.dt[1], self.cc[1] = cross_correlation.xcorr_max(b)
            self.dt[2], self.cc[2] = cross_correlation.xcorr_max(c)
            self.dt /= sps
            if np.amax(self.cc) >= 0.5 and abs(np.sum(self.dt)) <=50:
                alpha = np.linalg.lstsq(self.Amat, self.dt,rcond=None)[0]
                self.alpha = alpha
                self.apvel = 1.0/np.sqrt(np.sum(alpha**2))
                azm = np.arctan2(alpha[0],alpha[1]) * 180.0/np.pi + 180
                if azm > 180:
                    azm -= 360.0
                self.azm = azm
                self.detection = True
        else:
            pass
    def beamform(self):
        if self.detection:
            tr = self.data[0].copy()
            sps = self.data[0].stats.sampling_rate
            t1 = np.ceil((self.CX[0]*self.alpha[0] + self.CY[0]*self.alpha[1]) * sps) 
            t2 = np.ceil((self.CX[1]*self.alpha[0] + self.CY[1]*self.alpha[1]) * sps)
            t3 = np.ceil((self.CX[2]*self.alpha[0] + self.CY[2]*self.alpha[1]) * sps)
            tr.stats.station = "BFS"
            tr.data = np.roll(self.data[0].data, int(t1)) + np.roll(self.data[1].data, int(t2)) + np.roll(self.data[2].data, int(t3))
            ii = np.argmax(tr.data)
            self.beam = [tr.stats.starttime + tr.stats.delta * ii, tr.data[ii]]
            self.BFS = tr
    def get_waveform(self,start, end, sps,fl,fh, method):
        st = Stream()
        streams = [(sta[0].split('.')[0], sta[0].split('.')[1],'*',"?HZ", start, end) for sta in self.stations]
        if method[0] == "SDS":
            sds = SDS(method[1])
            try:
                st += sds.get_waveforms_bulk(streams)
            except:
                pass
        if method[0] == "FDSNWS":
            fdsnws = Client(method[1])
            try:
                st += fdsnws.get_waveforms_bulk(streams)
            except:
                pass
        if st.count() > 0:
            st.merge(fill_value=0)
            st.detrend("demean")
            st.detrend("linear")
            st.filter('bandpass', freqmin=fl, freqmax=fh)
            st.resample(sps, window='cosine')
            sss = [tr.stats.starttime for tr in st]
            eee = [tr.stats.endtime for tr in st]
            sss.sort(reverse=True)
            eee.sort()
            t1 = sss[0] + (1000000-sss[0].microsecond)/1000000
            t2 = eee[0] - eee[0].microsecond/1000000
            if (t2-t1) < 100 :
                st.clear()
            else:
                st.trim(t1,t2)

        tr = Stream()
        for sta in self.stations:
            s = sta[0].split('.')[1]
            try:
                tr +=  st.select(station=s).sort(reverse=True)[0]
            except:
                print("No data for stations %s" %sta)
                return (-1,-1)
        npts = [x.stats.npts for x in tr]
        st = tr.copy()
        if st.count() == 3 and np.mean(npts) == npts[0]:
            st.normalize()
            for tr in st:
                tr.data = envelope(tr.data)
            self.data = st
        else:
            self.data = None

    def reset(self):
        self.cc = np.zeros(3)
        self.dt = np.zeros(3)
        self.alpha = 0.0
        self.apvel = 0.0
        self.azm = 0.0
        self.detection = False
        self.have_data = False
        self.data = None
        self.beam = None
        self.BFS = None

class Triads(object):
    """
    A list object of Triad with functions for processing and manipulation.
    """
    def __init__(self, triads=None):
        self.triads = []
        if isinstance(triads, Triad):
            triads = [triads]
        if triads:
            self.triads.extend(triads)
    def __iter__(self):
        """
        """
        return list(self.triads).__iter__()

    def __nonzero__(self):
        """
        """
        return bool(len(self.triads))
    def __len__(self):
        """
        Return the number of Triads in the Stream object.
        """
        return len(self.triads)
    
    count = __len__

    def __setitem__(self, index, triad):
        """
        __setitem__ method of obspy.Stream objects.
        """
        self.traces.__setitem__(index, triad)

    def __getitem__(self, index):
        """
        __getitem__ method of obspy.Stream objects.
        :return: Trace objects
        """
        if isinstance(index, slice):
            return self.__class__(triads=self.triads.__getitem__(index))
        else:
            return self.triads.__getitem__(index)

    def __delitem__(self, index):
        """
        Passes on the __delitem__ method to the underlying list of traces.
        """
        return self.triads.__delitem__(index)

    def __getslice__(self, i, j, k=1):
        """
        __getslice__ method of obspy.Stream objects.
        :return: Stream object
        """
        # see also https://docs.python.org/3/reference/datamodel.html
        return self.__class__(triads=self.triads[max(0, i):max(0, j):k])

    def append(self, triad):
        if isinstance(triad, Triad):
            self.triads.append(triad)
        else:
            msg = 'Append only supports a single Triad object as an argument.'
            raise TypeError(msg)
        return self
    
    def copy(self):
        return copy.deepcopy(self)
    
    def extend(self, triad_list):
        if isinstance(triad_list, list):
            for _i in triad_list:
                # Make sure each item in the list is a triad.
                if not isinstance(_i, Triad):
                    msg = 'Extend only accepts a list of Triad objects.'
                    raise TypeError(msg)
            self.triads.extend(triad_list)
        elif isinstance(triad_list, Triads):
            self.triads.extend(triad_list.triads)
        else:
            msg = 'Extend only supports a list of Triad objects as argument.'
            raise TypeError(msg)
        return self
    
    def select_active(self, apvel=None, azm=None):
        if not apvel:
            apvel=[0,100]
        if not azm:
            azm=[0,360] 
        active_triads = []
        for td in self.triads:
            if td.detection and td.apvel >= apvel[0] and td.apvel <= apvel[1] \
                and td.azm >= azm[0] and td.azm<= azm[1]:
                active_triads.append(td)
        return self.__class__(triads=active_triads)

    def get_waveform__(self, start=None, end=None, source=None, target_sps=1):
        """
        1. Request data from source.
        2. Preprocess: filter, resample and align.
        """
        if not start:
            end = UTCDateTime.now() -10
            start = end - 300
        if not source:
            raise ValueError("No data source defined")
        if source[0] == "sds":
            sds_arch = source[1]
        sds = SDS(sds_arch)
        for td in self.triads:
            st = Stream()
            for sta in td.stations:
                s = sta[0].split('.')
                try:
                    st += sds.get_waveforms(s[0],s[1],"*", "?HZ", start, end)[0]
                except:
                    print("No data for ", s[0],s[1],"*", "?HZ", start, end)
                    continue
            if st.count() > 0:
                st.merge(fill_value=0)
                st.detrend("demean")
                st.detrend("linear")
                st.filter('bandpass', freqmin=1/250, freqmax=1/20)
                st.resample(target_sps,window='cosine')
                sss = [tr.stats.starttime for tr in st]
                eee = [tr.stats.endtime for tr in st]
                sss.sort(reverse=True)
                eee.sort()
                t1 = sss[0] + (1000000-sss[0].microsecond)/1000000
                t2 = eee[0] - eee[0].microsecond/1000000
                if (t2-t1) < 100 :
                    continue
                else:
                    st.trim(t1,t2)
            if st.count() == 3:
                st.normalize()
                #for tr in st:
                #    tr.data = envelope(tr.data)
                td.data = st
                td.have_data = True
            else:
                continue
    def get_waveform(self, start=None, end=None, source=None, cfg=None):
        """
        1. Request data from source.
        2. Preprocess: filter, resample and align.
        """
        import concurrent.futures
        from tqdm import tqdm

        if cfg == None:
            raise ValueError("No configuration")
        fl = cfg['Filter'][0]
        fh = cfg['Filter'][1]

        if not start:
            end = UTCDateTime.now() -10
            start = end - 300
        processes = []
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for i, td in enumerate(self.triads):
               
                p = executor.submit(request_data,i, td, start, end, cfg['target_sps'],fl,fh, cfg['Source'])
                processes.append(p)
            pbar = tqdm(total=len(processes), desc='Overall progress')
            for f in concurrent.futures.as_completed(processes):
                jj = f.result()[0]
                if jj == -1:
                    continue
                else:
                    self.triads[jj].data = f.result()[1]
                    self.triads[jj].have_data = True
                pbar.update(n=1)
    def map_plot(self,center=[0.0,0.0],x=20,y=20):
        import pygmt
        import pandas as pd
        proj = "A%f/%f/5i" %(center[0],center[1])
        fig = pygmt.Figure()
        fig.basemap(region="g", projection=proj, frame=True)
        fig.coast(shorelines=True)
        sta= []
        for td in self.triads:
            if td.detection and (td.apvel > 2.5 and td.apvel<5.5):
                sta.append([td.clon,td.clat,td.azm, 1.0])
        s = pd.DataFrame(sta)
        fig.plot(x=center[0],y=center[1],style="a0.8c",color="red")
        if len(s) > 1:
            fig.plot(x=s[0],y=s[1], style="V0.2c+ea+bc",direction=[s[2],s[3]],pen="0.5p",color="cyan" )
        fig.savefig('triad.jpg')
        fig.show()