import numpy as np
import copy
from obspy import UTCDateTime, Stream
from obspy.clients.filesystem.sds import Client as SDS
from geographiclib.geodesic import Geodesic
class Triad(object):
    def __init__(self, stations=[]):
        if len(stations) != 3:
            return None
        self.stations = stations
        self.dist, self.ang = self.calc_dis_ang()
        self.clat, self.clon, self.CX, self.CY, self.Amat = self.centroid()
        self.cc = np.zeros(3)
        self.dt = np.zeros(3)
        self.apvel = 0.0
        self.azm = 0.0
        self.detection = False
        self.have_data = False
        self.data = None
    
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
        from obspy.signal import cross_correlation
        if self.have_data:
            a = cross_correlation.correlate(self.data[0].data,self.data[1].data,shift=shift)
            b = cross_correlation.correlate(self.data[0].data,self.data[2].data,shift=shift)
            c = cross_correlation.correlate(self.data[1].data,self.data[2].data,shift=shift)
            self.dt[0], self.cc[0] = cross_correlation.xcorr_max(a)
            self.dt[1], self.cc[1] = cross_correlation.xcorr_max(b)
            self.dt[2], self.cc[2] = cross_correlation.xcorr_max(c)
            if np.amax(self.cc) >= 0.5 and np.sum(self.dt) <=50:
                alpha = np.linalg.lstsq(self.Amat, self.dt,rcond=None)[0]
                self.apvel = 1.0/np.sqrt(np.sum(alpha**2))
                self.azm = np.arctan2(alpha[0],alpha[1]) * 180.0/np.pi
                self.detection = True
        else:
            pass



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
    
    def select_active(self):
        active_triads = []
        for td in self.triads:
            if td.detection:
                active_triads.append(td)
        return self.__class__(triads=active_triads)

    def get_waveform(self, start=None, end=None, source=None, target_sps=1):
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
            if st.count() == 3:
                st.detrend("demean")
                st.detrend("linear")
                st.filter('bandpass', freqmin=1/250, freqmax=1/20)
                st.resample(target_sps,window='cosine')
                st.normalize()
                td.data = st
                td.have_data = True
            else:
                continue