import numpy as np
import pandas as pd
from geographiclib.geodesic import Geodesic
class Triad(object):
    def __init__(self, stations=[]):
        if len(stations) != 3:
            return None
        self.stations = stations
        self.dist, self.ang = self.calc_dis_ang()
        self.cc = np.zeros(3)
        self.dt = np.zeros(3)

    
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

class Triads(object):
    """
    A list object of Triad
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
    def append(self, triad):
        if isinstance(triad, Triad):
            self.triads.append(triad)
        else:
            msg = 'Append only supports a single Triad object as an argument.'
            raise TypeError(msg)
        return self
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