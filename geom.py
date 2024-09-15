#!/usr/bin/env python
# coding=utf8
#
# Copyright 2024 Artem Padokhin <padokhin@physics.msu.ru>
#
# This file is part of vTEC.
#
# vTEC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# vTEC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with vTEC.  If not, see <http://www.gnu.org/licenses/>.



import numpy as np
from phys_consts import RE, IPPh

def xyz2spher(x,y,z):
    '''
    ECEF X,Y,Z [m] to spherical L(longitude), B(latitude) [rad], H (height above sphere) [m]
    '''
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    
    Q = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    # H - height
    H = Q - RE

    # L - longitude
    L = np.arctan2(y, x)
    L[L < 0] = L[L < 0] + 2. * np.pi

    # B - spherical latitude
    B =  np.pi / 2. - np.arccos(z / Q)
    return L, B, H


def MF(el, IPPh=IPPh):
    """
    :param el: elevation angle in [rad]
    :param IPPh: height of IPP in [m]
    """
    return 1./np.sqrt(1 - (RE * np.cos(el) / (RE + IPPh)) ** 2)


def ipp(obs_x, obs_y, obs_z, sat_x, sat_y, sat_z, IPPh=IPPh):
    """
    :param obs_x: observer x coord in ECEF [m]
    :param obs_y: observer y coord in ECEF [m]
    :param obs_z: observer z coord in ECEF [m]
    :param sat_x: satellite x coord in ECEF [m]
    :param sat_y: satellite y coord in ECEF [m]
    :param sat_z: satellite z coord in ECEF [m]
    :param IPPh: height of IPP in meters
    """
    A = (sat_x - obs_x) ** 2 + (sat_y - obs_y) ** 2 + (sat_z - obs_z) ** 2 
    obs2 = obs_x **2 + obs_y **2 + obs_z **2    
    B = 2 * (obs_x * sat_x + obs_y * sat_y + obs_z * sat_z - obs2)       
    C = obs2 - (RE + IPPh)**2       
    D = B**2 - 4*A*C

    t1 = (- B + np.sqrt(D)) / (2. * A)
    t2 = (- B - np.sqrt(D)) / (2. * A)

    idx1 = np.where((t1<=1)&(t1>=0))
    idx2 = np.where((t2<=1)&(t2>=0))

    t = np.empty(len(sat_x))

    t[idx1] = t1[idx1]
    t[idx2] = t2[idx2]


    ipp_x = obs_x * (1 - t) + sat_x * t
    ipp_y = obs_y * (1 - t) + sat_y * t
    ipp_z = obs_z * (1 - t) + sat_z * t
    
    ipp_l, ipp_b, ipp_h = xyz2spher(ipp_x,ipp_y,ipp_z)
    return ipp_l, ipp_b, ipp_h


def elevation(obs_x, obs_y, obs_z, sat_x, sat_y, sat_z):
    '''
    :param obs_x: observer x coord in ECEF [m]
    :param obs_y: observer y coord in ECEF [m]
    :param obs_z: observer z coord in ECEF [m]
    :param sat_x: satellite x coord in ECEF [m]
    :param sat_y: satellite y coord in ECEF [m]
    :param sat_z: satellite z coord in ECEF [m]
    '''   
    A = (sat_x - obs_x) ** 2 + (sat_y - obs_y) ** 2 + (sat_z - obs_z) ** 2 
    C = obs_x **2 + obs_y **2 + obs_z **2    
    B = obs_x * sat_x + obs_y * sat_y + obs_z * sat_z - C

    el = np.arcsin(B / (np.sqrt (A * C)))
    return el



