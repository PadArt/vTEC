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




import os
import numpy as np
from itertools import chain
from datetime import datetime
from geom import ipp, elevation
from time_utils import sec_of_day, sec_of_interval


def load_data(filepath):

    FIELDS = ['datetime', 'el', 'ipp_lat', 'ipp_lon', 'tec']
    DTYPE = (object, float, float, float, float)

    tec_suite_FIELDS = ['datetime', 'sat_x', 'sat_y', 'sat_z', 'tec']
    tec_suite_DTYPE = (object, float, float, float, float)
    dformat = "%Y-%m-%dT%H:%M:%S"
    convert = lambda x: datetime.strptime(x, dformat)



    with open(filepath, 'r') as fp:
        for l_no, line in enumerate(fp):
            if '(L, B, H)' in line:
                obs_lon, obs_lat, obs_h  = map(float, line.split(":")[1].strip().split(','))
            if '(X, Y, Z)' in line:
                obs_x, obs_y, obs_z  = map(float, line.split(":")[1].strip().split(','))

    
    data = np.genfromtxt(filepath, 
                         comments='#', 
                         names=tec_suite_FIELDS, 
                         dtype=tec_suite_DTYPE,
                         converters={"datetime": convert},  
                         )


    ip = ipp(obs_x, obs_y, obs_z, data['sat_x'], data['sat_y'], data['sat_z'])
    el = elevation(obs_x, obs_y, obs_z, data['sat_x'], data['sat_y'], data['sat_z'])


    R = np.empty((len(data['tec']),), 
                    list(zip(FIELDS,DTYPE)))

    
    R['datetime'] = data['datetime']
    R['el'] = el

    R['ipp_lat'] = ip[1]
    R['ipp_lon'] = ip[0]
    R['tec'] = data['tec']

    
    return R, obs_lat, obs_lon



def getContInt(time, tec, lon, lat, el,  el_cutoff=30, maxgap=30, maxjump=1):
    r = np.array(range(len(time)))
    idx = np.isfinite(tec) & np.isfinite(lon) & np.isfinite(lat) & np.isfinite(el) & (el > np.deg2rad(el_cutoff)) & (tec != 0.)
    r = r[idx]
    intervals = []
    if len(r) == 0:
        return intervals
    beginning = r[0]
    last = r[0]
    last_time = time[last]
    for i in r[1:]:
        if abs(time[i] - last_time) > maxgap or abs(tec[i] - tec[last]) > maxjump:
            intervals.append((beginning, last))
            beginning = i
        last = i
        last_time = time[last]
        if i == r[-1]:
            intervals.append((beginning, last))
    return idx, intervals


def prepare_data(paths, el_cutoff = 30, short = 600, sparse = 90):

    Atec = Atime = Along = Alat = Ael = Atime_ref = Along_ref = Alat_ref = Ael_ref = np.array([])


    for subdir, dirs, files in os.walk(paths):

        for file in files:

            filepath = subdir + os.sep + file

            if filepath.endswith(".dat"):
                print (filepath)

                try:            
                    data, obs_lat, obs_lon = load_data(filepath)
                    tt = sec_of_day(data['datetime'])                
                    idx, intervals = getContInt(tt, data['tec'], data['ipp_lon'], data['ipp_lat'], data['el'], el_cutoff=el_cutoff,  maxgap=35., maxjump=0.5)



                    for ii in intervals:
                        if (tt[ii[1]] - tt[ii[0]]) >  short:     
                            tec_out = data['tec'][ii[0]:ii[1]] 
                            time_out = data['datetime'][ii[0]:ii[1]]
                            ipp_lon_out = data['ipp_lon'][ii[0]:ii[1]]
                            ipp_lat_out = data['ipp_lat'][ii[0]:ii[1]]
                            el_out = data['el'][ii[0]:ii[1]]

                            ind_sparse = (tt[ii[0]:ii[1]] % sparse == 0)


                            tec_out = tec_out[ind_sparse]

                            time_out = time_out[ind_sparse]
                            ipp_lon_out = ipp_lon_out[ind_sparse]
                            ipp_lat_out = ipp_lat_out[ind_sparse]
                            el_out = el_out[ind_sparse]

                            dtec = tec_out[1:] - tec_out[0:-1]
                            time_out_ref = time_out[0:-1]
                            time_out = time_out[1:]
                            ipp_lon_out_ref = ipp_lon_out[0:-1]
                            ipp_lon_out = ipp_lon_out[1:]
                            ipp_lat_out_ref = ipp_lat_out[0:-1]
                            ipp_lat_out = ipp_lat_out[1:]
                            el_out_ref = el_out[0:-1]
                            el_out = el_out[1:]
                       
                            Atec = np.append(Atec, dtec)
                            Atime = np.append(Atime, time_out)
                            Along = np.append(Along, ipp_lon_out)
                            Alat = np.append(Alat, ipp_lat_out)
                            Ael = np.append(Ael, el_out)
                            Atime_ref = np.append(Atime_ref, time_out_ref)
                            Along_ref = np.append(Along_ref, ipp_lon_out_ref)
                            Alat_ref = np.append(Alat_ref, ipp_lat_out_ref)
                            Ael_ref = np.append(Ael_ref, el_out_ref)



                      
                        else: 
                            print('too short interval')

                except:
                    print('warning')

    t = np.min(Atime)
    time0 = t.replace(hour=0, minute=0, second=0, microsecond=0)
    print ('number of observations', len(Atec))
    return Atec, sec_of_interval(Atime, time0), Along, Alat, Ael, sec_of_interval(Atime_ref, time0), Along_ref, Alat_ref, Ael_ref, np.deg2rad(obs_lat), np.deg2rad(obs_lon), time0

