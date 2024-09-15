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
from scipy.linalg import inv
from scipy.sparse import lil_matrix, csr_matrix
from lemkelcp.lemkelcp import lemkelcp
from data import prepare_data
from geom import MF
from phys_consts import secs_in_day
import datetime



def observation(nT, ncoefs):
    nT_add = 1
    matrix = lil_matrix((nT + nT_add, (nT + nT_add) * ncoefs)) 
    
    for i in np.arange(nT + nT_add):
        matrix[i, i*ncoefs] = 1.
    return matrix


def LCPCorrection(res, Ninv, nT, ncoefs):

    G = observation(nT, ncoefs)
    
    w = G.dot(res)
    idx = (w<0)
 
    if np.any(idx):
                            
        Gnew = G[idx,:]
        wnew = w[idx]
   
        print('constructing M')

        NGT = Ninv * Gnew.transpose()
        M = Gnew.dot(NGT)
        print('solving LCP')

        sol = lemkelcp(M,wnew,1000000)
        try:        
            res = res + NGT.dot(sol[0])
            print ('lcp adjustment done')
        except:
            print('no lcp solution found')

    return res


def normal_system(theta_station, phi_station, nT, ndays,
                            time, theta, phi, el, 
                            time_ref, theta_ref, phi_ref, el_ref, rhs):
    """
    :param nT: number of time intervals
    :param ndays: number of days
    :param time: array of times of IPPs in secs
    :param theta: array of longitudes of IPPs in rads
    :param phi: array of latitudes of IPPs in rads
    :param el: array of elevation angles in rads
    :param time_ref: array of ref times of IPPs in sec
    :param theta_ref: array of ref longitudes (LTs) of IPPs in rads
    :param phi_ref: array of ref co latitudes of IPPs in rads
    :param el_ref: array of ref elevation angles in rads
    :param rhs: array of rhs (measurements TEC difference on current and ref rays)
    """

    print('constructing normal system for series')
    tmc = time
    tmr = time_ref
    SF = MF(el)
    SF_ref = MF(el_ref)
 
    # Construct weight matrix for the observations
    len_rhs = len(rhs)
    P = lil_matrix((len_rhs, len_rhs))
    el_sin = np.sin(el)
    elr_sin = np.sin(el_ref)
    #diagP =  1. #  / (SF + SF_ref)**2 
    diagP = (el_sin ** 2) * (elr_sin ** 2) / (el_sin ** 2 + elr_sin ** 2) 

    P.setdiag(diagP)
    P = P.tocsr()
 
 
   # Construct matrix of the problem (A)
    n_coefs = 10
 
    tic = np.round(tmc * nT / (ndays * secs_in_day)).astype('int')
    tir = np.round(tmr * nT / (ndays * secs_in_day)).astype('int')
    dt = ndays * secs_in_day / nT
    
    ti_unique, ti_num_obs  = np.unique(np.append(tic, tir), return_counts=True)
    
    print (dt, tic.min(), tic.max())

    abs_tics = np.sort(np.array([i for i in range(nT + 1) if i not in ti_unique]))[::-1]
    print(abs_tics, len(ti_unique))
    
    matrix = lil_matrix((len_rhs, len(ti_unique) * n_coefs)) 

    for i in range(0, len_rhs, 1): 

        t_i_c = np.argwhere(ti_unique==tic[i])[0][0]

        matrix[i, t_i_c * n_coefs + 0] += SF[i] * 1
        matrix[i, t_i_c * n_coefs + 1] += SF[i] * (theta[i] - theta_station)
        matrix[i, t_i_c * n_coefs + 2] += SF[i] * (phi[i] - phi_station)
        matrix[i, t_i_c * n_coefs + 3] += SF[i] * (time[i] - tic[i] * dt)
        matrix[i, t_i_c * n_coefs + 4] += SF[i] * (theta[i] - theta_station)**2
        matrix[i, t_i_c * n_coefs + 5] += SF[i] * (phi[i] - phi_station)**2
        matrix[i, t_i_c * n_coefs + 6] += SF[i] * (time[i] - tic[i] * dt)**2

        matrix[i, t_i_c * n_coefs + 7] += SF[i] * (theta[i] - theta_station) * (phi[i] - phi_station)
        matrix[i, t_i_c * n_coefs + 8] += SF[i] * (phi[i] - phi_station) * (time[i] - tic[i] * dt)
        matrix[i, t_i_c * n_coefs + 9] += SF[i] * (time[i] - tic[i] * dt) * (theta[i] - theta_station)


        t_i_r = np.argwhere(ti_unique==tir[i])[0][0]

        matrix[i, t_i_r * n_coefs + 0] += -SF_ref[i] * 1
        matrix[i, t_i_r * n_coefs + 1] += -SF_ref[i] * (theta_ref[i] - theta_station)
        matrix[i, t_i_r * n_coefs + 2] += -SF_ref[i] * (phi_ref[i] - phi_station)
        matrix[i, t_i_r * n_coefs + 3] += -SF_ref[i] * (time_ref[i] - tir[i] * dt)
        matrix[i, t_i_r * n_coefs + 4] += -SF_ref[i] * (theta_ref[i] - theta_station)**2
        matrix[i, t_i_r * n_coefs + 5] += -SF_ref[i] * (phi_ref[i] - phi_station)**2
        matrix[i, t_i_r * n_coefs + 6] += -SF_ref[i] * (time_ref[i] - tir[i] * dt)**2

        matrix[i, t_i_r * n_coefs + 7] += -SF_ref[i] * (theta_ref[i] - theta_station) * (phi_ref[i] - phi_station)
        matrix[i, t_i_r * n_coefs + 8] += -SF_ref[i] * (phi_ref[i] - phi_station) * (time_ref[i] - tir[i] * dt)
        matrix[i, t_i_r * n_coefs + 9] += -SF_ref[i] * (time_ref[i] - tir[i] * dt) * (theta_ref[i] - theta_station)


    matrix = matrix.tocsr()

    # define normal system
    AP = matrix.transpose().dot(P)
    N = AP.dot(matrix).todense()    
    b = AP.dot(rhs)
    print('normal matrix (N) done')

    # # solve normal system
    Ninv = np.linalg.inv(N)     
    res = np.dot(Ninv.A, b)
    print('normal system solved')

    return res, ti_unique, ti_num_obs//2, Ninv.A



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Solve raw TECs to ')
    parser.add_argument('--in_path', 
                        default = '/home/teufel/sgpo/altboc',  
                        help='Path to data, after prepare script')
    parser.add_argument('--out_file', 
                        default = 'result.npz', 
                        help='Path to data, after prepare script')
    parser.add_argument('--short',
                        default = 150, 
                        help='minimum length of continuity interval[s]')
    parser.add_argument('--sparse',
                        default = 30, 
                        help='sparsing [s]')
    parser.add_argument('--el_cutoff',
                        default = 5, 
                        help='elevation cutoff in deg')
    parser.add_argument('--intervals',
                        default = 4 * 24, 
                        help='number of intervals per day')
    args = parser.parse_args()


    paths = args.in_path
    outputfile = args.out_file
    short = args.short 
    sparse = args.sparse 
    intervals = args.intervals
    el_cutoff = args.el_cutoff

    rhs, time, lon, lat, el, time_ref, lon_ref, lat_ref, el_ref, st_lat, st_lon, time0 = prepare_data(paths, el_cutoff=el_cutoff, short = short, sparse = sparse)

    ndays = np.ceil((np.max(time) - np.min(time)) / secs_in_day).astype('int') # number of days in input file

    nT = intervals * ndays
    ncoefs = 10
   
    res, ti_unique, ti_num_obs, Ninv = normal_system(st_lon, st_lat, nT, ndays, 
                                                     time, lon, lat, el, 
                                                     time_ref, lon_ref, lat_ref, el_ref, rhs)


    res_lcp = np.array(LCPCorrection(res, Ninv, nT=len(ti_unique)-1, ncoefs=ncoefs))

 
    time_out = np.array([time0 + datetime.timedelta(hours=h) for h in ti_unique*24.* ndays/nT])
    res_out = res_lcp[0::ncoefs]
    nobs_out = ti_num_obs

    np.savez(outputfile, time = time_out, res=res_out, nobs = nobs_out, st_lat = np.rad2deg(st_lat), st_lon = np.rad2deg(st_lon))


