# vTEC

Single GNSS Station Absolute Vertical Total Electron Content Estimator. 


## About

This code implements phase-difference approach to reconstruct absolute vertical TEC from single station GNSS dual frequency phase observations and does not require satellite/receiver DCBs estimations and thus makes it easier to combine different GNSS systems within single reconstruction algorithm. This  method uses the representation of the ionosphere as a thin layer with TEC variations in the vicinity of the station given by truncated Taylor series up to second order in space and time. Slant to vertical TEC conversion is done by SLM mapping function. The expansion coefficients are determined by least squares with inequality constrains representing the positivity of TEC, implemented by solving the corresponding linear complementarity problem (LCP). 

## Prerequirements 

Code and launch were tested for Linux (Mint) and `python3.6.9`
Code uses `numpy`, `scipy` and `lemkelcp`  packages.

## Use test data

Input data for `vTEC` is the output data (in TXT format) from the [`tec-suite`](https://github.com/gnss-lab/tec-suite) package developed by [`SIMuRG`](https://simurg.space/) team. To use with `vTEC` you need to change default output of `tec-suite` in its configuration file `tecs.cfg` as following:

    recFields = 'datetime, sat.x, sat.y, sat.z, tec.l1l2'

## Run processing

To run processing script simply execute: 

    python3 vTECv.py --in_path `data_path` --out_file `res_file` --short 150 --sparse 30 --el_cutoff 5 --intervals 96

`data_path` is the directory of `tec-suite` output, `res_file` is a file in *.npz format to store results, with result.npz set as a default, `short` key sets the minimum accepted length of data continuity interval in seconds with default value 150, `sparse` key sets sparsing of data in seconds with default 30,  `el_cutoff` is minimal accepted satellite elevation in degrees with default 5, and `intervals` is number of estimated vTEC values per day.
