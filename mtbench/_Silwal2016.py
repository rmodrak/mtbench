#!/usr/bin/env python

import numpy as np

from os.path import abspath, dirname, join
from mtuq.misfit import Misfit
from mtuq.process_data import ProcessData


names = [
    '20070911234634153', #0
    '20070919112226549', #1
    '20071003140612444', #2
    '20071010180326301', #3
    '20071128235703849', #4
    '20080314093821771', #5
    '20080327230745201', #6
    '20080828231418631', #7
    '20080918194353069', #8
    '20081228071310738', #9
    '20090124180950811', #10
    '20090215193500098', #11
    '20090223000427175', #12
    '20090317011333066', #13
    '20090407201255351', #14
    '20090414171427415', #15
    '20090430045457938', #16
    '20090524094004552', #17
    '20090622192805162', #18
    '20090626164820729', #19
    '20090730223910267', #20
    ]

depths = [
    94000.,
    47000.,
    32000.,
    27000.,
    69000.,
   139000.,
    65000.,
    54000.,
    81000.,
    82000.,
   105000.,
    43000.,
    81000.,
    96000.,
    39000.,
   111000.,
    40000.,
   109000.,
    62000.,
    56000.,
    60000.,
    ]

magnitudes = [
    4.40,
    4.50,
    5.00,
    4.20,
    4.90,
    5.10,
    5.10,
    4.20,
    4.60,
    4.60,
    5.80,
    4.50,
    4.90,
    4.30,
    4.50,
    4.30,
    4.80,
    4.60,
    5.40,
    4.20,
    4.60,
    ]

selected_events = [
   0,    
   1, 
   2,
  #3,
   4,
   5,
  #6,
   7,
   8,
   9,
   12,
   13,
   14,
  #16,
   20,
    ]


expected_results = None


def data_processing(
        path_greens, path_weights):

    process_bw = ProcessData(
        filter_type='Bandpass',
        freq_min= 0.25,
        freq_max= 0.667,
        pick_type='taup',
        taup_model='ak135',
        window_type='body_wave',
        window_length=15.,
        capuaf_file=path_weights,
        ) 

    process_sw = ProcessData(
        filter_type='Bandpass',
        freq_min=0.025,
        freq_max=0.0625,
        pick_type='taup',
        taup_model='ak135',
        window_type='surface_wave',
        window_length=120.,
        capuaf_file=path_weights,
        ) 

    return (process_bw, process_sw)


def data_processing_FK(
        path_greens, path_weights):

    process_bw = ProcessData(
        filter_type='Bandpass',
        freq_min= 0.25,
        freq_max= 0.667,
        pick_type='FK_metadata',
        FK_database=path_greens,
        window_type='body_wave',
        window_length=15.,
        capuaf_file=path_weights,
        )

    process_sw = ProcessData(
        filter_type='Bandpass',
        freq_min=0.025,
        freq_max=0.0625,
        pick_type='FK_metadata',
        FK_database=path_greens,
        window_type='surface_wave',
        window_length=120.,
        capuaf_file=path_weights,
        )

    return process_bw, process_sw


def misfit_functions():
    misfit_bw = Misfit(
        time_shift_min=-2.,
        time_shift_max=+2.,
        time_shift_groups=['ZR'],
        )

    misfit_sw = Misfit(
        time_shift_min=-10.,
        time_shift_max=+10.,
        time_shift_groups=['ZR','T'],
        )

    return misfit_bw, misfit_sw


def basepath():
    import mtbench
    return abspath(join(dirname(mtbench.__file__), '..'))


def fullpath(*args):
    return join(basepath(), 'WAVEFORMS/Silwal2016', *args)


