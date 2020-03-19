#!/usr/bin/env python

import numpy as np
from os.path import abspath, dirname, join
from mtuq.misfit import Misfit
from mtuq.process_data import ProcessData


def moment(Mw):
    return 10.**(1.5*Mw + 9.1)



names = [
    '20061009013528000',
    '20090525005443124',
    '20130212025751273',
    '20160106013000964',
    '20160909003001386',
    '20170903033001760',
    '20170903033831810',
    '20160912113255770',
    '20171115052932820',
    ]

depths = [
    1000.0,
    1000.0,
    1000.0,
    1000.0,
    1000.0,
    1000.0,
    1000.0,
   13000.0,
    3000.0,
    ]

magnitudes = [
    3.84,
    4.38,
    4.45,
    4.25,
    4.48,
    5.26,
    4.26,
    5.50,
    5.50,
    ]


magnitudes_Fig2 = [
    3.83,
    4.39,
    4.58,
    4.48,
    4.73,
    5.26,
    ]


magnitudes_Fig2_caption = [
    3.84,
    4.38,
    4.45,
    4.25,
    4.48,
    5.18,
    ]


selected_events = [
    0,    
  ##1, 
    2,
    3,
  ##4,
    5,
  ##6,
  ##7,
  ##8,
    ]


expected_results = [
    {
      'kappa': 338.,
      'theta': 9.,
      'sigma': 65.,
      #'Mw': 3.83,
      'M0': moment(3.83),
      'gamma': -23.,
      'delta': 59.,
    },
    {
      'kappa': 145.,
      'theta': 9.,
      'sigma': 65.,
      #'Mw': 4.39,
      'M0': moment(4.39),
      'gamma': -9.,
      'delta': 48.,
    },
    {
      'kappa': 116.,
      'theta': 21.,
      'sigma': 55.,
      #'Mw': 4.58,
      'M0': moment(4.58),
      'gamma': -14.,
      'delta': 59.,
    },
    {
      'kappa': 239.,
      'theta': 25.,
      'sigma': 45.,
      #'Mw': 4.48,
      'M0': moment(4.48),
      'gamma': -14.,
      'delta': 76.,
    },
    {
      'kappa': 116.,
      'theta': 21.,
      'sigma': 75.,
      #'Mw': 4.73,
      'M0': moment(4.73),
      'gamma': -14.,
      'delta': 66.,
    },
    {
      'kappa': 136.,
      'theta': 31.,
      'sigma': 75.,
      #'Mw': 5.26,
      'M0': moment(5.26),
      'gamma': -5.,
      'delta': 76.,
    },
    ]




def data_processing_handles(
        path_greens, path_weights):

    process_bw = ProcessData(
        filter_type='Bandpass',
        period_min=1.,
        period_max=5.,
        pick_type='taup',
        taup_model='ak135',
        window_type='body_wave',
        window_length=15.,
        capuaf_file=path_weights,
        )

    process_sw = ProcessData(
        filter_type='Bandpass',
        period_min=20.,
        period_max=50.,
        pick_type='taup',
        taup_model='ak135',
        window_type='surface_wave',
        window_length=400.,
        capuaf_file=path_weights,
        apply_statics=True,
        ) 

    return process_bw, process_sw


def data_processing_handles_FK(
        path_greens, path_weights):

    process_bw = ProcessData(
        filter_type='Bandpass',
        period_min=1.,
        period_max=5.,
        pick_type='FK_metadata',
        FK_database=path_greens,
        window_type='body_wave',
        window_length=15.,
        capuaf_file=path_weights,
        )

    process_sw = ProcessData(
        filter_type='Bandpass',
        period_min=20.,
        period_max=50.,
        pick_type='FK_metadata',
        FK_database=path_greens,
        window_type='surface_wave',
        window_length=400.,
        capuaf_file=path_weights,
        apply_statics=True,
        )

    return process_bw, process_sw


def misfit_handles():
    misfit_bw = Misfit(
        time_shift_min=-2.,
        time_shift_max=+2.,
        time_shift_groups=['ZR'],
        )

    misfit_sw = Misfit(
        time_shift_min=-5.,
        time_shift_max=+5.,
        time_shift_groups=['ZR','T'],
        )

    return misfit_bw, misfit_sw


def basepath():
    import mtbench
    return abspath(dirname(mtbench.__file__))


def fullpath(*args):
    return join(basepath(), 'input/Alvizuri2018', *args)


