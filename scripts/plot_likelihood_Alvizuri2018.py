#!/usr/bin/env python

import numpy as np
import subprocess

from copy import deepcopy
from likelihood_analysis import likelihood_analysis, progress
from _Alvizuri2018 import fullpath, names, depths, magnitudes,\
    data_processing_handles_FK, misfit_handles, selected_events, expected_results
from mtuq.grid import FullMomentTensorGridRegular


if __name__=='__main__':
    _i, _n = 1, len(selected_events)

    #
    # loop over selected events from Alvizuri2018
    #

    for index in selected_events:
        progress(_i, _n)

        event_id = names[index]
        depth = depths[index]
        magnitude = magnitudes[index]

        model = "MDJ2"
        solver = "FK"

        path_data, path_weights, path_greens = ( 
            fullpath(event_id, '*BH.[zrt]'), 
            fullpath(event_id, 'weights.dat'),
            "/home/rmodrak/data/FK/MDJ2",
            )

        grid = FullMomentTensorGridRegular(
            npts_per_axis=10,
            magnitudes=[magnitude],
            )

        process_bw, process_sw = data_processing_handles_FK(
            path_greens, path_weights,
            )

        misfit_bw, misfit_sw = misfit_handles(
            )

        misfit_bw_ZR = misfit_bw
        misfit_bw_ZR.time_shift_groups = ['ZR']

        misfit_sw_ZR = deepcopy(misfit_sw)
        misfit_sw_ZR.time_shift_groups = ['ZR']

        misfit_sw_T = deepcopy(misfit_sw)
        misfit_sw_T.time_shift_groups = ['T']

        likelihood_analysis(
            event_id,
            path_data,
            path_greens,
            path_weights,
            solver,
            model,
            [process_bw, process_sw, process_sw],
            [misfit_bw_ZR, misfit_sw_ZR, misfit_sw_T],
            grid,
            magnitude,
            depth)

        try:
            subprocess.call(
                'gm montage -tile 3x2 -geometry 1750x3000+0+75 ./%s_*.png ./%s.png' %
                (event_id, event_id), shell=True)
        except:
            pass

        _i += 1

