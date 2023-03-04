#!/usr/bin/env python

import numpy as np
from mtbench import bench, progress
from _Alvizuri2018 import fullpath, names, depths, magnitudes,\
    data_processing_FK, misfit_functions, selected_events, expected_results
from mtuq.grid import FullMomentTensorGridSemiregular


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

        grid = FullMomentTensorGridSemiregular(
            npts_per_axis=10,
            magnitudes=[magnitude],
            )

        process_bw, process_sw = data_processing_FK(
            path_greens, path_weights,
            )

        bench(
            event_id,
            path_data,
            path_greens,
            path_weights,
            solver,
            model,
            grid,
            magnitude,
            depth,
            process_bw,
            process_sw,
            include_bw=False,
            include_rayleigh=True,
            include_love=True,
            include_mt=True,
            include_force=False,
            plot_waveforms=True,
            lune_misfit=True,
            )

        _i += 1

