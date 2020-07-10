#!/usr/bin/env python

import numpy as np
from mtbench import bench, progress
from _Alvizuri2018 import fullpath, names, depths, magnitudes,\
    data_processing, misfit_functions, selected_events, expected_results
from mtuq.grid import FullMomentTensorGridRandom


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

        model = "mdj2_ak135f_celso"
        solver = "AxiSEM"

        path_data, path_weights, path_greens = ( 
            fullpath(event_id, '*BH.[zrt]'), 
            fullpath(event_id, 'weights.dat'),
            "/home/rmodrak/data/axisem/mdj2_ak135f_celso-2s",
            )

        grid = FullMomentTensorGridRandom(
            npts=2000000,
            magnitudes=[magnitude],
            )

        process_data_functions = data_processing(
            path_greens, path_weights,
            )

        bench(
            event_id,
            path_data,
            path_greens,
            path_weights,
            solver,
            model,
            process_data_functions,
            misfit_functions(),
            grid,
            magnitude,
            depth)

        _i += 1

