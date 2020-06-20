#!/usr/bin/env python

import numpy as np
from mtbench import run_grid_search, progress
from _Alvizuri2018 import fullpath, names, depths, magnitudes,\
    data_processing_handles, misfit_handles, selected_events, expected_results
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

        model = "ak135"
        solver = "syngine"

        path_data, path_weights, path_greens = ( 
            fullpath(event_id, '*BH.[zrt]'), 
            fullpath(event_id, 'weights.dat'),
            "http://service.iris.edu/irisws/syngine/1",
            )

        grid = FullMomentTensorGridSemiregular(
            npts_per_axis=15,
            magnitudes=[magnitude],
            )

        process_bw, process_sw = data_processing_handles(
            path_greens, path_weights,
            )

        misfit_bw, misfit_sw = misfit_handles(
            )

        run_grid_search(
            event_id,
            path_data,
            path_greens,
            path_weights,
            solver,
            model,
            process_bw,
            process_sw,
            misfit_bw,
            misfit_sw,
            grid,
            magnitude,
            depth)

        _i += 1

