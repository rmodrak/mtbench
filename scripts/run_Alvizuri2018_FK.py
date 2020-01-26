#!/usr/bin/env python

import numpy as np
from mtbench import benchmark, progress
from mtbench.Alvizuri2018 import fullpath, names, depths, magnitudes_Fig2,\
    data_processing_handles_FK, misfit_handles, selected_events, expected_results
from mtuq.grid import FullMomentTensorGridRegular

magnitudes = magnitudes_Fig2


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

        sources = FullMomentTensorGridRegular(
            npts_per_axis=15,
            magnitudes=[magnitude],
            )

        process_bw, process_sw = data_processing_handles_FK(
            path_greens, path_weights,
            )

        misfit_bw, misfit_sw = misfit_handles(
            )

        benchmark(
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
            sources,
            magnitude,
            depth)

        _i += 1

