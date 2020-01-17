#!/usr/bin/env python

import numpy as np
from mtbench import benchmark, progress
from mtbench.Silwal2016 import fullpath, names, depths, magnitudes,\
    data_processing_handles, misfit_handles, selected_events, expected_results
from mtuq.grid import DoubleCoupleGridRegular


if __name__=='__main__':
    _i, _n = 1, len(selected_events)

    #
    # loop over selected events from Silwal2016
    #

    for index in selected_events:
        progress(_i, _n)

        event_id = names[index]
        depth = depths[index]
        magnitude = magnitudes[index]

        model = "ak135f_scak"
        solver = "AxiSEM"

        path_data, path_weights, path_greens = ( 
            fullpath(event_id, '*BH.[zrt]'), 
            fullpath(event_id, 'weights.dat'),
            "/home/rmodrak/data/axisem/ak135f_scak-2s",
            )

        sources = DoubleCoupleGridRegular(
            npts_per_axis=15,
            magnitude=magnitude,
            )

        process_bw, process_sw = data_processing_handles(
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

