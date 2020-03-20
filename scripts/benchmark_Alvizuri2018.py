#!/usr/bin/env python

import numpy as np
from mtbench import run_grid_search, progress
from _Alvizuri2018 import fullpath, names, depths, magnitudes,\
    data_processing_handles_FK, misfit_handles, selected_events, expected_results
from mtuq.grid import UnstructuredGrid


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

        source = expected_results[index]

        model = "MDJ2"
        solver = "FK"

        path_data, path_weights, path_greens = (
            fullpath(event_id, '*BH.[zrt]'),
            fullpath(event_id, 'weights.dat'),
            "/home/rmodrak/data/FK/MDJ2",
            )

        #from mtbench.Alvizuri2018 import to_mij
        from mtpar import tt2cmt
        to_mij = lambda kappa, theta, sigma, M0, gamma, delta : tt2cmt(gamma, delta, M0, kappa, theta, sigma)

        grid = UnstructuredGrid(
            source.items(),
            callback=to_mij)

        process_bw, process_sw = data_processing_handles_FK(
            path_greens, path_weights,
            )

        misfit_bw, misfit_sw = misfit_handles(
            )

        grid_search(
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

