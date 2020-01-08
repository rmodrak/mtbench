#!/usr/bin/env python

import numpy as np
from mtbench import benchmark, progress
from mtbench.AlvizuriTape2018 import names, depths, magnitudes, selected,    data_processing_handles, misfit_handles, fullpath
from mtuq.grid import FullMomentTensorGridRegular


if __name__=='__main__':
    _i, _n = 1, len(selected)

    #
    # loop over selected events from AlvizuriTape2018
    #

    for index in selected:
        progress(_i, _n)

        event_id = names[index]
        depth = depths[index]
        magnitude = magnitudes[index]

        model = "ak135f_mdj2"
        solver = "AxiSEM"

        path_data, path_weights, path_greens = ( 
            fullpath(event_id, '*BH.[zrt]'), 
            fullpath(event_id, 'weights.dat'),
            "/home/rmodrak/data/axisem/ak135f_mdj2-2s",
            )

        sources = FullMomentTensorGridRegular(
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

