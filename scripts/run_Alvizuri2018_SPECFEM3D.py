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

        #model = "1D_ak135f_no_mud"
        model = "s40rts_crust1.0"

        solver = "SPECFEM3D"

        path_data, path_weights, path_greens = ( 
            fullpath(event_id, '*BH.[zrt]'), 
            fullpath(event_id, 'weight_celso.dat'),
            "/Users/rmodrak/Downloads/greens/output/NKT/"+model,
            )

        grid = FullMomentTensorGridRandom(
            npts=500000,
            magnitudes=np.linspace(magnitude-1, magnitude+1, 5),
            )

        process_bw, process_sw = data_processing(
            path_greens, path_weights,
            )

        print(process_sw.__dict__)

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
            )

        _i += 1

