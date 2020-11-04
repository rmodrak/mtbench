#!/usr/bin/env python

#
# This scripts generates the contents of  mtbench/scripts
#

Imports="""#!/usr/bin/env python

import numpy as np
from mtbench import bench, progress
from _REFERENCE import fullpath, names, depths, magnitudes,\\
    data_processing, misfit_functions, selected_events, expected_results
from mtuq.grid import DoubleCoupleGridRegular

"""


Docstring="""
if __name__=='__main__':
    _i, _n = 1, len(selected_events)

    #
    # loop over selected events from REFERENCE
    #
"""


Main="""
    for index in selected_events:
        progress(_i, _n)

        event_id = names[index]
        depth = depths[index]
        magnitude = magnitudes[index]

        model = MODEL
        solver = SOLVER

        path_data, path_weights, path_greens = ( 
            fullpath(event_id, '*BH.[zrt]'), 
            fullpath(event_id, 'weights.dat'),
            PATH_GREENS,
            )

        grid = DoubleCoupleGridRegular(
            npts_per_axis=15,
            magnitudes=[magnitude],
            )

        process_bw, process_sw = data_processing(
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
            )

        _i += 1

"""


if __name__=='__main__':
    import re

    #
    # Silwal2016
    #

    with open('scripts/run_Silwal2016_AxiSEM.py', 'w') as file:
        lines = ''
        lines += Imports
        lines += Docstring
        lines += Main

        lines = re.sub(
            'REFERENCE',
            'Silwal2016',
            lines)

        lines = re.sub(
            'MODEL',
            '"scak_ak135f"',
            lines)

        lines = re.sub(
            'SOLVER',
            '"AxiSEM"',
            lines)

        lines = re.sub(
            'PATH_GREENS',
            '"/home/rmodrak/data/axisem/scak_ak135f-2s"',
            lines)

        lines = re.sub(
            'include_bw=False',
            'include_bw=True',
            lines)

        file.write(lines)


    with open('scripts/run_Silwal2016_FK.py', 'w') as file:
        lines = ''
        lines += Imports
        lines += Docstring
        lines += Main

        lines = re.sub(
            'REFERENCE',
            'Silwal2016',
            lines)

        lines = re.sub(
            'MODEL',
            '"scak"',
            lines)

        lines = re.sub(
            'PATH_GREENS',
            '"/store/wf/FK_synthetics/scak"',
            lines)

        lines = re.sub(
            'SOLVER',
            '"FK"',
            lines)

        lines = re.sub(
            'data_processing',
            'data_processing_FK',
            lines)

        lines = re.sub(
            'include_bw=False',
            'include_bw=True',
            lines)

        file.write(lines)


    with open('scripts/run_Silwal2016_syngine.py', 'w') as file:
        lines = ''
        lines += Imports
        lines += Docstring
        lines += Main

        lines = re.sub(
            'REFERENCE',
            'Silwal2016',
            lines)

        lines = re.sub(
            'MODEL',
            '"ak135"',
            lines)

        lines = re.sub(
            'SOLVER',
            '"syngine"',
            lines)

        lines = re.sub(
            'PATH_GREENS',
            '"http://service.iris.edu/irisws/syngine/1"',
            lines)

        lines = re.sub(
            'include_bw=False',
            'include_bw=True',
            lines)

        file.write(lines)


    #
    # Alvizuri2018
    #

    with open('scripts/run_Alvizuri2018_AxiSEM.py', 'w') as file:
        lines = ''
        lines += Imports
        lines += Docstring
        lines += Main

        lines = re.sub(
            'REFERENCE',
            'Alvizuri2018',
            lines)

        lines = re.sub(
            'MODEL',
            '"mdj2_ak135f_celso"',
            lines)

        lines = re.sub(
            'SOLVER',
            '"AxiSEM"',
            lines)

        lines = re.sub(
            'PATH_GREENS',
            '"/home/rmodrak/data/axisem/mdj2_ak135f_celso-2s"',
            lines)

        lines = re.sub(
            'DoubleCoupleGridRegular',
            'FullMomentTensorGridRandom',
            lines)

        lines = re.sub(
            'npts_per_axis=15',
            'npts=2000000',
            lines)

        file.write(lines)


    with open('scripts/run_Alvizuri2018_FK.py', 'w') as file:
        lines = ''
        lines += Imports
        lines += Docstring
        lines += Main

        lines = re.sub(
            'REFERENCE',
            'Alvizuri2018',
            lines)

        lines = re.sub(
            'MODEL',
            '"MDJ2"',
            lines)

        lines = re.sub(
            'PATH_GREENS',
            '"/home/rmodrak/data/FK/MDJ2"',
            lines)

        lines = re.sub(
            'SOLVER',
            '"FK"',
            lines)

        lines = re.sub(
            'data_processing',
            'data_processing_FK',
            lines)

        lines = re.sub(
            'DoubleCoupleGridRegular',
            'FullMomentTensorGridRandom',
            lines)

        lines = re.sub(
            'npts_per_axis=15',
            'npts=1000000',
            lines)

        file.write(lines)


    with open('scripts/run_Alvizuri2018_syngine.py', 'w') as file:
        lines = ''
        lines += Imports
        lines += Docstring
        lines += Main

        lines = re.sub(
            'REFERENCE',
            'Alvizuri2018',
            lines)

        lines = re.sub(
            'MODEL',
            '"ak135"',
            lines)

        lines = re.sub(
            'SOLVER',
            '"syngine"',
            lines)

        lines = re.sub(
            'PATH_GREENS',
            '"http://service.iris.edu/irisws/syngine/1"',
            lines)

        lines = re.sub(
            'DoubleCoupleGridRegular',
            'FullMomentTensorGridRandom',
            lines)

        lines = re.sub(
            'npts_per_axis=15',
            'npts=1000000',
            lines)

        file.write(lines)

