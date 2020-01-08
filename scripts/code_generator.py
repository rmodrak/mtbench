# This scripts generates all the other scripts found in mtbench/scripts


Imports="""#!/usr/bin/env python

import numpy as np
from mtbench import benchmark, progress
from mtbench.REFERENCE import names, depths, magnitudes, selected,\
    data_processing_handles, misfit_handles, fullpath
from mtuq.grid import DoubleCoupleGridRegular

"""


Docstring="""
if __name__=='__main__':
    _i, _n = 1, len(selected)

    #
    # loop over selected events from REFERENCE
    #
"""


Main="""
    for index in selected:
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

"""


if __name__=='__main__':
    import re

    #
    # SilwalTape2016
    #

    with open('SilwalTape2016/run_AxiSEM.py', 'w') as file:
        lines = ''
        lines += Imports
        lines += Docstring
        lines += Main

        lines = re.sub(
            'REFERENCE',
            'SilwalTape2016',
            lines)

        lines = re.sub(
            'MODEL',
            '"ak135f_scak"',
            lines)

        lines = re.sub(
            'SOLVER',
            '"AxiSEM"',
            lines)

        lines = re.sub(
            'PATH_GREENS',
            '"/home/rmodrak/data/axisem/ak135f_scak-2s"',
            lines)

        file.write(lines)


    with open('SilwalTape2016/run_FK.py', 'w') as file:
        lines = ''
        lines += Imports
        lines += Docstring
        lines += Main

        lines = re.sub(
            'REFERENCE',
            'SilwalTape2016',
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
            'data_processing_handles',
            'data_processing_handles_FK',
            lines)

        file.write(lines)


    with open('SilwalTape2016/run_syngine.py', 'w') as file:
        lines = ''
        lines += Imports
        lines += Docstring
        lines += Main

        lines = re.sub(
            'REFERENCE',
            'SilwalTape2016',
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

        file.write(lines)


    #
    # AlvizuriTape2018
    #

    with open('AlvizuriTape2018/run_AxiSEM.py', 'w') as file:
        lines = ''
        lines += Imports
        lines += Docstring
        lines += Main

        lines = re.sub(
            'REFERENCE',
            'AlvizuriTape2018',
            lines)

        lines = re.sub(
            'MODEL',
            '"ak135f_mdj2"',
            lines)

        lines = re.sub(
            'SOLVER',
            '"AxiSEM"',
            lines)

        lines = re.sub(
            'PATH_GREENS',
            '"/home/rmodrak/data/axisem/ak135f_mdj2-2s"',
            lines)

        lines = re.sub(
            'DoubleCouple',
            'FullMomentTensor',
            lines)

        file.write(lines)


    with open('AlvizuriTape2018/run_FK.py', 'w') as file:
        lines = ''
        lines += Imports
        lines += Docstring
        lines += Main

        lines = re.sub(
            'REFERENCE',
            'AlvizuriTape2018',
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
            'data_processing_handles',
            'data_processing_handles_FK',
            lines)

        lines = re.sub(
            'DoubleCouple',
            'FullMomentTensor',
            lines)

        file.write(lines)


    with open('AlvizuriTape2018/run_syngine.py', 'w') as file:
        lines = ''
        lines += Imports
        lines += Docstring
        lines += Main

        lines = re.sub(
            'REFERENCE',
            'AlvizuriTape2018',
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
            'DoubleCouple',
            'FullMomentTensor',
            lines)

        file.write(lines)


