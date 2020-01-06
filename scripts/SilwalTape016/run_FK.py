#!/usr/bin/env python

import numpy as np

from mtbench import benchmark
from mtbench.SilwalTape2016 import names, depths, magnitudes, selected,\
    data_processing_handles_FK, misfit_handles, fullpath

from mtuq.grid import DoubleCoupleGridRegular


if __name__=='__main__':
    count = 0
    for _i in selected:
        count += 1
        print('\nEVENT %d of %d\n' % (count, len(selected)))

        event_name = names[_i]

        path_data = fullpath(event_name, '*.[zrt]')
        path_weights = fullpath(event_name, 'weights.dat')
        path_greens = '/store/wf/FK_synthetics/scak'

        solver = 'FK'
        model = 'scak'

        sources = DoubleCoupleGridRegular(
            npts_per_axis=50,
            magnitude=magnitudes[_i])

        process_bw, process_sw = data_processing_handles_FK(
            path_greens, path_weights)

        misfit_bw, misfit_sw = misfit_handles()

        benchmark(
            event_name, path_data, path_greens, path_weights, solver, model,
            process_bw, process_sw, misfit_bw, misfit_sw, 
            sources, magnitudes[_i], depths[_i])

