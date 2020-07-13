#!/usr/bin/env python

import os
import numpy as np
import warnings

from pandas import DataFrame
from xarray import DataArray

from mtuq import read, open_db, download_greens_tensors
from mtuq.graphics import plot_data_greens
from mtuq.grid_search import grid_search
from mtuq.util.cap import parse_station_codes, Trapezoid
from mtuq.util.math import list_intersect_with_indices
from mtuq.util.signal import get_components



def bench(
    event_id,
    path_data,
    path_greens,
    path_weights,
    solver,
    model,
    process_data_functions,
    misfit_functions,
    grid,
    magnitude,
    depth,
    include_mt=True,
    include_force=False,
    estimate_sigma=False,
    plot_waveforms=True,
    verbose=True):

    nn = len(process_data_functions)

    if len(process_data_functions)!=len(misfit_functions):
        raise Exception()


    if verbose:
        print('event:   %s' % event_id)
        print('data:    %s' % path_data)
        print('weights: %s\n' % path_weights)


    #
    # The main I/O work starts now
    #

    print('Reading data...\n')

    data = read(path_data, format='sac',
        event_id=event_id,
        station_id_list=parse_station_codes(path_weights),
        tags=['units:cm', 'type:velocity']) 

    data.sort_by_distance()
    stations = data.get_stations()

    origin = data.get_origins()[0]
    origin.depth_in_m = depth

    processed_data = []
    for process_data in process_data_functions:
        processed_data += [data.map(process_data)]


    print('Reading Green''s functions...\n')

    db = open_db(path_greens, format=solver,
        model=model, include_mt=include_mt, include_force=include_force)

    greens = db.get_greens_tensors(stations, origin, model)

    greens.convolve(Trapezoid(magnitude=magnitude))

    processed_greens = []
    for process_data in process_data_functions:
        processed_greens += [greens.map(process_data)]


    #
    # The main computational work starts nows
    #

    print('Evaluating misfit...\n')

    results = []
    for _i, misfit in enumerate(misfit_functions):
        task(_i, nn)
        results += [grid_search(processed_data[_i], processed_greens[_i], misfit, origin, grid)]

    # what index corresponds to minimum misfit?
    results_sum = sum(results)

    idx = results_sum.source_idxmin()

    best_source = grid.get(idx)
    source_dict = grid.get_dict(idx)


    if estimate_sigma:
        print('  estimating data variance...\n')

        for _i, misfit in enumerate(misfit_functions):

            groups = misfit.time_shift_groups
            if len(groups) > 1:
               print('Too many time shift groups. Skipping...')
               continue

            components = []
            for component in groups[0]:
               components += [component]

            sigma = _estimate_sigma(processed_data[_i], processed_greens[_i],
                best_source, misfit.norm, components,
                misfit.time_shift_min, misfit.time_shift_max)

            _write(event_id+'_'+str(_i)+'.sigma', sigma)


    #
    # Saving results
    #

    print('Saving results...\n')

    if plot_waveforms:
        print('  plotting waveforms...\n')

        plot_data_greens(event_id+'.png', 
            processed_data[0], processed_data[1], 
            processed_greens[0], processed_greens[1], 
            process_data_functions[0], process_data_functions[1], 
            misfit_functions[0], misfit_functions[1], 
            stations, origin, best_source, source_dict)


    for _i, ds in enumerate(results):
        task(_i, nn)
        _save(event_id+'_'+str(_i), ds)


    print('\nFinished\n')



#
# variance estimation
#

def _estimate_sigma(data, greens, best_source, norm, components,
    time_shift_min, time_shift_max):
    """ A posteriori standard deviation estimate
    """

    # error checking
    assert norm in ('L1', 'L2')

    residuals = []
    for _j, d in enumerate(data):
        _components, indices = list_intersect_with_indices(
            components, get_components(d))

        if not indices:
            continue

        # generate synthetics
        greens[_j]._set_components(_components)
        s = greens[_j].get_synthetics(best_source)

        # time sampling scheme
        npts = d[0].data.size
        dt = d[0].stats.delta

        padding_left = int(+time_shift_max/dt)
        padding_right = int(-time_shift_min/dt)
        npts_padding = padding_left + padding_right

        # array to hold cross correlations
        corr = np.zeros(npts_padding+1)

        #
        # calculate residuals
        # 

        corr[:] = 0.
        for _k in indices:
            corr += np.correlate(s[_k].data, d[_k].data, 'valid')

        npts_shift = padding_left - corr.argmax()
        time_shift = npts_shift*dt

        # what start and stop indices will correctly shift synthetics
        # relative to data?
        start = padding_left - npts_shift
        stop = start + npts

        for _k in indices:

            # substract data from shifted synthetics
            r = s[_k].data[start:stop] - d[_k].data

            # sum the resulting residuals
            if norm=='L1':
                residuals += [np.sum(np.abs(r))*dt]

            elif norm=='L2':
                residuals += [np.sum(r**2)*dt]


    return np.mean(residuals)**0.5



#
# utility functions
#

def _save(filename, results):
    if issubclass(type(results), DataArray):
        results.save(filename+'.nc')

    elif issubclass(type(results), DataFrame):
        results.save(filename+'.h5')


def _write(filename, value):
    np.savetxt(filename, np.array([value]))


def progress(_i, _n):
    print('\nEVENT %d of %d\n' % (_i, _n))


def task(_i, _n):
    if _n > 1:
        print('  task %d of %d' % (_i+1, _n))


def _override(misfit_functions):
    # overrides whatever norms are specified 
    # (because estimates_sigma only allows L2 norm)
    for misfit in misfit_functions:
        if misfit.norm!='L2':
            warnings.warn('Forcing L2 norm')
            misfit.norm='L2'
    return misfit_functions

