#!/usr/bin/env python

import os
import numpy as np
import warnings

from pandas import DataFrame
from xarray import DataArray

from mtuq import read, open_db, download_greens_tensors
from mtuq.graphics import plot_data_greens1, plot_data_greens2, plot_misfit_lune
from mtuq.grid_search import grid_search
from mtuq.misfit import Misfit
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
    grid,
    magnitude,
    depth,
    process_bw,
    process_sw,
    minmax_bw=[-2., +2.],
    minmax_sw=[-5., +5.],
    include_bw=False,
    include_rayleigh=True,
    include_love=True,
    include_mt=True,
    include_force=False,
    write_sigma=False,
    write_norm_data=False,
    plot_waveforms=True,
    verbose=True):

    """ Carries out a separate grid search for each chosen data type and
    peforms simple statistical analyses
    """

    data_processing = []
    if include_bw:
        data_processing += [process_bw]
    if include_rayleigh:
        data_processing += [process_sw]
    if include_love:
        data_processing += [process_sw]

    misfit_functions = []
    if include_bw:
        misfit_functions += [_get_misfit_bw(minmax_bw)]
    if include_rayleigh:
        misfit_functions += [_get_misfit_rayleigh(minmax_sw)]
    if include_love:
        misfit_functions += [_get_misfit_love(minmax_sw)]


    if verbose:
        print('event:   %s' % event_id)
        print('data:    %s' % path_data)
        print('weights: %s\n' % path_weights)

    ntasks = int(include_bw)+\
             int(include_rayleigh)+\
             int(include_love)


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
    for process_data in data_processing:
        processed_data += [data.map(process_data)]


    print('Reading Green''s functions...\n')

    db = open_db(path_greens, format=solver,
        model=model, include_mt=include_mt, include_force=include_force)

    greens = db.get_greens_tensors(stations, origin, model)

    greens.convolve(Trapezoid(magnitude=magnitude))

    processed_greens = []
    for process_data in data_processing:
        processed_greens += [greens.map(process_data)]


    #
    # The main computational work starts nows
    #

    print('Evaluating misfit...\n')

    results = []
    for _i, misfit in enumerate(misfit_functions):
        task(_i, ntasks)
        results += [grid_search(processed_data[_i], processed_greens[_i], 
                                misfit, origin, grid)]

    # what index corresponds to minimum misfit?
    results_sum = sum(results)

    idx = results_sum.idxmin('source')

    best_source = grid.get(idx)
    source_dict = grid.get_dict(idx)


    if write_sigma:
        print('  estimating variance...\n')

        for _i, misfit in enumerate(misfit_functions):

            groups = misfit.time_shift_groups
            if len(groups) > 1:
               print('Too many time shift groups. Skipping...')
               continue

            components = []
            for component in groups[0]:
               components += [component]

            sigma = calculate_sigma(processed_data[_i], processed_greens[_i],
                best_source, misfit.norm, components,
                misfit.time_shift_min, misfit.time_shift_max)

            _write(event_id+'_'+str(_i)+'.sigma', sigma)


    if write_norm_data:
        print('  calculating data norm...\n')

        for _i, misfit in enumerate(misfit_functions):

            groups = misfit.time_shift_groups
            if len(groups) > 1:
               print('Too many time shift groups. Skipping...')
               continue

            components = []
            for component in groups[0]:
               components += [component]

            norm_data = calculate_norm_data(processed_data[_i], misfit.norm, components)

            _write(event_id+'_'+str(_i)+'.norm_data', norm_data)



    #
    # Saving results
    #

    print('Saving results...\n')

    if plot_waveforms:
        print('  plotting waveforms...\n')

        _plot_waveforms(event_id+'_waveforms.png', 
            processed_data,
            processed_greens,
            include_bw,
            bool(include_rayleigh or include_love),
            process_bw,
            process_sw,
            minmax_bw,
            minmax_sw,
            stations,
            origin,
            best_source,
            source_dict)

    for _i, ds in enumerate(results):
        task(_i, ntasks)
        _save(event_id+'_'+str(_i), ds)


    print('\nFinished\n')



def calculate_norm_data(data, norm, components):
    # error checking
    assert norm in ('L1', 'L2')

    norm_data = 0.
    for _j, d in enumerate(data):
        _components, indices = list_intersect_with_indices(
            get_components(d), components)

        if not indices:
            continue

        # time sampling scheme
        npts = d[0].data.size
        dt = d[0].stats.delta

        for _k in indices:
            r = d[_k].data

            if norm=='L1':
                norm_data += np.sum(np.abs(r))*dt

            elif norm=='L2':
                norm_data += np.sum(r**2)*dt

    return norm_data



def calculate_sigma(data, greens, best_source, norm, components,
    time_shift_min, time_shift_max):

    # error checking
    assert norm in ('L1', 'L2')

    residuals = []
    for _j, d in enumerate(data):
        _components, indices = list_intersect_with_indices(
            components, get_components(d))

        if not indices:
            continue

        # generate synthetics
        greens[_j]._set_components(get_components(d))
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


def _plot_waveforms(filename,
    data,
    greens,
    include_bw,
    include_sw,
    process_bw,
    process_sw,
    minmax_bw,
    minmax_sw,
    stations,
    origin,
    best_source,
    source_dict):

    if include_bw and include_sw:
        plot_data_greens2(filename,
            data[0], data[1], 
            greens[0], greens[1],
            process_bw, process_sw, 
            _get_misfit_bw(minmax_bw), _get_misfit_sw(minmax_sw),
            stations, 
            origin, 
            best_source,
            source_dict)

    elif include_sw:
        plot_data_greens1(filename,
            data[0],
            greens[0],
            process_sw,
            _get_misfit_sw(minmax_sw),
            stations, 
            origin, 
            best_source,
            source_dict)


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


def _get_misfit_rayleigh(minmax):
    return Misfit(
        norm='L2',
        time_shift_min=minmax[0],
        time_shift_max=minmax[1],
        time_shift_groups=['ZR'],
        )

def _get_misfit_love(minmax):
    return Misfit(
        norm='L2',
        time_shift_min=minmax[0],
        time_shift_max=minmax[1],
        time_shift_groups=['T'],
        )

def _get_misfit_bw(minmax):
    return Misfit(
        norm='L2',
        time_shift_min=minmax[0],
        time_shift_max=minmax[1],
        time_shift_groups=['ZR'],
        )

def _get_misfit_sw(minmax):
    return Misfit(
        norm='L2',
        time_shift_min=minmax[0],
        time_shift_max=minmax[1],
        time_shift_groups=['ZR','T'],
        )

