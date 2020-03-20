#!/usr/bin/env python

import os
import numpy as np

from mtuq import read, open_db, download_greens_tensors
from mtuq.graphics import plot_data_greens, plot_beachball
from mtuq.graphics.uq_vw import plot_likelihood_vw, plot_misfit_vw
from mtuq.grid_search import grid_search
from mtuq.util.cap import parse_station_codes, Trapezoid
from mtuq.util.math import list_intersect_with_indices
from mtuq.util.signal import get_components



def progress(_i, _n):
    print('\nEVENT %d of %d\n' % (_i, _n))


def run_grid_search(
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
    depth,
    verbose=True):


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

    data_bw = data.map(process_bw)
    data_sw = data.map(process_sw)


    print('Reading Green''s functions...\n')
    db = open_db(path_greens, format=solver)
    greens = db.get_greens_tensors(stations, origin, model)

    greens.convolve(Trapezoid(magnitude=magnitude))
    greens_bw = greens.map(process_bw)
    greens_sw = greens.map(process_sw)


    #
    # The main computational work starts nows
    #

    print('Evaluating body wave misfit...\n')
    results_bw = grid_search(data_bw, greens_bw, misfit_bw, origin, grid)

    print('Evaluating surface wave misfit...\n')
    results_sw = grid_search(data_sw, greens_sw, misfit_sw, origin, grid)

    results_sum = results_bw + results_sw

    best_misfit = (results_sum).min()
    best_source = grid.get(results_sum.argmin())
    lune_dict = grid.get_dict(results_sum.argmin())



    #
    # Saving results
    #

    print('Saving results...\n')

    plot_data_greens(event_id+'.png', 
        data_bw, data_sw, greens_bw, greens_sw, process_bw, process_sw, 
        misfit_bw, misfit_sw, stations, origin, best_source, lune_dict)

    plot_beachball(event_id+'_beachball.png', best_source)

    grid.save(event_id+'.nc', results_sum)

    print('Finished\n')



def plot_likelihoods(
    event_id,
    path_data,
    path_greens,
    path_weights,
    solver,
    model,
    process_data_handles,
    misfit_function_handles,
    grid,
    magnitude,
    depth,
    verbose=True):


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

    db = open_db(path_greens, format=solver)
    greens = db.get_greens_tensors(stations, origin, model)
    greens.convolve(Trapezoid(magnitude=magnitude))

    processed_data = []
    processed_greens = []
    misfit_values = []

    N = len(misfit_function_handles)
    _i = 0

    for process_data, misfit_function in zip(
        process_data_handles, misfit_function_handles):

        processed_data += [data.map(process_data)]
        processed_greens += [greens.map(process_data)]

        misfit_values += [grid_search(processed_data[-1], processed_greens[-1],
            misfit_function, origin, grid)]

        grid.save(event_id+'_values%d.nc' % _i, misfit_values[-1])

        _i += 1


    #
    # What is the best source?
    #

    misfit_sum = np.zeros(misfit_values[0].shape)
    for _i in range(N):
        misfit_sum = misfit_values[_i]

    best_misfit = misfit_sum.min()
    best_source = grid.get(misfit_sum.argmin())
    lune_dict = grid.get_dict(misfit_sum.argmin())


    #
    # Likelihood analysis
    #

    sigma = []
    for _i in range(N):
        time_shift_min = misfit_function_handles[_i].time_shift_min
        time_shift_max = misfit_function_handles[_i].time_shift_max
        components = misfit_function_handles[_i].time_shift_groups

        sigma += [estimate_sigma(processed_data[_i], processed_greens[_i],
            best_source, 'L2', time_shift_min, time_shift_max, components)]


    likelihood_values = []
    for _i in range(N):
        likelihood_values += [np.exp(-misfit_values[_i].copy()/(2.*sigma[_i]**2))]


    likelihood_product = np.ones(likelihood_values[0].shape)
    for _i in range(N):
        likelihood_product *= likelihood_values[_i]

    print('Finished\n')




def estimate_sigma(data, greens, source, norm='L2', time_shift_min=0., time_shift_max=0., components=['ZR']):

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
        s = greens[_j].get_synthetics(source)

        # time sampling scheme
        npts = d[0].data.size
        dt = d[0].stats.delta

        padding_left = int(+time_shift_max/dt)
        padding_right = int(-time_shift_min/dt)
        npts_padding = padding_left + padding_right

        # array to hold cross correlations
        corr = np.zeros(npts_padding+1)

        #
        # evaluate misfit
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
