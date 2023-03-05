#e!/usr/bin/env python

import os
import numpy as np
import warnings
from os.path import join

from mtuq import read, open_db, download_greens_tensors
from mtuq.graphics import plot_data_greens1, plot_data_greens2,\
    plot_misfit_lune, plot_likelihood_lune, plot_marginal_lune,\
    plot_misfit_vw, plot_likelihood_vw, plot_marginal_vw,\
    plot_misfit_dc, plot_likelihood_dc, plot_marginal_dc,\
    plot_variance_reduction_lune, plot_time_shifts, plot_amplitude_ratios
from mtuq.graphics import plot_beachball as _plot_beachball
from mtuq.grid_search import DataArray, DataFrame, grid_search
from mtuq.misfit import Misfit
from mtuq.misfit.waveform._stats import estimate_sigma, calculate_norm_data
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
    calculate_sigma=False,
    calculate_norm_data=False,
    save_misfit=False,
    plot_beachball=True,
    plot_waveforms=True,
    lune_misfit=False,
    lune_likelihood=False,
    lune_marginal=False,
    lune_variance_reduction=False,
    vw_misfit=False,
    vw_likelihood=False,
    vw_marginal=False,
    dc_misfit=False,
    dc_likelihood=False,
    dc_marginal=False,
    verbose=True):

    """ Carries out a separate grid search for each chosen data type and
    performs simple statistical analyses
    """

    #
    # parameter checking
    #
    if any((
        lune_likelihood,
        lune_marginal,
        vw_likelihood,
        vw_marginal,
        dc_likelihood,
        dc_marginal,
        )):
        calculate_sigma = True


    if any((
        lune_variance_reduction,
        )):
        calculate_norm_data = True

    labels = []
    if include_bw:
        labels += ['bw']
    if include_rayleigh:
        labels += ['rayleigh']
    if include_love:
        labels += ['love']

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

    print('SOLVER:', solver)
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

    idx = results_sum.source_idxmin()

    best_source = grid.get(idx)
    source_dict = grid.get_dict(idx)


    if calculate_sigma:
        print('  estimating variance...\n')

        devs = []
        vars = []

        for _i, misfit in enumerate(misfit_functions):

            groups = misfit.time_shift_groups
            if len(groups) > 1:
               print('Too many time shift groups. Skipping...')
               continue

            components = []
            for component in groups[0]:
               components += [component]

            devs += [estimate_sigma(processed_data[_i], processed_greens[_i],
                best_source, misfit.norm, components,
                misfit.time_shift_min, misfit.time_shift_max)]

            vars += [devs[-1]**2]

            _write(event_id+'_'+str(_i)+'.sigma', devs[-1])


    if calculate_norm_data:
        print('  calculating data norm...\n')

        norm_data = []
        for _i, misfit in enumerate(misfit_functions):

            groups = misfit.time_shift_groups
            if len(groups) > 1:
               print('Too many time shift groups. Skipping...')
               continue

            components = []
            for component in groups[0]:
               components += [component]

            norm_data += [calculate_norm_data(processed_data[_i], misfit.norm, components)]

            _write(event_id+'_'+str(_i)+'.norm_data', norm_data[-1])


    #
    # Generating figures
    #
    print('Generating figures...\n')

    if plot_beachball:
        print('  plotting beachball...\n')

        _plot_beachball(event_id+'_beachball.png',
            best_source, stations, origin)

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

    if lune_misfit:
        print('  plotting misfit...')
        _map(event_id+'_misfit_lune', labels, plot_misfit_lune, results)

    if lune_likelihood:
        print('  plotting maximum likelihoods...')
        _map(event_id+'_likelihood_lune', labels, plot_likelihood_lune, results, vars)

    if lune_marginal:
        print('  plotting marginal likelihoods...')
        _map(event_id+'_marginal_lune', labels, plot_marginal_lune, results, vars)

    if lune_variance_reduction:
        pass


    if vw_misfit:
        print('  plotting misfit...')
        _map(event_id+'_misfit_vw', labels, plot_misfit_vw, results)

    if vw_likelihood:
        print('  plotting maximum likelihoods...')
        _map(event_id+'_likelihood_vw', labels, plot_likelihood_vw, results, vars)

    if vw_marginal:
        print('  plotting marginal likelihoods...')
        _map(event_id+'_marginal_vw', labels, plot_marginal_vw, results, vars)


    if dc_misfit:
        print('  plotting misfit...')
        _map(event_id+'_misfit_dc', labels, plot_misfit_dc, results)

    if dc_likelihood:
        print('  plotting maximum likelihoods...')
        _map(event_id+'_likelihood_dc', labels, plot_likelihood_dc, results, vars)

    if dc_marginal:
        print('  plotting maximum likelihoods...')
        _map(event_id+'_marginal_dc', labels, plot_marginal_dc, results, vars)



    #
    # Saving results
    #

    if save_misfit:
        print('Saving results...\n')

        for _i, ds in enumerate(results):
            task(_i, ntasks)
            _save(event_id+'_'+str(_i), ds)


    print('\nFinished\n')



#
# graphics
#

def _map(dirname, labels, func, *sequences, **kwargs):
    """ Used to map plotting function onto a sequence of misfit or likelihood
    surfaces
    """

    os.makedirs(dirname, exist_ok=True)

    for _i, arg_list in enumerate(zip(*sequences)):
        filename = join(dirname, '%s.png' % labels[_i])

        # call plotting function
        func(filename, *arg_list, **kwargs)


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

