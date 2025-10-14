#e!/usr/bin/env python

import os
import numpy as np
import warnings
from os.path import exists, join
from math import prod

from mtuq import read, open_db, download_greens_tensors
from mtuq.graphics import plot_data_greens1, plot_data_greens2,\
    plot_misfit_lune, plot_likelihood_lune, plot_marginal_lune,\
    plot_misfit_vw, plot_likelihood_vw, plot_marginal_vw,\
    plot_misfit_dc, plot_likelihood_dc, plot_marginal_dc,\
    plot_variance_reduction_lune, plot_time_shifts, plot_amplitude_ratios,\
    plot_cdf, plot_pdf, plot_screening_curve,\
    plot_beachball
from mtuq.grid import UnstructuredGrid
from mtuq.grid_search import DataArray, DataFrame, grid_search, MTUQDataArray
from mtuq.misfit import Misfit
from mtuq.misfit.waveform._stats import estimate_sigma, calculate_norm_data
from mtuq.util.cap import parse_station_codes, Trapezoid
from mtuq.util.math import list_intersect_with_indices
from mtuq.util.signal import get_components


# workaround name conflicts
_plot_beachball = plot_beachball
_calculate_norm_data = calculate_norm_data



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
    lune_misfit=True,
    lune_likelihood=False,
    lune_marginal=False,
    lune_variance_reduction=False,
    vw_misfit=False,
    vw_likelihood=False,
    vw_marginal=False,
    dc_misfit=False,
    dc_likelihood=False,
    dc_marginal=False,
    omega_pdfs=False,
    omega_cdfs=False,
    screening_curves=False,
    station_contributions=True,
    path_output='.',
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
        omega_cdfs,
        omega_pdfs,
        screening_curves,
        )):
        calculate_sigma = True


    if any((
        lune_variance_reduction,
        )):
        calculate_norm_data = True

    if omega_pdfs or omega_cdfs:
        try:
            assert type(grid)==UnstructuredGrid
        except:
            print('Angular distance CDFs and PDFs require randomly-spaced grid')
            omega_pdfs = False
            omega_cdfs = False


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


    if not path_output:
        path_output = './'

    if not exists(path_output):
        os.makedirs(path_output, exist_ok=True)

    if verbose:
        print('event:   %s' % event_id)
        print('data:    %s' % path_data)
        print('weights: %s' % path_weights)
        print('output:  %s\n'% path_output)

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

    # holds misfit surfaces from each individual stations and data type
    station_array = []

    # sum over stations to obtain misfit surfaces from each individual data type
    results_sum = []

    for _i, misfit in enumerate(misfit_functions):
        task(_i, ntasks)
        station_array += [[]]

        for _j, station in enumerate(stations):
            print(f'\n  {station.id}\n')

            station_array[-1] += [grid_search(
                processed_data[_i].select(station), processed_greens[_i].select(station), 
                misfit, origin, grid, verbose=0)]

        results_sum += [MTUQDataArray(**{
            'data': np.sum(station_array[-1], axis=0)/len(stations),
            'coords': station_array[-1][0].coords,
            'dims': station_array[-1][0].dims,
            })]

    if calculate_norm_data:
        print('  calculating data norm...\n')
    
        norms = []
        for _i, misfit in enumerate(misfit_functions):
    
            groups = misfit.time_shift_groups
            if len(groups) > 1:
               print('Too many time shift groups. Skipping...')
               continue

            components = []
            for component in groups[0]:
               components += [component]

            norms += [_calculate_norm_data(processed_data[_i], misfit.norm, components)]

            _write(event_id+'_'+str(_i)+'.norm_data', norms[-1])


    if include_rayleigh and include_love:
        idx_rayleigh = labels.index('rayleigh')
        idx_love = labels.index('love')

        results_sum += [results_sum[idx_rayleigh] + results_sum[idx_love]]
        norms += [2.]
        labels += ['rayleigh+love']

        #results_sum += [sum([results_sum[_i]*norms[_i] for _i in range(len(results_sum))])]
        #norms += [norms[idx_rayleigh] + norms[idx_love]]
        #labels += ['rayleigh+love_2']


    # what index corresponds to minimum misfit?
    results_weighted = sum([results_sum[_i]*norms[_i] for _i in range(len(results_sum))])
    idx = results_weighted.source_idxmin()
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


    #
    # Generating figures
    #
    print('Generating figures...\n')

    if plot_beachball:
        print('  plotting beachball...\n')

        _plot_beachball(path_output+'/'+event_id+'_beachball.png',
            best_source, stations, origin)

    if plot_waveforms:
        print('  plotting waveforms...\n')

        _plot_waveforms(path_output+'/'+event_id+'_waveforms.png',
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
        _map(path_output+'/'+event_id+'_misfit_lune', labels, plot_misfit_lune, results_sum)

    if lune_likelihood:
        print('  plotting maximum likelihoods...')
        _map(path_output+'/'+event_id+'_likelihood_lune', labels, plot_likelihood_lune, results_sum, vars)

    if lune_marginal:
        print('  plotting marginal likelihoods...')
        _map(path_output+'/'+event_id+'_marginal_lune', labels, plot_marginal_lune, results_sum, vars)

    if lune_variance_reduction:
        print('  plotting variance reduction...')
        _map(path_output+'/'+event_id+'_variance_reduction', labels, plot_variance_reduction_lune, results_sum, [1. for _ in range(len(results_sum))])


    if vw_misfit:
        print('  plotting misfit...')
        _map(path_output+'/'+event_id+'_misfit_vw', labels, plot_misfit_vw, results_sum)

    if vw_likelihood:
        print('  plotting maximum likelihoods...')
        _map(path_output+'/'+event_id+'_likelihood_vw', labels, plot_likelihood_vw, results_sum, vars)

    if vw_marginal:
        print('  plotting marginal likelihoods...')
        _map(path_output+'/'+event_id+'_marginal_vw', labels, plot_marginal_vw, results_sum, vars)


    if dc_misfit:
        print('  plotting misfit...')
        _map(path_output+'/'+event_id+'_misfit_dc', labels, plot_misfit_dc, results_sum)

    if dc_likelihood:
        print('  plotting maximum likelihoods...')
        _map(path_output+'/'+event_id+'_likelihood_dc', labels, plot_likelihood_dc, results_sum, vars)

    if dc_marginal:
        print('  plotting maximum likelihoods...')
        _map(path_output+'/'+event_id+'_marginal_dc', labels, plot_marginal_dc, results_sum, vars)


    if omega_pdfs:
        print('  plotting angular distance PDFs...')
        _map(path_output+'/'+event_id+'_omega', [label+'_pdf' for label in labels], plot_pdf, results_sum, vars)

    if omega_cdfs:
        print('  plotting angular distance CDFs...')
        _map(path_output+'/'+event_id+'_omega', [label+'_cdf' for label in labels], plot_cdf, results_sum, vars)

    if screening_curves:
        print('  plotting explosion screening curves...')
        _map(path_output+'/'+event_id+'_curves', labels, plot_screening_curve, results_sum, vars)


    if station_contributions:
        os.makedirs(path_output+'/'+event_id+f'_station_contributions',exist_ok=True)

        for _i, station in enumerate(stations):
            plot_variance_reduction_lune(
                path_output+'/'+event_id+f'_station_contributions/rayleigh_{station.id}.png',
                station_array[0][_i],[1.], title=station.id)
    
            plot_variance_reduction_lune(
                path_output+'/'+event_id+f'_station_contributions/love_{station.id}.png',
                station_array[1][_i],[1.], title=station.id)


    #
    # Saving results
    #

    if save_misfit:
        print('Saving results...\n')

        for _i, ds in enumerate(results_sum):
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
        normalize=True,
        verbose=0,
        )

def _get_misfit_love(minmax):
    return Misfit(
        norm='L2',
        time_shift_min=minmax[0],
        time_shift_max=minmax[1],
        time_shift_groups=['T'],
        normalize=True,
        verbose=0,
        )

def _get_misfit_bw(minmax):
    return Misfit(
        norm='L2',
        time_shift_min=minmax[0],
        time_shift_max=minmax[1],
        time_shift_groups=['ZR'],
        normalize=True,
        verbose=0,
        )

def _get_misfit_sw(minmax):
    return Misfit(
        norm='L2',
        time_shift_min=minmax[0],
        time_shift_max=minmax[1],
        time_shift_groups=['ZR','T'],
        normalize=True,
        verbose=0,
        )

