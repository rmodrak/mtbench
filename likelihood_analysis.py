#!/usr/bin/env python

import os
import numpy as np
import subprocess

from mtuq import read, open_db, download_greens_tensors
from mtuq.graphics import plot_data_greens, plot_beachball
from mtuq.graphics.uq import plot_likelihood, plot_marginal, plot_misfit
from mtuq.grid_search import grid_search
from mtuq.misfit.level0 import estimate_sigma
from mtuq.util.cap import parse_station_codes, Trapezoid
from mtuq.util.math import list_intersect_with_indices
from mtuq.util.signal import get_components


titles = [
    'Body wave',
    'Rayleigh',
    'Love',
    ]


def progress(_i, _n):
    print('\nEVENT %d of %d\n' % (_i, _n))


def likelihood_analysis(
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
    _i = -1

    for process_data, misfit_function in zip(
        process_data_handles, misfit_function_handles):

        _i += 1

        processed_data += [data.map(process_data)]
        processed_greens += [greens.map(process_data)]

        misfit_values += [grid_search(processed_data[-1], processed_greens[-1],
            misfit_function, origin, grid)]

        struct = grid.to_dataarray(misfit_values[_i])
        title = titles[_i]+' misfit'
        plot_misfit(event_id+'_misfit_'+str(_i)+'.png', struct, title=title)



    #
    # What is the best source?
    #

    misfit_sum = np.zeros((grid.size, 1))
    for _i in range(N):
        misfit_sum += misfit_values[_i]

    struct = grid.to_dataarray(misfit_sum)
    title = 'Total misfit'
    plot_misfit(event_id+'_misfit_sum.png', struct, title=title)

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
            best_source, 'L2', components[0], time_shift_min, time_shift_max)]


    normalized_values = []
    for _i in range(N):
        if np.isnan(sigma[_i]):
            normalized_values += [[]]
            continue

        normalized_values += [misfit_values[_i].copy()/sigma[_i]]

        struct = grid.to_dataarray(normalized_values[_i])

        title = titles[_i]+' maximum likelihood'+u'\n@~s@~ = %.1e' % sigma[_i]
        plot_likelihood(event_id+'_likelihood_'+str(_i)+'.png', struct, title=title)

        #title = 'Marginal likelihood\n'+subtitles[_i]
        #plot_marginal(event_id+'_marginal'+str(_i)+'.png', struct, title=title)


    normalized_sum = np.zeros((grid.size, 1))
    for _i in range(N):
        if np.isnan(sigma[_i]):
            continue
        normalized_sum += normalized_values[_i]

    struct = grid.to_dataarray(normalized_sum)

    title = 'Total maximum likelihood'
    plot_likelihood(event_id+'_likelihood_sum.png', struct, title=title)

    #title = 'Marginal likelihood\nTotal'
    #plot_marginal(event_id+'_marginal.png', struct, title=title)

    print('Finished\n')




