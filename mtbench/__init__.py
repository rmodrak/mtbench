#!/usr/bin/env python

import os
import numpy as np

from mtuq import read, open_db, download_greens_tensors
from mtuq.graphics import plot_data_greens, plot_beachball
from mtuq.grid_search import grid_search
from mtuq.util.cap import parse_station_codes, Trapezoid


def progress(_i, _n):
    print('\nEVENT %d of %d\n' % (_i, _n))


def benchmark(
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

    if process_bw:
        data_bw = data.map(process_bw)
    if process_sw:
        data_sw = data.map(process_sw)


    print('Reading Green''s functions...\n')

    if solver.lower() == 'syngine':
        greens = download_greens_tensors(stations, origin, model)
    else:
        db = open_db(path_greens, format=solver)
        greens = db.get_greens_tensors(stations, origin, model)

    greens.convolve(Trapezoid(magnitude=magnitude))
    if process_bw:
        greens_bw = greens.map(process_bw)
    if process_sw:
        greens_sw = greens.map(process_sw)


    #
    # The main computational work starts nows
    #

    if misfit_bw:
        print('Evaluating body wave misfit...\n')

        results_bw = grid_search(
            data_bw, greens_bw, misfit_bw, origin, sources)

    else:
        results_bw = np.zeros(len(sources))


    if misfit_sw:
        print('Evaluating surface wave misfit...\n')

        results_sw = grid_search(
            data_sw, greens_sw, misfit_sw, origin, sources)

    else:
        results_sw = np.zeros(len(sources))


    best_misfit = (results_bw + results_sw).min()
    best_source = sources.get((results_bw + results_sw).argmin())


    #
    # Saving results
    #

    print('Saving results...\n')

    plot_data_greens(event_id+'.png', 
        data_bw, data_sw, greens_bw, greens_sw, process_bw, process_sw, 
        misfit_bw, misfit_sw, stations, origin, best_source)

    plot_beachball(event_id+'_beachball.png', best_source)

    print('Finished\n')


