import os

import datablock


def create_input(geometry, plasma_background, params):
    files = input_file_names('..', params)

    t = {}
    t['geometry'] = datablock.create_geometry(geometry)
    t['background'] = datablock.create_plasma_background(plasma_background)

    t_flx, flx = plasma_background['influx']
    t['influx'] = datablock.create_influx_datablock(t_flx, flx)

    rho = plasma_background['rho']
    D = plasma_background['D']
    v = plasma_background['v']
    params['transport_datablock'] = datablock.create_transport_datablock(rho, D, v)

    t['main'] = datablock.create_param_file(params)

    for key, value in t.iteritems():
        open(files[key], 'w').write(value)


def input_directory_names(maindir):
    """
    >>> o = input_directory_names('xxx')
    >>> print o['nete']
    xxx/nete
    >>> print o['param_files']
    xxx/param_files
    >>> print o['results']
    xxx/results
    """
    out = {}
    out['nete'] = os.path.join(maindir, 'nete')
    out['param_files'] = os.path.join(maindir, 'param_files')
    out['results'] = os.path.join(maindir, 'results')
    return out


def input_file_names(maindir, params):
    """
    >>> import defaults
    >>> o = input_file_names('xxx', defaults.main)

    >>> print o['geometry']
    xxx/nete/grid_99999.0
    >>> print o['background']
    xxx/nete/pp99999.0
    >>> print o['main']
    xxx/param_files/main_99999.0
    >>> print o['influx']
    xxx/nete/Arflx99999.dat
    """
    dirnames = input_directory_names(maindir)
    casename = '%(shot)s.%(index)d' % params

    out = {}
    out['geometry'] = os.path.join(dirnames['nete'], 'grid_%s' % casename)
    out['background'] = os.path.join(dirnames['nete'], 'pp%s' % casename)
    out['main'] = os.path.join(dirnames['param_files'], 'main_%s' % casename)
    out['influx'] = os.path.join(dirnames['nete'],
            'Arflx%(shot)s.dat' % params)

    return out


import numpy as np
def rectangular_pulse(length, max_value, time_offset=0.0):
    """
    >>> x,y = rectangular_pulse(1.0, 1.0)
    >>> print x
    [ -1.00000000e-03   0.00000000e+00   1.00000000e+00   1.00100000e+00]
    >>> print y
    [ 0.  1.  1.  0.]
    """
    delta_change = 1e-3 * length
    t = np.array([-delta_change, 0.0, length, length + delta_change])
    t += time_offset
    y = np.array([0.0, max_value, max_value, 0.0])

    return (t,y)

if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=False)

