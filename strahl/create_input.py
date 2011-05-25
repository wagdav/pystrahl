import os

import datablocks
import parameters


def create_input(rc, working_directory):
    parameters.sanity_check(rc)
    files = input_file_names(working_directory, rc)

    t = {}
    t['geometry'] = datablocks.geometry(rc)
    t['background'] = datablocks.plasma_background(rc)
    t['influx'] = datablocks.impurity_influx(rc)
    t['main'] = datablocks.main_parameter_file(rc)

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
    >>> from parameters import defaultParams
    >>> o = input_file_names('xxx', defaultParams())

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
            '%(impurity.element)sflx%(shot)s.dat' % params)

    return out


if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=False)

