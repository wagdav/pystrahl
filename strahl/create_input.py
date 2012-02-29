import os

import datablocks
import parameters
from viz import read_results

def run(rc, working_directory):
    create_input(rc, working_directory)

    _go_to_wk_and_run(rc, working_directory)
    of = output_filename(rc, working_directory)
    return read_results(of)


def _go_to_wk_and_run(rc, working_directory):
    curdir = os.getcwd()
    os.chdir(working_directory)
    os.system(_cmd_write_run_file(rc))
    os.system(_cmd_run_case(rc))
    os.chdir(curdir)


def _cmd_write_run_file(rc):
    r"""
    >>> input_dict = {'shot':12345, 'index':0, 'numerical.time.final': 1.23}
    >>> _cmd_write_run_file(input_dict)
    'echo -e "main_12345.0\n1.230\nE" > param_files/run_main_12345.0'
    """
    timeFinal = rc['numerical.time.final']
    mainName = _main(rc)
    cmd = 'echo -e "%s\n%1.3f\nE" > param_files/run_%s' % (mainName,
            timeFinal, mainName)
    return cmd


def _cmd_run_case(rc):
    """
    >>> input_dict = {'shot':12345, 'index':0, 'numerical.time.final': 1.23}
    >>> _cmd_run_case(input_dict)
    'strahl < param_files/run_main_12345.0'
    """
    return './strahl < param_files/run_%s' % _main(rc)

def output_filename(rc, working_directory):
    """
    >>> output_filename({'impurity.element':'Ar'}, '/home/user/strahl_wk')
    '/home/user/strahl_wk/result/Arstrahl_result.dat'
    >>> output_filename({'impurity.element':'Ne'}, '/home/user/strahl_wk')
    '/home/user/strahl_wk/result/Nestrahl_result.dat'
    """
    basename = rc['impurity.element'] + 'strahl_result.dat'
    return os.path.join(working_directory, 'result', basename)


def create_input(rc, working_directory):
    parameters.sanity_check(rc)
    files = input_file_names(working_directory, rc)

    t = {}
    t['geometry'] = datablocks.geometry2(rc)
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
    casename = _casename(params)

    out = {}
    out['geometry'] = os.path.join(dirnames['nete'], 'grid_%s' % casename)
    out['background'] = os.path.join(dirnames['nete'], 'pp%s' % casename)
    out['main'] = os.path.join(dirnames['param_files'], _main(params))
    out['influx'] = os.path.join(dirnames['nete'],
            '%(impurity.element)sflx%(shot)s.dat' % params)

    return out


def _main(rc):
    """
    >>> _main(dict(shot=12345, index=0))
    'main_12345.0'
    """
    return 'main_%s' % _casename(rc)


def _casename(rc):
    """
    >>> _casename(dict(shot=12345, index=0))
    '12345.0'
    """
    return '%(shot)s.%(index)d' % rc


if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=False)

