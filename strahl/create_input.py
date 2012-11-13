import os
import shutil

import datablocks
import result

def run(rc, out=None, wk='./wk'):
    """
    rc : dict
        Input parameters for STRAHL.
    out : str
        Name of the output file.
    wk : str
        Working directory of STRAHL.
    """
    create_input(rc, wk)

    _go_to_wk_and_run(rc, wk)
    of = output_filename(rc, wk)
    if out != None:
        shutil.copyfile(of, out)
    return result.load(of)


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
    """
    out = {}
    out['nete'] = os.path.join(maindir, 'nete')
    out['param_files'] = os.path.join(maindir, 'param_files')
    out['results'] = os.path.join(maindir, 'results')
    return out


def input_file_names(maindir, params):
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
