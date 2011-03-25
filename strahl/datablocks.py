import numpy as np

import templates


def create_plasma_background(rc):
    """
    >>> import defaults
    >>> o = create_plasma_background(defaults.defaultParams)
    """
    rho = rc['background.rho_poloidal']
    ne = rc['background.electron_density']
    te = rc['background.electron_temperature']

    out = [ 'DENSITY',
            create_pp_datablock(rho, ne),

            '\nELECTRON TEMPERATURE',
            create_pp_datablock(rho, te),

            '\nION TEMPERATURE',
            create_pp_datablock(rho, te)]

    return ''.join(out)


def create_geometry(rc):
    """
    >>> import defaults
    >>> o = create_geometry(defaults.defaultParams)
    """
    geom = {}
    for key in rc.keys():
        if not key.startswith('geometry'): continue
        newkey = key.replace('geometry.', '')
        geom[newkey] = rc[key]

    for key in ['rho_pol', 'rho_vol', 'R_lfs', 'R_hfs']:
        geom[key] = array2text(geom[key])
    out = templates.geometry % geom

    dummy = array2text(np.zeros(geom['n_sep']))
    dummy = '\ncv\n%s\n' % dummy
    out += 12*dummy

    return out


def create_param_file(rc):
    """
    >>> import defaults
    >>> o = create_param_file(defaults.defaultParams)
    """
    p = {}
    p['transport_datablock'] = create_transport_datablock(rc)
    p.update(rc)

    return templates.param_file % p


def array2text(a, scale=False, cols_per_row=6):
    """
    >>> array2text([20, 30, 40])
    '2.0000e+01 3.0000e+01 4.0000e+01'
    >>> array2text([20, 30, 40], scale=True)
    '2.0000e+01 1.0000e+00 1.5000e+00 2.0000e+00'
    """
    format_str = '%1.4e'

    s = []
    if scale:
        s.append(format_str % a[0])

    for count, i in enumerate(a, len(s) + 1):
        if scale:
            oo=i/float(a[0])
        else:
            oo=i
        s.append(format_str % oo)
        if not count % cols_per_row:
            s[-1] += '\n'
    return ' '.join(s)


def create_pp_datablock(x, y, decay_length=1.0):
    """
    >>> o = create_pp_datablock([0, 1, 2],[1,2,3])
    >>> print o
    <BLANKLINE>
    cv time vector
      1
      0.0
    <BLANKLINE>
    cv ne function
      'interp'
    <BLANKLINE>
    cv radial coordinate
      'poloidal rho'
    <BLANKLINE>
    cv # of interpolation points
      3
    <BLANKLINE>
    cv radial grid for interpolation
      0.0000e+00 1.0000e+00 2.0000e+00
    <BLANKLINE>
    cv ne[cm^-3]/ Te[eV] / Ti [eV]
      1.0000e+00 1.0000e+00 2.0000e+00 3.0000e+00
    <BLANKLINE>
    cv decay length[cm] in rho_volume
      1.0
    <BLANKLINE>
    """
    ntimes = 1
    time_vector = 0.0
    npoints = len(x)

    x = np.asarray(x)
    y = np.asarray(y)

    assert x.shape == y.shape, 'x and y must have the same shape.'
    rgrid=array2text(x)
    ygrid=array2text(y, scale=True)

    dd=dict(ntimes=ntimes, time_vector=time_vector, npoints=npoints,
            rgrid=rgrid, ygrid=ygrid, decay_length=decay_length)

    return templates.pp_datablock % dd


def create_transport_datablock(rc):
    """
    >>> r = [0, 1, 2, 3, 4, 5]
    >>> D = [1, 2, 3, 4, 5, 6]
    >>> v = [1, 2, 3, 4, 5, 6]
    >>> rc = {  'background.rho_poloidal' : r,
    ...         'impurity.diffusion_coefficient' : D,
    ...         'impurity.convection_velocity' : v}
    >>> o = create_transport_datablock(rc)
    """
    r = rc['background.rho_poloidal']
    D = rc['impurity.diffusion_coefficient']
    v = rc['impurity.convection_velocity']

    n_points = len(r)
    x = array2text(r)
    D = array2text(D)
    v = array2text(v)

    o = ''
    o += templates.transport_datablock % dict(n_points=n_points, x=x, y=D,
            drift='')
    o += templates.transport_datablock % dict(n_points=n_points, x=x, y=v,
            drift="'velocity'")

    return o


def create_influx_datablock(rc):
    """
    >>> t = [0.5, 1.0]
    >>> flx = [1.25e23, 2.5e23]
    >>> p = {'impurity.influx.time' : t, 'impurity.influx.flux' : flx}
    >>> o = create_influx_datablock(p)
    >>> print o
    2
    0.500 1.250e+23
    1.000 2.500e+23
    <BLANKLINE>
    """
    t = np.asarray(rc['impurity.influx.time'])
    flx = np.asarray(rc['impurity.influx.flux'])
    assert len(t) == len(flx), 'Flux error: shape mismatch'

    o = ''
    o += '%d\n' % len(t)
    for i,j in zip(t, flx):
        o += '%1.3f %1.3e\n' % (i, j)

    return o


def casedir_init (casedir=None):
    import os, errno

    if casedir:
        os.mkdir(casedir)
    else:
        casedir = os.getcwd()

    directories=['nete', 'param_files', 'result']

    for d in directories:
        os.mkdir(os.path.join(casedir, d))


if __name__ == '__main__':
    import doctest
    doctest.testmod()
