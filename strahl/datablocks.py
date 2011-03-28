import numpy as np

import templates


def plasma_background(rc):
    """
    >>> from parameters import defaultParams
    >>> o = plasma_background(defaultParams())
    """
    rho = rc['background.rho_poloidal']
    ne = rc['background.electron_density']
    te = rc['background.electron_temperature']

    out = [ 'DENSITY',
            plasma_profile(rho, ne),

            '\nELECTRON TEMPERATURE',
            plasma_profile(rho, te),

            '\nION TEMPERATURE',
            plasma_profile(rho, te)]

    return ''.join(out)


def geometry(rc):
    """
    >>> from parameters import defaultParams
    >>> o = geometry(defaultParams())
    """
    geom = {}

    radial_grid = rc['geometry.rho_volume']
    geom.update(rc)
    geom['n_grid'] = len(radial_grid)
    geom['n_sep'] = np.where(radial_grid <= 1.0, 1, 0).sum()
    geom['R_lfs'] = np.zeros_like(radial_grid)
    geom['R_hfs'] = np.zeros_like(radial_grid)

    for key in ['background.rho_poloidal', 'geometry.rho_volume', 'R_lfs',
            'R_hfs']:
        geom[key] = array2text(geom[key])

    out = templates.geometry % geom

    dummy = array2text(np.zeros(geom['n_sep']))
    dummy = '\ncv\n%s\n' % dummy
    out += 12*dummy

    return out


def main_parameter_file(rc):
    """
    >>> from parameters import defaultParams
    >>> o = main_parameter_file(defaultParams())
    """
    p = {}
    p['transport_datablock'] = transport_properties(rc)
    p.update(rc)

    return templates.param_file % p


def plasma_profile(x, y, decay_length=1.0):
    """
    >>> o = plasma_profile([0, 1, 2],[1,2,3])
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


def transport_properties(rc):
    """
    >>> r = [0, 1, 2, 3, 4, 5]
    >>> D = [1, 2, 3, 4, 5, 6]
    >>> v = [1, 2, 3, 4, 5, 6]
    >>> rc = {  'background.rho_poloidal' : r,
    ...         'impurity.diffusion_coefficient' : D,
    ...         'impurity.convection_velocity' : v}
    >>> o = transport_properties(rc)
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


def impurity_influx(rc):
    """
    >>> t = [0.5, 1.0]
    >>> flx = [1.25e23, 2.5e23]
    >>> p = {'impurity.influx' : (t, flx)}
    >>> o = impurity_influx(p)
    >>> print o
    2
    0.500 1.250e+23
    1.000 2.500e+23
    <BLANKLINE>
    """
    t, flx = np.asarray(rc['impurity.influx'])
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



if __name__ == '__main__':
    import doctest
    doctest.testmod()
