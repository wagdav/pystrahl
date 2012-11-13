import numpy as np

import input_templates


def plasma_background(rc):
    rho = rc['background.rhopol']
    ne = rc['background.electron_density']
    te = rc['background.electron_temperature']
    decay_length= rc['background.decay_length']

    out = [ 'DENSITY',
            plasma_profile(rho, ne, decay_length),

            '\nELECTRON TEMPERATURE',
            plasma_profile(rho, te, decay_length),

            '\nION TEMPERATURE',
            plasma_profile(rho, te, decay_length)]

    return ''.join(out)


def geometry(rc):
    geom = {}
    rhovol = rc['geometry.rhovol']
    geom['R_lfs'] = np.zeros_like(rhovol)
    geom['R_hfs'] = np.zeros_like(rhovol)
    geom['geometry.rhovol'] = rhovol
    geom['geometry.rhopol'] = rc['geometry.rhopol']

    for key in geom.keys():
        geom[key] = array2text(geom[key])

    geom['n_grid'] = len(rhovol)
    geom['n_sep'] = np.where(rhovol <= 1.0, 1, 0).sum()
    geom['geometry.rvol_lcfs'] = rc['geometry.rvol_lcfs']
    geom['geometry.major_radius'] = rc['geometry.major_radius']

    out = input_templates.geometry % geom

    dummy = array2text(np.zeros(geom['n_sep']))
    dummy = '\ncv\n%s\n' % dummy
    out += 12*dummy

    return out


def geometry2(rc):
    rhopol =  rc['background.rhopol']

    r0 = rc['geometry.major_radius'] / 100.
    return strahlgrid_circular(rhopol, r0)


import subprocess
import inspect, os
def strahlgrid_circular(rho, r0=0.88, z0=0.23, uloop=0., bt0=1.44):
    """
    Writes circular geometry parameters to 'starhlgrid.in' then calls
    the fortran code strahlgrid.

    strahlgrid.in => [strahlgrid] => output_file
    """
    NANGLE = 103
    shot = 99999
    time = 0
    #r0, z0, uloop = 0.88, 0.23, 0.
    #bt0 = 1.44

    nr = np.alen(rho)
    nsep = nr

    a = 0.22 * rho
    rmaj_in = r0 - a
    rmaj_out = r0 + a

    theta = np.linspace(0, 2*np.pi, NANGLE)
    amin = 0.22 * np.linspace(0, 1, nsep)
    theta, amin = np.meshgrid(theta,amin)
    r = r0 + amin * np.cos(theta)
    z = z0 + amin * np.sin(theta)
    eps = amin / r0

    q0, qedge = 1, 6.
    q = (qedge - q0) * rho[:,None]**2 + q0
    bt = bt0 / (1 + eps * np.cos(theta)) # toroidal field

    bp0 = bt0 * eps / (q * np.sqrt(1-eps**2))
    bp = bp0 / (1 + eps * np.cos(theta)) # poloidal field

    def a2s(a):
        txt = ' '.join('%1.3f' % i for i in a)
        txt += '\n'
        return txt

    print 'Writing strahlgrid.in'
    f = open('strahlgrid.in', 'w')
    f.write('%s\n' % shot)
    f.write('%s\n' % time)
    f.write('%s %s %s\n' % (r0, z0, uloop))
    f.write('%s %s\n' % (nr, nsep))
    f.write(a2s(rho))
    f.write(a2s(rmaj_in))
    f.write(a2s(rmaj_out))
    f.write(a2s(r.ravel('F')))
    f.write(a2s(z.ravel('F')))
    f.write(a2s(bp.ravel('F')))
    f.write(a2s(bt.ravel('F')))
    f.close()

    curdir = os.path.dirname(inspect.getfile(inspect.currentframe()))
    dir_ = os.path.join(curdir, '..')
    strahlgrid = os.path.join(dir_, 'strahlgrid')
    if not os.path.isfile(strahlgrid):
        raise IOError('strahlgrid is not found in %s '
            'Have you compiled it?' % dir_)

    p = subprocess.Popen(strahlgrid, stdout=subprocess.PIPE)
    output, error = p.communicate()

    return output


def main_parameter_file(rc):
    p = {}
    p['transport_datablock'] = transport_properties(rc)
    p.update(rc)

    return input_templates.param_file % p


def plasma_profile(x, y, decay_length):
    """
    >>> o = plasma_profile([0, 1, 2],[1,2,3], 1.0)
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
    cv decay length in rvol [cm]
      1.00
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

    return input_templates.pp_datablock % dd


def transport_properties(rc):
    """
    >>> r = [0, 1, 2, 3, 4, 5]
    >>> D = [1, 2, 3, 4, 5, 6]
    >>> v = [1, 2, 3, 4, 5, 6]
    >>> rc = {  'impurity.rhopol' : r,
    ...         'impurity.diffusion_coefficient' : D,
    ...         'impurity.convection_velocity' : v}
    >>> o = transport_properties(rc)
    """
    r = rc['impurity.rhopol']
    D = rc['impurity.diffusion_coefficient']
    v = rc['impurity.convection_velocity']

    n_points = len(r)
    x = array2text(r)
    D = array2text(D)
    v = array2text(v)

    o = ''
    o += input_templates.transport_datablock % dict(n_points=n_points, x=x, y=D,
            drift='')
    o += input_templates.transport_datablock % dict(n_points=n_points, x=x, y=v,
            drift="'velocity'")

    #o += "cv\n const_c 'parameter'"

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
