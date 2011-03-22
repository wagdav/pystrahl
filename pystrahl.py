import numpy as np

import templates


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
      1.0
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
    time_vector = 1.0
    npoints = len(x)

    x = np.asarray(x)
    y = np.asarray(y)

    assert x.shape == y.shape, 'x and y must have the same shape.'
    rgrid=array2text(x)
    ygrid=array2text(y, scale=True)

    dd=dict(ntimes=ntimes, time_vector=time_vector, npoints=npoints,
            rgrid=rgrid, ygrid=ygrid, decay_length=decay_length)

    return templates.pp_datablock % dd


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
