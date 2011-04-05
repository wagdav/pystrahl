import numpy as np

from line import line


def chord_endpoints(chord, xx, yy):
    """
    >>> x, y, z = _test_profile()
    >>> chord = [(0,0.5), (0.0, 0.6)]
    >>> chord_endpoints(chord, x, y)
    (0, 0, 2, 5)
    """
    x_coords, y_coords = chord

    x_range = xx[0]
    y_range = yy[:,0]
    ix0 = find_nearest_index(x_range, x_coords[0])
    ix1 = find_nearest_index(x_range, x_coords[1])
    iy0 = find_nearest_index(y_range, y_coords[0])
    iy1 = find_nearest_index(y_range, y_coords[1])

    p = (ix0 + iy0) + (ix1 + iy1)
    return p


def get_values_between_indices(a, p):
    l = line(*p)
    ret = []
    for i in l:
        try:
            value = a[i]
            ret.append(value)
        except IndexError:
            print 'Coord outside the range of the array'

    return np.array(ret)


def profile_along_chord(chord, x, y, f):
    """
    >>> x, y, f = _test_profile()
    >>> chord = [(0,0.5), (0.0, 0.6)]
    >>> d, p = profile_along_chord(chord, x, y, f)
    """
    p = chord_endpoints(chord, x, y)

    xl = get_values_between_indices(x, p)
    yl = get_values_between_indices(y, p)
    f = get_values_between_indices(f, p)

    distance = np.hypot(xl - xl[0], yl - yl[0])
    return distance, f


def integrate_along_chord(chord, x, y, f):
    """
    >>> x, y, f = _test_profile()
    >>> chord = [(0,0.5), (0.0, 0.6)]
    >>> o = integrate_along_chord(chord, x, y, f)
    >>> print round(o, 3)
    0.661
    """
    d, p = profile_along_chord(chord, x, y, f)
    return np.trapz(p, d)


def find_nearest_index(arr, value):
    """
    >>> import numpy as np
    >>> x = np.linspace(0,1,10)
    >>> find_nearest_index(x, 0.5)
    (4,)
    """
    idx=(np.abs(arr - value)).argmin()
    return np.unravel_index(idx, arr.shape)


def _test_profile():
    x, y = np.meshgrid(np.linspace(0,2.0,10), np.linspace(0,1,10))
    z = 1 - x**2

    return x,y,z


def _test_figures():
    import matplotlib.pyplot as plt
    plt.figure()
    x, y, z = _test_profile()
    plt.contour(x, y, z)


if __name__ == '__main__':
    import doctest
    doctest.testmod()


