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


def rectangular_pulse_with_decay(length, max_value, tau, t_offset=0.0):
    r"""
    Create a rectangular pulse with a *length* flattop followed by an
    exponential decay that has a characteristic time *tau*

      |-------|
      |         -
      |          \
      |           -
    __|            \_
       <------><---->
        length  tau
    """
    t, y = rectangular_pulse(length, max_value, t_offset)

    t_decay = np.linspace(0, 5 * tau, 10)
    y_decay = np.exp(-t_decay/tau)
    t_decay += t_offset + length
    y_decay *= max_value

    y_decay[-1] = 0
    t = np.hstack((t[:-1], t_decay))
    y = np.hstack((y[:-1], y_decay))

    return t, y


class modified_gauss(object):
    """
    Modified Gauss function based on STRAHL manual p. 35.

    Example
    -------
    >>> f = modified_gauss(0.1, 1.0, 2.5, 0.25, 0.08, 0.8)
    >>> x = np.array([0.0, 0.8, 1.0])
    >>> print np.around(f(x), 2)
    [ 0.1  2.5  1. ]
    """
    def __init__(self, y0, ysol, ymax, beta, gamma, rmax):
        left = lambda x: y0 + (ymax - y0) * np.exp(-(x - rmax)**2/beta**2)
        right = lambda x: ysol + (ymax - ysol) * np.exp(-(x - rmax)**2/gamma**2)

        self.left = left
        self.right = right
        self.rmax = rmax

    def __call__(self, r):
        condlist = [r <= self.rmax, self.rmax < r]
        funclist = [self.left, self.right]
        return np.piecewise(r, condlist, funclist)


def constant(a):
    return np.vectorize(lambda x: a)


def velocity_from_zero_flux(rho, D, ne):
    """
    >>> rho = np.linspace(0,1,20)
    >>> D = rho**2 + 0.02
    >>> ne = 1e20 * (1 -rho**2) + 1e18
    >>> v = velocity_from_zero_flux(rho, D, ne)
    """
    dne = np.gradient(ne)/np.gradient(rho)
    v = 100 * D * dne/ne

    v[rho > 0.9] = 0
    return v


if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=False)
