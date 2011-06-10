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


def velocity_from_zero_flux(rho_vol, rho_pol, D, ne):
    """
    >>> rho_pol = np.linspace(0,1,100)
    >>> rho_vol = 88 * rho_pol # rho_vol for circular cross section
    >>> D = modified_gauss(6, 2, 1.9, 0.4, 0.05, 0.8)(rho_pol)
    >>> ne = 1e20 * (1 - rho_pol**2) + 1e18
    >>> v = velocity_from_zero_flux(rho_vol, rho_pol, D, ne)

    Illustration
    ------------
    >>> import matplotlib.pyplot as plt
    >>> plt.figure(); plt.clf()
    >>> plt.subplot(211)
    >>> plt.plot(rho_pol, D,'.-')
    >>> plt.subplot(212)
    >>> plt.plot(rho_pol, v, '.-')
    >>> plt.show()
    """
    dne = np.gradient(ne) / np.gradient(rho_vol)
    v = D * dne/ne

    outside_mask = rho_pol >= 0.9

    x = rho_vol[outside_mask] - rho_vol[outside_mask][0]
    x /= 3
    vmax = v[outside_mask][0]
    v[outside_mask] = vmax * np.exp(-x**2)
    return v


if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=False)
