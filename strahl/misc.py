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

if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=False)
