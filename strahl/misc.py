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

if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=False)

    length = 1e-3
    t, y = rectangular_pulse_with_decay(length, 1, tau=length)
