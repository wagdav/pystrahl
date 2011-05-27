import os
import time

import h5py
import numpy as np
import matplotlib.pyplot as plt

def inverted_data(shot):
    """
    Load the cached inverted data for shot.  Return an InversionData object.
    """
    cache_file = _temp_filename(shot)
    gti = h5py.File(cache_file, 'r')

    inverted = gti['res/inverted']
    rho = gti['res/rho'][:,0]
    time = gti['res/time'][0,:]

    inverted = np.array(inverted)
    return InversionData(shot, rho, time, inverted)


def _temp_filename(shot):
    username = os.environ['USER']
    filename = 'gti_%d.h5' % shot

    return os.path.join('/tmp', username, filename)


def example_data():
    rho = np.linspace(0, 1, 41)
    time = np.linspace(0, 0.6, 100)

    xx, yy = np.meshgrid(rho, time)
    zz = (1 - yy) * np.exp(-xx**2)
    return InversionData(-1, rho, time, zz)


class InversionData(object):
    def __init__(self, shot, rho, time, emissivity):
        self.shot = shot
        self.rho = rho
        self.time = time
        self.emissivity = emissivity

    def select_time(self, begin, end):
        time_mask = (begin <= self.time) & (self.time < end)
        new_time = self.time[time_mask]
        new_emissivity = self.emissivity[time_mask]
        return type(self)(self.shot, self.rho, new_time, new_emissivity)

    def plot_contour(self):
        """
        >>> d = example_data()
        >>> d.plot_contour()
        """
        ax = plt.gca()
        ax.contour(self.time, self.rho, self.emissivity.T)
        ax.set_xlabel(r'$t\ \mathrm{[s]}$')
        ax.set_ylabel(r'$\rho_\psi$')

        if plt.isinteractive():
            ax.figure.canvas.draw()

    def plot_profiles(self, n=10):
        """
        >>> d = example_data()
        >>> d.plot_profiles()
        """
        time_points = len(self.time)
        time_mask = slice(None, None, time_points // n)

        ax = plt.gca()
        ax.plot(self.rho, self.emissivity[time_mask,:].T)
        ax.set_xlabel(r'$\rho_\psi$')

        if plt.isinteractive():
            ax.figure.canvas.draw()

    def animate_profiles(self, n=10, loop=1, dt=0.2):
        """
        >>> d = example_data()
        >>> d.animate_profiles(n=2, dt=0)
        """
        ax = plt.gca()

        time_points = len(self.time)
        time_mask = slice(None, None, time_points // n)
        data = self.emissivity[time_mask,:]

        ax.cla()
        line, = ax.plot(self.rho, data[0,:])
        ax.set_xlabel(r'$\rho_\psi$')
        ax.set_ylim(ymax=data.max())
        for i in xrange(loop):
            for d in data:
                line.set_ydata(d)
                ax.figure.canvas.draw()
                time.sleep(dt)

    def plot_time_evolution(self, n=10):
        rho_points = len(self.rho)
        rho_mask = slice(None, None, rho_points // n)
        time = self.time

        ax = plt.gca()
        ax.plot(self.time, self.emissivity[:, rho_mask])
        ax.set_xlabel(r'$t\ [\mathrm{s}]$')
        ax = plt.gca()

if __name__ == '__main__':
    import doctest
    doctest.testmod()

    shot = 42661
    d = inverted_data(shot)
    plt.figure(1); plt.clf()
    d.select_time(0.5, 1.0).plot_contour()
    plt.draw()
    plt.show()


