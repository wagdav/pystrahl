import numpy as np

import matplotlib.pyplot as plt

from strahl.geometry.chord import integrate_along_chord
from strahl.viz import read_results
from strahl.equilibrium import MockUpEquilibrium, plot_equilibrium
from strahl.geometry.tcv import dmpx_chords

def line_integrated_measurements(res, eq):
    chords = dmpx_chords()
    x, y, rho = eq.get_rho()

    E_sxr = []
    r = res['rho_poloidal_grid']
    for time_index in xrange(len(res['time'])):
        E_sxr.append(np.interp(rho, r, res['sxr_radiation'][time_index,-1,:], right=0.0))

    E_sxr = np.array(E_sxr)

    profiles = []
    for chord in chords:
        p = []
        for time_index in xrange(len(res['time'])):
            d = E_sxr[time_index]
            p.append(integrate_along_chord(chord, x, y, d))

        profiles.append(np.array(p))

    return np.array(profiles)


def plot_geometry(eq):
    ax = plt.gca()

    chords = dmpx_chords()

    labels = []
    for i,c in enumerate(chords):
        labels.append(r'$%d$' % i)

    for chord, label in zip(chords, labels):
        ax.plot(*chord, color='black', ls='--')
        ax.annotate(label, (chord[0][1], chord[1][1]))


if __name__ == '__main__':
    eq = MockUpEquilibrium()
    of = '/home/dwagner/work/strahl/result/strahl_result.dat'
    res = read_results(of)

    profiles = line_integrated_measurements(res, eq)

    plt.figure(1); plt.clf()
    ax = plt.gcf().add_subplot(121)
    ax.plot(res['time'], profiles.T)

    chords = dmpx_chords()
    labels = []
    for i,c in enumerate(chords):
        labels.append(r'chord %d' % i)
    ax.legend(labels)
    ax.set_xlabel(r'$t\ [\mathrm{s}]$')
    plt.draw()

    ax = plt.gcf().add_subplot(122)
    plot_geometry(eq)
    plot_equilibrium(eq)
    plt.draw()

    plt.show()

