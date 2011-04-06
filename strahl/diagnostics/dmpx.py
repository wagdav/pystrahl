import numpy as np

import matplotlib.pyplot as plt

from strahl.geometry.chord import integrate_along_chord

def line_integrated_measurements(res, eq):
    chords = measurement_chords()
    x, y, radius = eq.get_radius_contours()

    f_rhovol = res['radius_grid']
    radius_grid = np.linspace(0, eq.minor_radius, np.alen(f_rhovol))
    eq_rhovol = eq.get_rhovol(radius_grid)


    E_sxr = []
    for time_index in xrange(len(res['time'])):
        f = res['sxr_radiation'][time_index,-1,:]
        f_as_radius = np.interp(eq_rhovol, f_rhovol, f)

        E_sxr.append(np.interp(radius, radius_grid, f_as_radius, right=0.0))

    E_sxr = np.array(E_sxr)

    profiles = []
    for chord in chords:
        p = []
        for time_index in xrange(len(res['time'])):
            d = E_sxr[time_index]
            p.append(integrate_along_chord(chord, x, y, d))

        profiles.append(np.array(p))

    return np.array(profiles)


def plot_geometry():
    ax = plt.gca()

    chords = measurement_chords()

    labels = []
    for i,c in enumerate(chords):
        labels.append(r'$%d$' % i)

    for chord, label in zip(chords, labels):
        ax.plot(*chord, color='black', ls='--')
        ax.annotate(label, (chord[0][1], chord[1][1]))


def measurement_chords():
    x_start = 0.88
    y_start = -0.78
    x_end_ = np.linspace(0.6, 1.2, 8)
    y_end = 0.78

    chords = []
    for x_end in x_end_:
        chords.append(((x_start, x_end), (y_start, y_end)))

    return chords


if __name__ == '__main__':
    from strahl.viz import read_results
    from strahl.equilibrium import MockUpEquilibrium, plot_equilibrium

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

