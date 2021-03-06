import numpy as np
import matplotlib.pyplot as plt

from geometry.chord import integrate_along_chord

def line_integrated_measurements(res, eq, chord):
    x, y, psi_norm = eq.get_psi_contours()

    f_rho = res['rho_poloidal_grid']

    E_sxr = []
    for time_index, time in enumerate(res['time']):
        f = res['sxr_radiation'][time_index,-1,:]
        E_sxr.append(np.interp(np.sqrt(psi_norm), f_rho, f, right=0.0))

    E_sxr = np.array(E_sxr)

    times, chord_intensity = [], []
    for time_index, t in enumerate(res['time']):
        d = E_sxr[time_index]
        chord_intensity.append(integrate_along_chord(chord, x, y, d))
        times.append(t)

    times = np.array(times)
    chord_intensity = np.array(chord_intensity)

    return times, chord_intensity


if __name__ == '__main__':
    import doctest
    doctest.testmod()
