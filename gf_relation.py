import os
import numpy as np
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt

import strahl
import gti
import ppfit
from scipy.interpolate import LSQUnivariateSpline, UnivariateSpline

def epsilon_prime(res):
    """
    Return \epsilon'.

    time_grid, rho_grid, epsilonp
    """
    rho_pol = res['rho_poloidal_grid']
    time = res['time']

    sxr = res['sxr_radiation'][:,-1,:]
    impdens = res['total_impurity_density']

    rho_mask = rho_pol < 0.5

    center = impdens[:, rho_mask]
    maxpos = center.argmax()
    maxpos = np.unravel_index(maxpos, center.shape)

    time_eq = maxpos[0]
    print 't_eq=', time[time_eq], '(%d)' % time_eq
    assert np.all(impdens[time_eq] !=0)
    epsilonp = sxr[time_eq] / impdens[time_eq]
    return rho_pol, epsilonp


def plot_epsilon_prime(epsilon):
    ax = plt.gca()
    e_rho, e_data = epsilon
    ax.plot(e_rho, e_data.T)
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r"$\epsilon'$")


def epsilon_on_new_radius_grid(epsilon, inverted_rho):
    e_rho, e_data = epsilon
    epsilon_new = np.interp(inverted_rho, e_rho, e_data)
    return inverted_rho, epsilon_new


def remove_offset(inversion, time_bbox):
    time = inversion.time
    time_mask = (time_bbox[0] <= time) & (time < time_bbox[1])
    offset = inversion.emissivity[time_mask, :]
    offset = offset.mean(0)
    offset = offset[np.newaxis]
    new_emissivity = inversion.emissivity - offset

    d = type(inversion)(inversion.shot, inversion.rho, inversion.time,
            new_emissivity)
    return d


class DataSmoother(object):
    def __init__(self, raw_data, dt=0.015):
        self.time = raw_data.time
        self.rho = raw_data.rho
        self.data = raw_data.emissivity

        self.test_rho = [0, 0.3, 0.6, 0.8]
        self.test_times = [0.72, 0.75]
        self.dt = dt


        dndt = np.zeros_like(self.data)
        n = np.zeros_like(self.data)
        for ir, r in enumerate(self.rho):
            y, dy = self._time_derivative(ir)
            n[:, ir] = y
            dndt[:, ir] = dy

        self.n = n
        self.dndt = dndt
        self.dndr = np.gradient(n)[1]/np.gradient(self.rho)

        self.dndr[:,0] = 0

    def _time_derivative(self, rho_index):
        time = self.time
        signal = self.data[:, rho_index]

        tt = np.array_split(time, self._pieces())
        t = np.array([np.mean(i) for i in tt])
        s = LSQUnivariateSpline(time, signal, t, k=3)
        return s(time), s(time,1)

    def _pieces(self):
        time = self.time
        dt = self.dt
        dt_sample = time[1] - time[0]
        k_dt = dt // dt_sample
        pieces = len(time) // k_dt
        return pieces

    def plot_time_evolution(self):
        ax = plt.gca()
        for i in np.searchsorted(self.rho, self.test_rho):
            y, dy = self._time_derivative(i)
            label = r'$\rho=%1.2f$' % self.rho[i]
            line, = ax.plot(self.time, self.data[:, i], '-', label=label)
            ax.plot(self.time, y, lw=2, color='black')

    def plot_dndt(self):
        ax = plt.gca()
        for i in np.searchsorted(self.rho, self.test_rho):
            y, dy = self._time_derivative(i)
            label = r'$\rho=%1.2f$' % self.rho[i]
            ax.plot(self.time, dy, label=label)

    def plot_profiles(self):
        ax = plt.gca()
        for i in np.searchsorted(self.time, self.test_times):
            label = r'$t=%1.2f\ \mathrm{s}$' % self.time[i]
            line, = ax.plot(self.rho, self.data[i,:], 'o', label=label)
            ax.plot(self.rho, self.n[i,:], color='black')

        ax.set_ylabel(r'$n_\mathrm{imp}$')
        ax.set_xlabel(r'$\rho_\mathrm{vol}$')

    def plot_dndr(self):
        ax = plt.gca()
        for i in np.searchsorted(self.time, self.test_times):
            label = r'$t=%1.2f\ \mathrm{s}$' % self.time[i]
            line, = ax.plot(self.rho, self.dndr[i,:], '-', label=label)

    def check_time_derivatives(self):
        f = plt.gcf()
        f.clf()

        ax = f.add_subplot(211)
        self.plot_time_evolution()
        ax.legend()

        ax = f.add_subplot(212, sharex=ax)
        self.plot_dndt()
        ax.legend()
        f.canvas.draw()

    def check_spatial_derivatives(self):
        f = plt.gcf()
        f.clf()

        ax = f.add_subplot(211)
        self.plot_profiles()
        ax.legend()

        ax = f.add_subplot(212)
        self.plot_dndr()
        ax.legend()
        f.canvas.draw()

    def d_dr(self, times):
        time_mask = (times[0] <= self.time) & (self.time < times[1])
        return self.dndr[time_mask]

    def d_dt(self, times):
        time_mask = (times[0] <= self.time) & (self.time < times[1])
        return self.dndt[time_mask]

    def smooth_n(self, times):
        time_mask = (times[0] <= self.time) & (self.time < times[1])
        selected_time = self.time[time_mask]
        return self.rho, selected_time, self.n[time_mask]


class GradientFlux(object):
    def __init__(self, smooth_data, time_bbox, rho_pol=None):
        rho, time, n = smooth_data.smooth_n(time_bbox)
        dndt = smooth_data.d_dt(time_bbox)
        dndr = smooth_data.d_dr(time_bbox)

        self.rho_vol = rho
        self.rho_pol = rho_pol
        self.time = time
        self.dndt = dndt
        self.dndr = dndr
        self.n = n
        self.time_bbox = time_bbox

        self.gradient = - dndr / n
        self.flux = self.Gamma() / n

        self.gradient *= 100 #[1/m]
        self.flux /= 100 #[m/s]

    def Gamma(self):
        """
        Calculate the integral \int_0^r{ r' dndt dr'}.
        """
        r = self.rho_vol
        dndt = self.dndt
        f = -cumtrapz(dndt * r, r)
        f /= r[:-1]

        ret = np.zeros_like(dndt)
        ret[:, 1:] = f
        return ret

    def plot_gf(self, rho=0):
        ax = plt.gca()

        rho_index = np.searchsorted(self.rho_pol, rho)
        label = r'$\rho_\mathrm{pol}=%1.2f$, ' % self.rho_pol[rho_index]
        label += r'$\rho_\mathrm{vol}=%1.2f$' % self.rho_vol[rho_index]
        ax.plot(self.gradient[:,rho_index], self.flux[:, rho_index], '.',
                label=label)

        D, v = self._fit_Dv(rho_index)
        x =  self.gradient[:,rho_index]
        ax.plot(x, D * x + v)

        ax.set_xlabel('$-(\mathrm{d}n/\mathrm{d}r)/n\ [\mathrm{1/m}]$')
        ax.set_ylabel('$\Gamma/n\ [\mathrm{m/s}]$')
        ax.plot([0],[0],'kx')

        indices = [0, -1]
        texts = [ r'$t_\mathrm{b}$', r'$t_\mathrm{e}$']
        for index, text in zip(indices, texts):
            xy = (self.gradient[index, rho_index], self.flux[index, rho_index])
            ax.annotate(text, xy)
        ax.legend()

    def _fit_Dv(self, rho_index):
        D, v = np.polyfit(self.gradient[:, rho_index],
                self.flux[:, rho_index], 1)
        return D, v

    def Dv_profile(self):
        r, D, v = [], [], []

        for i, rho in enumerate(self.rho_pol):
            try:
                D_, v_ = self._fit_Dv(i)
            except TypeError:
                continue

            r.append(rho)
            D.append(D_)
            v.append(v_)

        r = np.array(r)
        D = np.array(D)
        v = np.array(v)

        p = ppfit.ppolyfit(r, D, pieces=6, left=(0, [0,0]))
        Ds = ppfit.ppolyval(p, r)
        p = ppfit.ppolyfit(r, v, pieces=6, left=(0, [0,0]))
        vs = ppfit.ppolyval(p, r)
        return r, Ds, vs


def from_strahl_result(inversion, strahl_result, parameters):
    """
    Reconstruct D, v profiles based on the method presented in Sertoli et al.
    PPCF 2011.

    Parameters
    ----------
    inversion : InversionData object
        Tomographic inversion of soft X-ray emissivity on which the analysis
        is based on.  The offset before the impurity injection must be already
        removed.
    strahl_result : dictionary
        Result of a STRAHL simulation.
    """
    influx_bbox = parameters['influx_bbox']
    background_bbox = parameters['background_bbox']

    rho = inversion.rho
    rho_vol = np.interp(rho, strahl_result['rho_poloidal_grid'],
            strahl_result['radius_grid'])
    rho_pol = np.interp(rho, strahl_result['rho_poloidal_grid'],
            strahl_result['rho_poloidal_grid'])

    epsilon = epsilon_prime(strahl_result)
    epsilon = epsilon_on_new_radius_grid(epsilon, rho_pol)

    impurity_density = gti.InversionData(inversion.shot, rho_vol,
            inversion.time, inversion.emissivity / epsilon[1])
    impurity_density = remove_offset(impurity_density, background_bbox)

    s = DataSmoother(impurity_density)
    gf = GradientFlux(s, influx_bbox, rho_pol=rho_pol)
    return s, gf, epsilon


def save_profiles(r, D, v, limits=(0,1)):
    mask = (limits[0] < r) & (r < limits[1]) & (D > 0)
    np.savetxt('D_profile.txt', D[mask])
    np.savetxt('r_profile.txt', r[mask])
    np.savetxt('v_profile.txt', v[mask])


if __name__ == '__main__':
    working_directory = './wk'
    of = os.path.join(working_directory, 'result', 'Arstrahl_result.dat')
    res = strahl.viz.read_results(of)
    inversion = gti.inverted_data(42314).select_time(0.7, 1.0)
    parameters = dict(
        influx_bbox = (0.721, 0.75),
        background_bbox = (0.70, 0.705),
    )

    s, gf, epsilon = from_strahl_result(inversion, res, parameters)

    plt.figure(21); plt.clf()
    s.check_time_derivatives()
    plt.axvspan(*gf.time_bbox, color='black', alpha=0.2)
    plt.draw()

    plt.figure(22); plt.clf()
    s.check_spatial_derivatives()
    plt.draw()

    r, D, v = gf.Dv_profile()

    plt.figure(23); plt.clf()
    for rho in s.test_rho:
        gf.plot_gf(rho)
    plt.legend(loc='best')
    plt.draw()

    plt.figure(24); plt.clf()
    plt.subplot(211)
    plt.plot(r, D, '-o')
    plt.axhline(y=0, color='black')

    plt.subplot(212)
    plt.plot(r, v, '-o')
    plt.axhline(y=0, color='black')
    plt.draw()

    plt.figure(25); plt.clf()
    plot_epsilon_prime(epsilon)

    save_profiles(r, D, v)
    plt.show()

