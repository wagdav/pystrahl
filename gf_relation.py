import os
import numpy as np
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt

import strahl
import gti
import ppfit


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

    #if maxpos[1] != 0:
    #    print 'Cannot find the maximal impurity density in the center'
    #    print maxpos[1]
    #    raise AssertionError

    time_eq = maxpos[0]

    epsilonp = sxr[time_eq] / impdens[time_eq]
    return rho_pol, epsilonp


def plot_epsilon_prime():
    ax = plt.gca()
    e_time, e_rho, e_data = epsilon_prime(res)
    ax.plot(e_rho, e_data.T)
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r"$\epsilon'$")


def Gamma(dndt, r):
    """
    Calculate the integral \int_0^r{ r' dndt dr'}.
    """
    f = -cumtrapz(dndt * r, r)
    f /= cumtrapz(np.ones_like(r), r)
    return f


def epsilon_on_new_radius_grid(epsilon, inverted_rho):
    e_rho, e_data = epsilon
    epsilon_new = np.interp(inverted_rho, e_rho, e_data)
    return inverted_rho, epsilon_new


def calculate_impurity_density(epsilon, inversion):
    offset = inversion.emissivity[0:100,:]
    offset = offset.mean(0)
    offset = offset[np.newaxis]

    impurity_density = (inversion.emissivity - offset) / epsilon[1]
    return impurity_density


def smooth_impurity_density(impurity_density):
    time = impurity_density.time
    rho = impurity_density.rho

    dt_sample = time[1] - time[0]
    dt_smooth = 0.02

    k_dt = dt_smooth // dt_sample
    pieces = len(time) // k_dt
    print pieces

    dndt = []
    smooth = []
    for rho_index, r in enumerate(rho):
        signal = impurity_density.emissivity[:,rho_index]

        p = ppfit.ppolyfit(time, signal, pieces, deg=3)
        dp = ppfit.ppolyder(p)

        y = ppfit.ppolyval(p, time)
        dy = ppfit.ppolyval(dp, time)

        dndt.append(dy)
        smooth.append(y)

    dndt = np.array(dndt)
    smooth = np.array(smooth).T
    dndt = gti.InversionData(impurity_density.shot, rho, time, dndt)
    smooth = gti.InversionData(impurity_density.shot, rho, time, smooth)

    return smooth, dndt


class GradientFlux(object):
    def __init__(self, impurity_density_exp):
        self.impurity_density_exp = impurity_density_exp
        self.eval_gradient_flux()

    def _smooth(self):
        smooth, dndt = smooth_impurity_density(self.impurity_density_exp)
        return smooth, dndt

    def eval_gradient_flux(self):
        impurity_density, dndt = self._smooth()
        rho = self.impurity_density_exp.rho
        time = self.impurity_density_exp.time
        rho_vol = np.interp(rho, res['rho_poloidal_grid'], res['radius_grid'])
        G = Gamma(dndt.emissivity.T, rho_vol)

        dndr = np.zeros_like(impurity_density.emissivity)
        for i, t in enumerate(time):
            dndr[i,:] = np.gradient(impurity_density.emissivity[i])/np.gradient(rho_vol)

        x = dndr[:,:-1] / impurity_density.emissivity[:,:-1]
        y = G / impurity_density.emissivity[:, :-1]

        self.gradient = x
        self.flux = y
        return x, y

    def plot_gradient(self, rho_index=0):
        ax = plt.gca()
        time = self.impurity_density_exp.time
        rho = self.impurity_density_exp.rho

        ax.plot(time, self.gradient[:, rho_index])
        ax.text(0.95, 0.95, r'$\rho=%1.2f$' % rho[rho_index],
                transform=ax.transAxes, ha='right', va='top')
        ax.figure.canvas.draw()

    def plot_flux(self, rho_index=0):
        ax = plt.gca()
        time = self.impurity_density_exp.time
        rho = self.impurity_density_exp.rho
        ax.plot(time, self.flux[:, rho_index])

        ax.text(0.95, 0.95, r'$\rho=%1.2f$' % rho[rho_index],
                transform=ax.transAxes, ha='right', va='top')
        ax.figure.canvas.draw()

    def plot_time_evolution(self):
        ax = plt.gca()

        d = self.impurity_density_exp
        ax.plot(d.time, d.emissivity[:, ::5], '.')
        smooth, dndt = self._smooth()
        ax.plot(smooth.time, smooth.emissivity[:, ::5], lw=2)

    def plot_gf(self, rho_index=0):
        ax = plt.gca()
        ax.plot(self.gradient[:,rho_index], self.flux[:, rho_index], '.',
                label=str(inversion.rho[rho_index]))

        D, v = self._fit_Dv(rho_index)
        x =  self.gradient[:,rho_index]
        ax.plot(x, -D * x+ v)

    def _fit_Dv(self, rho_index):
        D, v = np.polyfit(self.gradient[:,rho_index], self.flux[:,rho_index], 1)
        D = -D
        return D, v

    def Dv_profile(self, left=0, right=0.6):
        D = []
        v = []
        r = []

        for i, rho in enumerate(self.impurity_density_exp.rho):
            if rho < left or rho > right: continue
            try:
                D_, v_ = self._fit_Dv(i)
            except TypeError:
                continue

            r.append(inversion.rho[i])
            D.append(D_)
            v.append(v_)

        r = np.array(r)
        D = np.array(D) / 1e4
        v = np.array(v) / 1e2
        return r, D, v


working_directory = './wk'
of = os.path.join(working_directory, 'result', 'Arstrahl_result.dat')
res = strahl.viz.read_results(of)

inversion = gti.inverted_data(42661).select_time(0.5, 1.0)

epsilon = epsilon_prime(res)
epsilon = epsilon_on_new_radius_grid(epsilon, inversion.rho)

impurity_density_exp = gti.InversionData(inversion.shot, inversion.rho,
        inversion.time, inversion.emissivity / epsilon[1])

gf = GradientFlux(impurity_density_exp.select_time(0.52, 0.56))

r,D,v = gf.Dv_profile(left=0.0, right=0.6)

plt.figure(10); plt.clf()
plt.subplot(211)
plt.plot(r,D)

plt.subplot(212)
plt.plot(r,v)

plt.draw()
plt.show()

np.savetxt('D_profile.txt', D)
np.savetxt('r_profile.txt', r)
np.savetxt('v_profile.txt', v)
