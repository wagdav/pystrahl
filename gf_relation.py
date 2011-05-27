import os
import numpy as np
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt

import strahl
import gti
import ppfit

working_directory = './wk'


def epsilon_prime(res):
    """
    Return \epsilon'.

    time_grid, rho_grid, epsilonp
    """
    rho_pol = res['rho_poloidal_grid']
    time = res['time']

    sxr = res['sxr_radiation'][:,-1,:]
    impdens = res['total_impurity_density']

    epsilonp = sxr / impdens
    return time, rho_pol, epsilonp


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
    e_time, e_rho, e_data = epsilon

    epsilon_new = np.zeros((len(e_time), len(inverted_rho)))
    for i in xrange(len(e_time)):
        epsilon_new[i] = np.interp(inverted_rho, e_rho, e_data[i])

    return e_time, inverted_rho, epsilon_new


def epsilon_on_new_time_grid(epsilon, inverted_time):
    e_time, e_rho, e_data = epsilon

    epsilon_new = np.zeros((len(inverted_time), len(e_rho)))
    for i in xrange(len(e_rho)):
        epsilon_new[:,i] = np.interp(inverted_time, e_time, e_data[:,i])

    return inverted_time, e_rho, epsilon_new


def plot_gf(rho_index=0):
    ax = plt.gca()
    x, y = select_influx_phase(rho_index)
    ax.plot(x, y, '.', label=str(inverted_rho[rho_index]))


def do_fit(rho_index):
    x, y = select_influx_phase(rho_index)
    D, v = np.polyfit(x, y, 1)
    D = -D
    return D, v


def select_influx_phase(rho_index):
    xsel = x[time_mask][:, rho_index]
    ysel = y[time_mask][:, rho_index]

    mask = np.gradient(xsel) < 0
    xsel = xsel[mask]
    ysel = ysel[mask]
    return xsel, ysel


def remap_epsilon(res, inversion):
    inverted_rho, inverted_time, inverted = inversion

    epsilon = epsilon_prime(res)
    epsilon = epsilon_on_new_radius_grid(epsilon, inverted_rho)
    epsilon = epsilon_on_new_time_grid(epsilon, inverted_time)

    return epsilon


def calculate_impurity_density(epsilon, inversion):
    inverted_rho, inverted_time, inverted = inversion

    offset = inverted[0:100,:]
    offset = offset.mean(0)
    offset = offset[np.newaxis]
    return (inverted - offset) / epsilon[2][991]


def smooth_impurity_density(time, impurity_density):
    dndt = []
    impurity_density = []
    for rho_index in xrange(40):
        signal = impurity_density_exp[:,rho_index]

        p = ppfit.ppolyfit(time, signal, 20, deg=3)
        dp = ppfit.ppolyder(p)

        y = ppfit.ppolyval(p, time)
        dy = ppfit.ppolyval(dp, time)

        dndt.append(dy)
        impurity_density.append(y)

    dndt = np.array(dndt)
    impurity_density = np.array(impurity_density).T

    return impurity_density, dndt


def plot_impurity_density(rho_index=0):
    ax = plt.gca()
    ax.plot(time, impurity_density_exp[:, rho_index])
    ax.plot(time, impurity_density[:, rho_index])
    ax.plot(time[time_mask], impurity_density[time_mask][:, rho_index])

    ax.text(0.95, 0.95, r'$\rho=%1.2f$' % inverted_rho[rho_index],
            transform=ax.transAxes, ha='right', va='top')
    ax.figure.canvas.draw()


def plot_impurity_gradient(rho_index=0):
    ax = plt.gca()
    ax.plot(time, x[:, rho_index])
    ax.plot(time[time_mask], x[time_mask][:, rho_index])

    ax.text(0.95, 0.95, r'$\rho=%1.2f$' % inverted_rho[rho_index],
            transform=ax.transAxes, ha='right', va='top')
    ax.figure.canvas.draw()


def plot_impurity_flux(rho_index=0):
    ax = plt.gca()
    ax.plot(time, y[:, rho_index])
    ax.plot(time[time_mask], y[time_mask][:, rho_index])

    ax.text(0.95, 0.95, r'$\rho=%1.2f$' % inverted_rho[rho_index],
            transform=ax.transAxes, ha='right', va='top')
    ax.figure.canvas.draw()


of = os.path.join(working_directory, 'result', 'Arstrahl_result.dat')
res = strahl.viz.read_results(of)

inversion = gti.inverted_data(42661) # rho, time, data
epsilon = remap_epsilon(res, inversion)

impurity_density_exp = calculate_impurity_density(epsilon, inversion)
impurity_density, dndt = smooth_impurity_density(inversion[1], impurity_density_exp)

rho_vol = np.interp(inversion[0], res['rho_poloidal_grid'], res['radius_grid'])
G = Gamma(dndt.T, rho_vol)

# gradient
dndr = np.zeros_like(impurity_density)
for i in xrange(len(inversion[1])):
    dndr[i,:] = np.gradient(impurity_density[i])/np.gradient(rho_vol)
    l = dndr[i,:]

x = dndr[:,:-1] / impurity_density[:,:-1]
y = G / impurity_density[:, :-1]

time_range = (0.5 + 0.02, 0.5 + 0.04)
time_mask = (inversion[1] > time_range[0]) & (inversion[1] < time_range[1])

D = []
v = []
r = []
for i in xrange(3, 24):
    try:
        D_, v_ = do_fit(i)
    except TypeError:
        continue

    if D_ > 0 and D_ < 200 and abs(v_) < 2000:
        r.append(inversion[0][i])
        D.append(D_)
        v.append(v_)

r = np.array(r)
D = np.array(D)
v = np.array(v)

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

