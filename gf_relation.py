import numpy as np
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt

import strahl
import gti
from savitzky_golay import savitzky_golay

def epsilon_prime(res):
    rho_pol = res['rho_poloidal_grid']

    sxr = res['sxr_radiation'][:,-1,:]
    impdens = res['total_impurity_density']

    epsilonp = sxr / impdens
    return rho_pol, epsilonp


def plot_epsilon_prime():
    ax = plt.gca()
    e_rho, e_data = epsilon_prime(res)
    ax.plot(e_rho, e_data.T)
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r"$\epsilon'$")


def ddt(time, signal):
    window_size = 2001
    order = 4

    derivative = np.zeros_like(signal)

    for i in xrange(signal.shape[1]):
        derivative[:,i] = savitzky_golay(signal[:,1], window_size, order, deriv=1)
        derivative[:,i] /= -np.gradient(time)

    return derivative


def Gamma_A(rho, time, n):
    dndt = ddt(time, n) # time derivative

    gamma = np.zeros_like(dndt)
    for i in xrange(len(time)):
        a = cumtrapz(dndt[i] * rho, rho)
        r = cumtrapz(np.ones_like(rho),rho)
        gamma[i,:-1] = - a / r

    return gamma


def Gamma(dndt, r):
    """
    Calculate the integral \int_0^r{ r' dndt dr'}.
    """
    f = -cumtrapz(dndt * r, r)
    f /= cumtrapz(np.ones_like(r), r)
    return f


of = '/home/dwagner/work/strahl/result/strahl_result.dat'
res = strahl.viz.read_results(of)

inverted_rho, inverted_time, inverted = gti.inverted_data(42661)

rho_vol_res = res['radius_grid']
rho_res = res['rho_poloidal_grid']
rho_vol = np.interp(inverted_rho, rho_res, rho_vol_res)

# Map epsilon prime to *rho_vol*
e_rho, e_data = epsilon_prime(res)

time_res = res['time']
epsilon = np.zeros((len(time_res), len(inverted_rho)))
for i in xrange(len(time_res)):
    epsilon[i] = np.interp(inverted_rho, e_rho, e_data[i])

epsilon_i = np.zeros_like(inverted)
for i in xrange(len(inverted_rho)):
    epsilon_i[:,i] = np.interp(inverted_time, time_res, epsilon[:,i])

# Smooth the inverted signals in time
inverted_filtered = np.zeros_like(inverted)
for i in xrange(40):
    inverted_filtered[:,i] = savitzky_golay(inverted[:,i], 2001, 4, deriv=0)


# impurity density
offset = inverted_filtered[0:100,:]
offset = offset.mean(0)
offset = offset[np.newaxis]
impdens_exp = (inverted_filtered - offset) / epsilon_i



# flux
G = Gamma_A(rho_vol, inverted_time, impdens_exp)

# gradient
dndr = np.zeros_like(impdens_exp)
for i in xrange(len(inverted_time)):
    dndr[i,:] = np.gradient(impdens_exp[i])/np.gradient(rho_vol)
    l = dndr[i,:]
    l[l>0] = 0

x = dndr / impdens_exp
y = G / impdens_exp


time_range = (0.5 + 0.02, 0.5 + 0.05)
time_mask = (inverted_time > time_range[0]) & (inverted_time < time_range[1])



D = []
v = []
r = []
for i in xrange(len(inverted_rho)):
    D_, v_ = np.polyfit(x[time_mask][:,i], y[time_mask][:,i],1)
    D_ = -D_

    if D_ > 0.02:
        r.append(inverted_rho[i])
        D.append(D_)
        v.append(v_)

r = np.array(r)
D = np.array(D)
v = np.array(v)

plt.figure(10); plt.clf()

for i in [5, 10, 25]:
    plt.plot( x[time_mask][:,i], y[time_mask][:,i], label=str(inverted_rho[i]))

plt.legend()

plt.figure(11); plt.clf()
plt.subplot(211)
plt.plot(r,D)

plt.subplot(212)
plt.plot(r,v)

plt.draw()
plt.show()

np.savetxt('D_profile.txt', D)
np.savetxt('r_profile.txt', r)



