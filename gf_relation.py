import numpy as np
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt

import strahl
import gti
from scipy import interpolate

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


def smooth_derivative(time, y, time_new):
    s = smoothing_parameter(time)

    tck = interpolate.splrep(time, y, s=s)
    ynew = interpolate.splev(time_new, tck, der=0)
    ynew1 = interpolate.splev(time_new, tck, der=1)

    return ynew, ynew1


def smoothing_parameter(time):
    m = np.alen(time)
    return int(0.2*m)


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

def epsilon_on_new_radius_grid(epsilon, inverted_rho):
    e_time, e_rho, e_data = epsilon

    epsilon_new = np.zeros((len(e_time), len(inverted_rho)))
    for i in xrange(len(e_time)):
        epsilon_new[i] = np.interp(inverted_rho, e_rho, e_data[i])

    return epsilon_new


def epsilon_on_new_time_grid(epsilon, inverted_time, res):
    time_res = res['time']
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



