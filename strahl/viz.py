import numpy as np
import matplotlib.pyplot as plt

from scipy.io import netcdf_file

_ax = dict(time=0, total_radiation=-1)

def read_results(filename):
    f = netcdf_file(filename, 'r')
    variables = f.variables
    out = {}
    for variableName, variable in variables.iteritems():
        out[variableName] = np.array(variable)

    out['large_radius'] = f.large_radius[0]
    out['small_radius'] = f.small_radius[0]

    out['anomal_diffusion'] /= 1e4
    out['anomal_drift'] /= 1e2
    return out


radial_coordinate = 'rho_poloidal'

def radial_grid(res):
    if radial_coordinate == 'rho_volume':
        return res['radius_grid']
    elif radial_coordinate == 'rho_poloidal':
        return res['rho_poloidal_grid']
    else:
        raise NotImplementedError('%s grid'% radial_grid)


def set_xaxis_rho():
    ax = plt.gca()
    if radial_coordinate == 'rho_volume':
        ax.set_xlabel(r'$\rho_\mathrm{vol}\ \mathrm{[cm]}$')
    elif radial_coordinate == 'rho_poloidal':
        ax.set_xlabel(r'$\rho_\mathrm{pol}$')
    else:
        raise NotImplementedError('%s grid'% radial_grid)



def plot_overview(res):
    shape_ = (4,4)

    # first row
    ax = plt.subplot2grid(shape_, (0,0), colspan=2)
    plot_electron_density(res)
    ax = plt.subplot2grid(shape_, (0,2), colspan=2)
    plot_electron_temperature(res)
    ax.yaxis.set_label_position('right')
    ax.yaxis.set_ticks_position('right')

    # second row
    ax = plt.subplot2grid(shape_, (1,0), colspan=2)
    plot_influx_through_valve(res)
    ax = plt.subplot2grid(shape_, (1,2), colspan=2)
    plot_volume(res)
    ax.yaxis.set_label_position('right')
    ax.yaxis.set_ticks_position('right')

    # third row
    ax = plt.subplot2grid(shape_, (2,0), colspan=2)
    plot_diffusion(res)

    ax = plt.subplot2grid(shape_, (2,2), colspan=2)
    plot_pinch(res)
    ax.yaxis.set_label_position('right')
    ax.yaxis.set_ticks_position('right')

    # fourth row
    ax = plt.subplot2grid(shape_, (3,0), colspan=2)
    plot_total_impurity_density(res)
    decimate_plotted_lines()

    ax = plt.subplot2grid(shape_, (3,2), colspan=2)
    plot_sxr(res)
    decimate_plotted_lines()
    ax.yaxis.set_label_position('right')
    ax.yaxis.set_ticks_position('right')
    yfmt = ax.yaxis.get_major_formatter()
    yfmt.set_powerlimits((-3,4))
    #plt.subplots_adjust(left=0.07, right=0.93, top=0.93, bottom=0.07,
    #        wspace=0.45, hspace=0.25)
 

def plot_background(res):
    plt.clf()
    plt.subplot(211)
    plot_electron_density(res)

    plt.subplot(212)
    plot_electron_temperature(res)

    plt.gcf().canvas.set_window_title('STRAHL: background plasma parameters')


def plot_transport_profiles(res):
    plt.clf()
    plt.subplot(211)
    plot_diffusion(res)

    plt.subplot(212)
    plot_pinch(res)

    plt.gcf().canvas.set_window_title('STRAHL: transport properties')


def plot_impurity(res):
    plt.clf()
    plt.subplot(211)
    plot_total_impurity_density(res)
    decimate_plotted_lines()

    plt.subplot(212)
    plot_sxr(res)
    decimate_plotted_lines()

    plt.gcf().canvas.set_window_title('STRAHL: impurity')


def plot_electron_density(res):
    ax = plt.gca()

    r = radial_grid(res)
    ne = res['electron_density']
    ax.plot(r, ne.mean(_ax['time']))

    ax.set_ylabel('$n_\mathrm{e}\ \mathrm{[cm^{-3}]}$')
    set_xaxis_rho()

    ax.grid(True)


def plot_electron_temperature(res):
    ax = plt.gca()

    r = radial_grid(res)
    te = res['electron_temperature']
    ax.plot(r, te.mean(_ax['time']))

    ax.set_ylabel('$T_\mathrm{e}\ \mathrm{[eV]}$')
    set_xaxis_rho()

    ax.grid(True)


def plot_diffusion(res):
    ax = plt.gca()

    r = radial_grid(res)
    ax.plot(r, res['anomal_diffusion'].mean(_ax['time']))
    ax.set_ylabel(r'$D\ [\mathrm{m^2/s}]$')
    set_xaxis_rho()

    ax.grid(True)


def plot_pinch(res):
    ax = plt.gca()

    r = radial_grid(res)
    pro = res['anomal_drift']
    ax.plot(r, res['anomal_drift'].mean(_ax['time']))
    ax.set_ylabel(r'$V\ [\mathrm{m/s}]$')
    ax.grid(True)
    set_xaxis_rho()


def plot_sxr(res):
    ax = plt.gca()

    r = radial_grid(res)
    sxr = res['sxr_radiation']
    ax.plot(r, sxr[:,_ax['total_radiation'],:].T, '-')
    ax.set_ylabel(r'$E_\mathrm{SXR}\ [\mathrm{W/cm^3}]$')
    set_xaxis_rho()

    ax.grid(True)


def plot_total_impurity_density(res):
    ax = plt.gca()

    r = radial_grid(res)
    impdens = res['total_impurity_density']
    ax.plot(r, impdens.T, '-')

    ax.set_ylabel('$n_\mathrm{imp}\ [\mathrm{cm^{-3}}]$')
    set_xaxis_rho()
    ax.grid(True)

    for txt, line in zip(legend_from_time(res['time']), ax.lines):
        line.set_label(txt)


def plot_influx_through_valve(res):
    ax = plt.gca()

    ax.plot(res['time'], res['influx_through_valve'])

    ax.set_xlabel('$t\ [\mathrm{s}$]')
    ax.set_ylabel('$\Gamma_\mathrm{valve}\ [\mathrm{cm^{-1}s^{-1}}]$')
    ax.grid(True)


def plot_volume(res):
    ax = plt.gca()

    r = radial_grid(res)
    rho_volume = res['radius_grid']
    R_axis = res['large_radius']
    volume = rho_volume**2 * 2* np.pi**2 * R_axis #cm3
    volume /= 1e6 #m3

    ax.plot(r, volume)

    set_xaxis_rho()
    ax.grid(True)
    ax.set_ylabel('$V\ [\mathrm{m^3}]$')


def decimate_plotted_lines():
    ax = plt.gca()

    for i, line in enumerate(ax.lines):
        if i % 10:
            line.set_visible(False)
            line.set_label('')


def legend_from_time(time_vector):
    text = '$t=%1.2f\ \mathrm{s}$'
    return [ text % t for t in time_vector]


if __name__ == '__main__':
    of = '/home/dwagner/work/strahl/result/strahl_result.dat'
    res = read_results(of)

    plt.figure(1); plt.clf()
    plot_overview(res)
    plt.draw()

    plt.show()

