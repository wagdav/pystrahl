import os
import numpy as np
import matplotlib.pyplot as plt

import strahl
from strahl.diagnostics import dmpx
from crpppy.tcv_geom import vessel_patch
import crpppy.diagnostics.dmpx as tcv_dmpx
import crpppy.diagnostics.thomson as thomson
import ppfit

import gf_relation
import gti

working_directory = './wk'

def get_rho_volume(eq):
    r = np.linspace(0, 1, 20)
    volume = eq.get_volume(r)
    rho_vol = np.sqrt(volume/(2 * np.pi**2 * eq.major_radius))
    rho_vol_LCFS = rho_vol[-1]

    return rho_vol, rho_vol_LCFS


def velocity_from_zero_flux(rho, D, ne, eq):
    """
    >>> rho = np.linspace(0,1,20)
    >>> D = rho**2 + 0.02
    >>> ne = 1e20 * (1 -rho**2) + 1e18
    >>> v = velocity_from_zero_flux(rho, D, ne)
    """
    grho = eq.get_grho(rho)
    dne = np.gradient(ne, np.diff(rho)[0])

    v = 100 * D * grho * dne/ne

    v[rho > 0.9] = 0
    return v


def plot_geometry():
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.add_patch(vessel_patch())
    dmpx.plot_geometry()
    eq.plot()
    plt.draw()


def syntetic_chords(strahl_result, equilibrim, chord_indices):
    res = {}
    chords = dmpx.measurement_chords()
    for ci in chord_indices:
        chord = chords[ci]
        time, profile = dmpx.line_integrated_measurements(strahl_result,
                equilibrim, chord)
        res[ci] = profile
    return time, res


def plot_chord_evolution(measured, simulated,
        chord_indices=[32], time_offset=0.0):
    t, yy = simulated
    time_offset += measured.time[0]
    scaling_factor = measured.central_max() / yy[32].max()

    ax = plt.gca()
    for c in chord_indices:
        label = str(c)
        y = yy[c]

        line, = ax.plot(measured.time, measured.data[c], lw=0.1)
        ax.plot(t + time_offset, scaling_factor * y, lw=2,
                color=line.get_color(), label=label)

        maximum = np.max(y)
        ind = np.argmax(y)
        ax.plot(t[ind] + time_offset, scaling_factor*maximum, 'ko')

    ax.set_xlabel(r'$t\ [\mathrm{s}]$')


class DMPX_data(object):
    def __init__(self, time, data):
        self.time = time
        self.data = data
        if data.shape[1] != np.alen(time):
            print self.time.shape, self.data.shape
            raise AssertionError('Size mismatch')

    def central_max(self):
        central_chord = 32
        t, y = self._smooth_chord(central_chord)
        return max(y)

    def _smooth_chord(self, chord_index):
        signal = self.data[chord_index]
        p = ppfit.ppolyfit(self.time, signal, 20, deg=3)
        smooth = ppfit.ppolyval(p, self.time)
        return self.time, smooth

    def select_time(self, begin, end):
        time_mask = (begin <= self.time) & (self.time < end)
        new_time = self.time[time_mask]
        new_data = self.data[:, time_mask]
        return type(self)(new_time, new_data)

    def remove_offset(self, begin, end):
        offset = self.select_time(begin, end)
        offset = offset.data.mean(axis=1)
        new_data = self.data - offset[:,np.newaxis]
        return type(self)(self.time, new_data)

    def plot(self):
        ax = plt.gca()
        ax.plot(self.time, self.data.T)


def dmpx_from_shot(shot):
    time, data = tcv_dmpx.get_data(shot)
    return DMPX_data(time, data)


class StrahlSimulation(object):
    def __init__(self, equilibrium, dmpx_data, thomson_data, inversion=None):
        self.equilibrium = equilibrium
        self.dmpx_data = dmpx_data
        self.thomson_data = thomson_data
        self.inversion = inversion

        self.setup()

    def setup(self, from_file=False):
        eq = self.equilibrium
        thomson_data = self.thomson_data
        params = strahl.defaultParams()

        t, flx = strahl.rectangular_pulse_with_decay(length=5e-3,
                max_value=5.0e24, tau=160e-3)

        params['impurity.influx'] = (t, flx)
        params['numerical.time.dt'] = 3e-4
        params['numerical.time.final'] = 0.6
        params['numerical.grid.k'] = 8

        rho_vol, rho_vol_LCFS = get_rho_volume(eq)
        rho_pol = np.linspace(0,1,20)

        _, ne = thomson_data.fit_density(rho_pol)
        _, Te = thomson_data.fit_temperature(rho_pol)
        ne /= 1e6

        if from_file:
            r = np.loadtxt('r_profile.txt')
            D = np.loadtxt('D_profile.txt')
            v = np.loadtxt('v_profile.txt')
            D = np.interp(rho_pol, r, D)
            v = np.interp(rho_pol, r, v)

            v[0] = 0.
        else:
            D = 1* rho_pol**2 + 0.001
            D = 10 * np.ones_like(rho_pol)
            v = np.zeros_like(D)
            v = velocity_from_zero_flux(rho_pol, D, ne, eq)
            v = -150 * np.ones_like(rho_pol)

        params['geometry.rho_volume'] = rho_vol * 100
        params['background.rho_poloidal'] = rho_pol
        params['background.electron_density'] = ne
        params['background.electron_temperature'] = Te
        params['geometry.rho_volume_at_lcfs'] = rho_vol_LCFS * 100
        params['background.decay_length'] = rho_vol_LCFS * 0.1

        params['impurity.sol_width'] = rho_vol_LCFS * 0.1
        params['impurity.convection_velocity'] = v
        params['impurity.diffusion_coefficient'] = D

        strahl.create_input(params, working_directory)
        self.params = params

    def run(self):
        curdir = os.getcwd()
        os.chdir(working_directory)
        os.system('./strahl a')
        os.chdir(curdir)

        of = os.path.join(working_directory,'result','Arstrahl_result.dat')
        self.result = strahl.viz.read_results(of)

    def viz(self, offset=0):
        offset *= 10

        strahl_result = self.result
        eq = self.equilibrium
        dmpx_data = self.dmpx_data

        plt.figure(offset + 1); plt.clf()
        strahl.viz.plot_overview(strahl_result)
        plt.draw()

        plt.figure(offset + 2); plt.clf()
        ax = plt.gcf().add_subplot(111)

        chords = [32, 36, 25]
        simulated_chords = syntetic_chords(strahl_result, eq, chords)
        plot_chord_evolution(dmpx_data, simulated_chords, chords,
                time_offset=0.015)
        plt.ylim(ymin=0)
        plt.legend()
        plt.draw()

    def gf_loop(self, loop=3, resume=False):
        assert self.inversion is not None, 'Inversion data is missing.'
        if not resume:
            plt.figure(100); plt.clf()

        inversion = self.inversion

        self.setup(from_file=resume)
        D = self.params['impurity.diffusion_coefficient']
        v = self.params['impurity.convection_velocity']
        r = self.params['background.rho_poloidal']
        plot_Dv(r, D, v)
        for i in xrange(loop): # GF-loop
            self.run()

            s, gf, epsilon = gf_relation.from_strahl_result(inversion,
                    self.result, self.gf_parameters)
            r, D, v = gf.Dv_profile()

            r = np.interp(r, self.result['radius_grid'],
                    self.result['rho_poloidal_grid'])
            gf_relation.save_profiles(r, D, v)

            self.setup(from_file=True)

            plt.figure(100)
            plot_Dv(r, D, v)
            plt.draw()

        self.gf = gf


def plot_Dv(r, D, v):
    fig = plt.gcf()
    ax = fig.add_subplot(211)
    ax.plot(r, D, '-o')
    ax.axhline(y=0, color='black')
    ax.set_title('D,v control panel')

    ax = fig.add_subplot(212)
    ax.plot(r, v, '-o')


def simulation_from_shot(shot):
    time_bbox = db[shot]['time_bbox']
    background_bbox = db[shot]['background_bbox']

    equilibrium_time = sum(time_bbox)/2.0
    print 'Creating LiqueEqulibirium...'
    eq = strahl.LiuqeEquilibrium(shot, equilibrium_time)

    print 'Creating DMPXData...'
    dmpx_data = dmpx_from_shot(shot).select_time(*time_bbox)
    dmpx_data = dmpx_data.remove_offset(*background_bbox)

    print 'Creating ThomsonData...'
    thomson_data = thomson.thomson_data(shot)
    thomson_data = thomson_data.select_time(*time_bbox)

    print 'Creating InversionData'
    inversion = gti.inverted_data(shot).select_time(*time_bbox)

    sim = StrahlSimulation(eq, dmpx_data, thomson_data, inversion)
    sim.gf_parameters = db[shot]
    return sim


db = {
    42661 : dict(
        time_bbox = (0.5, 1.0),
        background_bbox = (0.5, 0.505),
        influx_bbox = (0.515, 0.55),
    ),
    42314 : dict( #density ne = 1fr, Ip = 125 kA 
        time_bbox = (0.7, 1.0),
        background_bbox = (0.7, 0.705),
        influx_bbox = (0.715, 0.75),
    )
}

try:
    sim
except NameError:
    sim = simulation_from_shot(42661)

plt.show()
