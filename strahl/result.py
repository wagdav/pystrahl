__all__ = ['load', 'Result']


import numpy as np
from scipy.io import netcdf_file


def load_old(filename):
    f = netcdf_file(filename, 'r')
    variables = f.variables
    out = {}
    for name in variables:
        out[name] = variables[name].data

    out['large_radius'] = f.large_radius
    out['small_radius'] = f.small_radius
    out['maximum_charge'] = f.maximum_charge
    out['element'] = f.species

    out['anomal_diffusion'] = out['anomal_diffusion'] / 1e4
    out['anomal_drift'] = out['anomal_drift'] / 1e2
    return Result(out)


class Result_old(object):
    def __init__(self, result):
        """
        Parameters
        ----------
        result : dict or str
        """

        if isinstance(result, dict):
            self.output = result
        elif isinstance(result, basestring):
            self.output = load(result)

    def radial_grid(self, coords='rhopol'):
        if coords == 'rhopol':
            return self.output['rho_poloidal_grid']
        elif coords == 'rhovol':
            return self.output['radius_grid']
        else:
            raise ValueError

    @property
    def rhopol(self):
        return self.radial_grid('rhopol')

    @property
    def rhovol(self):
        return self.radial_grid('rhovol')

    @property
    def time(self):
        return self.output['time']

    @property
    def central_impurity_density(self):
        """ Impurity density at the center in 1/cm3 shape = (radial_grid)"""
        return self.impurity_density[:, 0]

    @property
    def electron_density(self):
        """ electron density in 1/cm3 shape = (time, radial_grid)"""
        return self.output['electron_density']

    @property
    def electron_temperature(self):
        """ electron density in 1/cm3 shape = (time, radial_grid)"""
        return self.output['electron_temperature']

    @property
    def impurity_density(self):
        """ Impurity density in 1/cm3 shape = (time, radial_grid)"""
        return self.output['total_impurity_density']

    @property
    def diffusion_profile(self):
        """
        Diffusion coefficient in m^2/s.
        """
        return self.output['anomal_diffusion'][0]

    @property
    def advection_profile(self):
        """
        Advection velocity profile in m/s
        """
        return self.output['anomal_drift'][0]

    @property
    def volume_integrated_density(self):
        d = self.output

        time = d['time']
        rhovol = d['radius_grid']  # cm
        rhopol = d['rho_poloidal_grid']
        Raxis = d['large_radius']  # cm
        impurity_density = d['total_impurity_density']

        dVdr = 4 * np.pi ** 2 * Raxis * rhovol

        core = rhopol <= 1.0
        integrand = impurity_density * dVdr
        int_imp = np.trapz(integrand[:, core], rhovol[core])

        return time, int_imp

    @property
    def sxr_radiation(self):
        """
        SXR radiadion in W/cm^3. shape = (time, radial_grid)
        """
        sxr = self.output['sxr_radiation'][:, -1, :]
        return sxr

    # neoclassical transport
    def neoclassical_diffusion(self):
        keys = ['classical_diff_coeff', 'banana_plateau_diff_coeff',
                'pfirsch_schlueter_diff_coeff']

        values = dict((k, self.output[k] / 1e4) for k in keys)
        values['total'] = reduce(lambda x, y: x + y, values.values())
        return values

    def neoclassical_drift(self):
        keys = ['classical_drift', 'banana_plateau_drift',
                'pfirsch_schlueter_drift']
        values = dict((k, self.output[k] / 1e2) for k in keys)
        values['total'] = reduce(lambda x, y: x + y, values.values())
        return values

    @property
    def neoclassical_diffusion_total(self):
        return self.neoclassical_diffusion()['total']

    @property
    def neoclassical_drift_total(self):
        return self.neoclassical_drift()['total']

    def get(self, attr):
        """
        attr : str or sequence of str
        """
        keys = attr.replace(',', '').split()
        return [self._get_by_key(k) for k in keys]

    def __getitem__(self, key):
        return self.output[key]

    def _get_by_key(self, k):
        additional_keys = {
            't': 'time',
            'sxr': 'sxr_radiation',
            'nimp': 'impurity_density',
            'nimp0': 'central_impurit_density',
            'te': 'electron_temperature',
            'ne': 'electron_density',
            'Dnc': 'neoclassical_diffusion_total',
            'D': 'diffusion_profile',
            'v': 'advection_profile',
            'vnc': 'neoclassical_drift_total',
        }

        if k in additional_keys:
            return getattr(self, additional_keys[k])
        if hasattr(self, k):
            return getattr(self, k)
        elif k in self.output:
            return self.output[k]
        else:
            raise ValueError('no key: %s' % k)


class Result:
    """
    This is inteded to replace the overcomplicated Result class.
    """
    time = ['time', 'influx_through_valve']
    grid = ['rho_poloidal_grid', 'radius_grid']
    time_grid = ['electron_density', 'electron_temperature',
                 'total_impurity_density', 'anomal_diffusion', 'anomal_drift']

    kind = {'time': ['time', 'influx_through_valve'],
            'grid': ['rho_poloidal_grid', 'radius_grid'],
            'time_grid': ['electron_density', 'electron_temperature',
                          'total_impurity_density', 'anomal_diffusion',
                          'anomal_drift', 'total_impurity_density'],
            'time_radiation-stages_grid': ['sxr_radiation']}

    def __init__(self):
        self.profiles = np.rec.fromarrays([[]])

    def load(self, name):
        if isinstance(name, basestring):
            p = self._load_netcdf(netcdf_file(name, 'r'))
        else:
            p = self._load_netcdf(name)

        self.profiles = SIUnits().apply(p)
        return self

    def _load_netcdf(self, netcdf):
        a = [(n, netcdf.variables[n].data) for n in self._allnames()]
        return dict(a)

    def _apply_name_convention(self, d):
        for old, new in self.rename.iteritems():
            d[new] = d.pop(old)
        return d

    def _allnames(self):
        names = []
        for v in self.kind.values():
            names.extend(v)
        return names

    def __getitem__(self, key):
        return self.profiles[key]


class Units:
    """Consistent naming and unit convention
    names : map from string to string
        Names to replace
    units : map from string to func
        Apply function to the values with certain keys.
    """

    names = {}
    units = {}

    def apply(self, d):
        d = self._rename(d)
        d = self._convert(d)
        return d

    def _convert(self, d):
        pass

    def _rename(self, d):
        pass


icm3_to_im3 = lambda x: x * 1e6  # cm^-3 to m^-3


class SIUnits(Units):

    names = {'rho_poloidal_grid': 'rhopsi',
             'radius_grid': 'rhovol'}

    units = {'anomal_diffusion': lambda x: x / 1e4,
             'anomal_drift': lambda x: x / 1e2,
             'electron_density': icm3_to_im3,
             'total_impurity_density': icm3_to_im3,
             'sxr_radiation': icm3_to_im3} 

    def _convert(self, d):
        for name, f in self.units.iteritems():
            d[name] = f(d[name])
        return d

    def _rename(self, d):
        for old, new in self.names.iteritems():
            d[new] = d.pop(old)
        return d


def load(name):
    return Result().load(name)


# test code
from nose.tools import assert_equal, set_trace
from os.path import basename, dirname, realpath, join


def sample_output_name():
    curdir = dirname(realpath(__file__))
    return join(curdir, '..', 'sample_output.nc')


def test_sample_name():
    fname = sample_output_name()
    assert_equal(basename(fname), 'sample_output.nc')


def test_load_netdf():
    res = Result()
    d = res._load_netcdf(netcdf_file(sample_output_name(), 'r'))
    assert_equal(set(d.keys()), set(res._allnames()))


def test_load():
    res = load(sample_output_name())
    assert_equal(res.profiles['anomal_diffusion'][0, 0], 0.5)  # [m2/s]
    assert 'rhovol' in res.profiles.keys()
