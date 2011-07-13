__all__ = [ 'NumericalParameters', 'ImpurityParameters',
            'BackgroundParameters', 'ThomsonProfilesfromShot', 'STRAHLConfig',
            'EquilibriumfromShot', 'TestProfiles', 'CircularGeometry']

import numpy as np

_mandatory_keys = {
    'numerical' :
        [ 'grid.k', 'grid.radial_points', 'internal_eps', 'time.dt',
          'time.final', 'max_internal_steps', 'iteration_type' ],
    'impurity' :
        [ 'decay_length', 'delta_source', 'atomic_weight', 'divertor_puff',
          'element', 'source_position', 'sol_width', 'energy_of_neutrals',
          'parallel_loss_time', 'diffusion_coefficient',
          'convection_velocity', 'influx' ],
    'background' :
        [ 'atomic_weight', 'charge', 'decay_length', 'rho_poloidal',
          'electron_density', 'electron_temperature' ],
    'geometry' :
        ['rho_volume_at_lcfs', 'major_radius', 'rho_volume' ]
}


def _has_all_mandatory_keys(d, section):
    assert (set(d.as_dict().keys()) == _set_of_mandatory_keys(section)),\
        '%s section contains wrong keys.' % section


def _set_of_mandatory_keys(section):
    keys_needed = []
    for key in _mandatory_keys[section]:
        keys_needed.append('.'.join((section, key)))
    return set(keys_needed)


class NumericalParameters(object):
    """
    >>> numpar = NumericalParameters(timeFinal=5, grid_k=10, timeStep=0.02)
    >>> _has_all_mandatory_keys(numpar, 'numerical')
    """
    def __init__(self, timeStep, timeFinal, grid_k=10, eps=0.02,
            radial_points=101, max_internal_steps=100):
        d = {}
        d['numerical.time.final'] = timeFinal
        d['numerical.grid.k'] = grid_k
        d['numerical.time.dt'] = timeStep
        d['numerical.grid.radial_points'] = radial_points
        
        d['numerical.internal_eps'] = eps
        d['numerical.max_internal_steps'] = max_internal_steps
        d['numerical.iteration_type'] = 1
        self._d = d

    def as_dict(self):
        return self._d


class ImpurityParameters(object):
    """
    >>> from strahl.misc import *

    >>> influx = rectangular_pulse(5e-3, 5e21)
    >>> D = constant(0.1)
    >>> v = constant(-1)
    >>> imp = ImpurityParameters(element='Ar', influx=influx, D=D, v=v)
    >>> _has_all_mandatory_keys(imp, 'impurity')
    """
    def __init__(self, element, influx, D, v):
        rho = np.linspace(0, 1, 20)
        self.D = D(rho)
        self.v = v(rho)
        self.influx = influx
        self.element = element

    def get_influx(self):
        return self._influx

    def set_influx(self, influx):
        t, flx = influx
        t, flx = np.asarray(t), np.asarray(flx)
        assert t.shape == flx.shape, 'Invalid shape t.shape != flx.shape'
        self._influx = (t, flx)

    def get_element(self):
        return self._element, self._atomic_weight

    def set_element(self, element):
        atomic_weights = dict(Ar=18)

        if not element in atomic_weights.keys():
            raise NotImplementedError('%s as impurity is not supported.'
                    % element)
        self._element = element
        self._atomic_weight = atomic_weights[element]

    def as_dict(self):
        d = {}
        element, atomic_weight = self.element
        d['impurity.element'] = element
        d['impurity.atomic_weight'] = atomic_weight

        # transport coefficients
        d['impurity.diffusion_coefficient'] = self.D
        d['impurity.convection_velocity'] = self.v

        # impurity influx
        d['impurity.influx'] = self.influx

        # FIXME some hardwired values
        d['impurity.energy_of_neutrals'] = 1.0
        d['impurity.parallel_loss_time'] = 2.5
        d['impurity.sol_width'] = 1.0
        d['impurity.source_position'] = 1000
        d['impurity.divertor_puff'] = False
        d['impurity.delta_source'] = 0
        d['impurity.decay_length'] = 0.5

        return d

    influx = property(get_influx, set_influx)
    element = property(get_element, set_element)


class BackgroundParameters(object):
    def __init__(self, element, rho, ne, Te):
        self.element = element
        self.profiles = (rho, ne, Te)

    def set_element(self, element):
        if element not in ['H', 'D']:
            raise NotImplementedError('%s as background plasma is not'
                ' supported.' % element)
        self._element = element

    def get_element(self):
        return self._element

    def get_profiles(self):
        """
        Retrun (rho, ne, Te) tuple.
        """
        return self._profiles

    def set_profiles(self, profiles):
        rho, ne, Te = profiles
        assert rho.shape == ne.shape, 'Incorrect shape rho.shape != ne.shape'
        assert rho.shape == Te.shape, 'Incorrect shape rho.shape != Te.shape'
        self._profiles = (rho, ne, Te)

    def as_dict(self):
        atomic_weights = dict(H=1, D=2)
        charges = dict(H=1, D=1)
        rho, ne, Te = self.profiles

        d = {}
        d['background.atomic_weight'] = atomic_weights[self.element]
        d['background.charge'] = charges[self.element]
        d['background.decay_length'] = 1 # FIXME: implicit parameter
        d['background.rho_poloidal'] = rho
        d['background.electron_density'] = ne
        d['background.electron_temperature'] = Te

        return d

    element = property(get_element, set_element)
    profiles = property(get_profiles, set_profiles)


class TestProfiles(BackgroundParameters):
    """
    Deuterium plasma background with parabolic temperature and density
    profiles.

    >>> bg = TestProfiles()
    >>> _has_all_mandatory_keys(bg, 'background')
    >>> bg.element
    'D'
    >>> rho, ne, Te = bg.profiles
    >>> rho.min(), rho.max(), rho.mean()
    (0.0, 1.0, 0.5)
    """
    def __init__(self):
        rho = np.linspace(0,1,20)
        ne = 1e13 * (1 - rho**2) + 1e11
        Te = 1e3 * (1 - rho**2) + 10
        super(TestProfiles, self).__init__('D', rho, ne, Te)


class ThomsonProfilesfromShot(BackgroundParameters):
    def __init__(self, shot, time_bbox):
        self.shot = shot
        self.time_bbox = time_bbox

    def as_dict(self):
        """MDS Acess downloading data, etc"""
        return {}


class EquilibriumfromShot(object):
    def __init__(self, shot, time_bbox):
        pass
    pass


class GeometryParameters(object):
    pass


class CircularGeometry(GeometryParameters):
    """
    Circular plasma with a given minor and major radius (both in meters).

    >>> geom = CircularGeometry(minor_radius=0.22, major_radius=0.88)
    >>> _has_all_mandatory_keys(geom, 'geometry')
    """
    def __init__(self, minor_radius, major_radius):
        self._minor_radius = minor_radius
        self._major_radius = major_radius

    def as_dict(self):
        d = {}
        d['geometry.rho_volume_at_lcfs'] = self._minor_radius * 1e2
        d['geometry.major_radius'] = self._major_radius * 1e2
        d['geometry.rho_volume'] = np.linspace(0, 1, 20)
        d['geometry.rho_volume'] *= d['geometry.rho_volume_at_lcfs']
        return d


class STRAHLConfig(object):
    def __init__(self, numerical, impurity, background, geometry):
        assert isinstance(numerical, NumericalParameters)
        assert isinstance(impurity, ImpurityParameters)
        assert isinstance(background, BackgroundParameters)
        assert isinstance(geometry, GeometryParameters)

        self.numerical = numerical
        self.impurity = impurity
        self.background = background
        self.geometry = geometry

    def get_params(self):
        d = {}
        d.update(self.numerical.as_dict())
        d.update(self.impurity.as_dict())
        d.update(self.background.as_dict())
        d.update(self.geometry.as_dict())

        # FIXME add some global parameters
        d['index'] = 0
        d['save_all'] = True
        d['shot'] = 99999

        return d

    params = property(get_params)


def read_config(filename):
    loc, glob = {}, {}
    execfile(filename, loc, glob)
    return glob['config']


if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=False)
