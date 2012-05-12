__all__ = [ 'NumericalParameters', 'ImpurityParameters',
            'BackgroundParameters', 'ThomsonProfilesfromShot', 'STRAHLConfig',
            'EquilibriumfromShot', 'TestProfiles', 'CircularGeometry']

import numpy as np
import periodictable

_mandatory_keys = {
    'numerical' :
        [ 'grid.k', 'grid.radial_points', 'internal_eps', 'time.dt',
          'time.final', 'max_internal_steps', 'iteration_type' ],
    'impurity' :
        [ 'decay_length', 'delta_source', 'atomic_weight', 'divertor_puff',
          'element', 'source_position', 'sol_width', 'energy_of_neutrals',
          'parallel_loss_time', 'diffusion_coefficient',
          'convection_velocity', 'influx', 'rho_poloidal'],
    'background' :
        [ 'atomic_weight', 'charge', 'decay_length', 'rho_poloidal',
          'electron_density', 'electron_temperature' ],
    'geometry' :
        ['rho_volume_at_lcfs', 'major_radius', 'rho_volume', 'sol_width',
            'limiter_position' ]
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
    >>> from misc import *

    >>> influx = rectangular_pulse(5e-3, 5e21)
    >>> D = constant(0.1)
    >>> v = constant(-1)
    >>> imp = ImpurityParameters(element='Ar', influx=influx, D=D, v=v)
    >>> _has_all_mandatory_keys(imp, 'impurity')
    """
    def __init__(self, element, influx, D, v, energy_of_neutrals=1.):
        self.rho = np.linspace(0, 1, 20)
        self.D = D(self.rho)
        self.v = v(self.rho)
        self.influx = influx
        self.element = element
        self.energy_of_neutrals = 1.

        self.set_recycling(R=0)

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
        isotope = periodictable.elements.isotope(element)
        self._element = isotope.symbol
        self._atomic_weight = isotope.mass

    def set_recycling(self, R, divsol=1e9, pump=1e9):
        self.recycling_switch = bool(R)
        self.recycling_R = R
        self.recycling_divsol = divsol
        self.recycling_pump = pump

    def set_transport_profiles(self, rho, D, v):
        self.rho = rho
        self.D = D
        self.v = v

    def as_dict(self):
        d = {}
        element, atomic_weight = self.element
        d['impurity.element'] = self._pad_with_underscore(element)
        d['impurity.atomic_weight'] = atomic_weight

        # transport coefficients
        d['impurity.diffusion_coefficient'] = self.D
        d['impurity.convection_velocity'] = self.v
        d['impurity.rho_poloidal'] = self.rho

        # impurity influx
        d['impurity.influx'] = self.influx

        d['impurity.energy_of_neutrals'] = self.energy_of_neutrals

        # recyling
        d['recycling.switch'] = self.recycling_switch
        d['recycling.wall_R'] = self.recycling_R
        d['recycling.tau_divsol'] = self.recycling_divsol
        d['recycling.tau_pump'] = self.recycling_pump

        # FIXME some hardwired values
        d['impurity.parallel_loss_time'] = 2.5
        d['impurity.sol_width'] = 1.0
        d['impurity.source_position'] = 6
        d['impurity.divertor_puff'] = False
        d['impurity.delta_source'] = 0
        #d['impurity.decay_length'] = 1000
        d['impurity.decay_length'] = 0.5

        return d

    def _pad_with_underscore(self, name):
        if len(name) == 2:
            ret = name
        elif len(name) == 1:
            ret = name + '_'
        else:
            raise ValueError('invalid element symbol %s' % name)
        return ret

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
        rho = np.linspace(0, 1.0, 20)
        r = rho / rho.max()
        ne = 3e13 * (1 - r**2) + 0.5e13
        Te = 1e3 * (1 - r**2) + 150
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
        d['geometry.sol_width'] = 1
        d['geometry.limiter_position'] = 1
        return d


class NeoclassicalTransport(object):
    def __init__(self, type_='off', limits=(0., 1.)):
        self.eta = 0
        self.set_type(type_)
        self.limits = limits

    def set_type(self, p, eta=100):
        type_ = 3
        if p in ['one stage']:
            type_ = 2
        elif p in ['all stages', 'on']:
            type_ = 3
        elif p in ['off']:
            eta = 0
        else:
            raise ValueError('unknown method for neoclassical transport: %s' %
                    p)

        self.type_ = type_
        self.eta = eta

    def as_dict(self):
        d = {}
        d['neoclassical.on'] = self.eta
        d['neoclassical.type'] = self.type_
        d['neoclassical.no_bp_contrib'] = 0
        d['neoclassical.high_limit'] = self.limits[1]
        d['neoclassical.low_limit'] = self.limits[0]

        return d


class STRAHLConfig(object):
    def __init__(self, numerical, impurity, background, geometry,
            neoclassical=None):
        assert isinstance(numerical, NumericalParameters)
        assert isinstance(impurity, ImpurityParameters)
        assert isinstance(background, BackgroundParameters)
        assert isinstance(geometry, GeometryParameters)

        self.numerical = numerical
        self.impurity = impurity
        self.background = background
        self.geometry = geometry

        if neoclassical == None:
            self.neoclassical = NeoclassicalTransport()
        else:
            self.neoclassical = neoclassical

    def get_params(self):
        d = {}
        d.update(self.numerical.as_dict())
        d.update(self.impurity.as_dict())
        d.update(self.background.as_dict())
        d.update(self.geometry.as_dict())
        d.update(self.neoclassical.as_dict())

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
