__all__ = [ 'NumericalParameters', 'ImpurityParameters',
            'BackgroundParameters', 'ThomsonProfilesfromShot', 'STRAHLConfig',
            'EquilibriumfromShot', 'TestProfiles', 'CircularGeometry']

import numpy as np

class NumericalParameters(object):
    """
    >>> numpar = NumericalParameters(timeFinal=5, grid_k=10, timeStep=0.02)
    >>> p = numpar.as_dict()

    >>> assert p.has_key('numerical.grid.k')
    >>> assert p.has_key('numerical.grid.radial_points')
    >>> assert p.has_key('numerical.internal_eps')
    >>> assert p.has_key('numerical.time.dt')
    >>> assert p.has_key('numerical.time.final')
    >>> assert p.has_key('numerical.max_internal_steps')
    >>> assert p.has_key('numerical.iteration_type')
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
    >>> p = imp.as_dict()

    >>> assert p.has_key('impurity.decay_length')
    >>> assert p.has_key('impurity.delta_source')
    >>> assert p.has_key('impurity.atomic_weight')
    >>> assert p.has_key('impurity.divertor_puff')
    >>> assert p.has_key('impurity.element')
    >>> assert p.has_key('impurity.source_position')
    >>> assert p.has_key('impurity.sol_width')
    >>> assert p.has_key('impurity.energy_of_neutrals')
    >>> assert p.has_key('impurity.parallel_loss_time')

    Transport parameters
    >>> assert p.has_key('impurity.diffusion_coefficient')
    >>> assert p.has_key('impurity.convection_velocity')

    Influx time trace
    >>> assert p.has_key('impurity.influx')
    """
    
    def __init__(self, element, influx, D, v):
        if not element in ['Ar']:
            raise NotImplementedError('Unknown element: %s')

        atomic_weights = dict(Ar=18)

        self.element = element
        self.atomic_weight = atomic_weights[element]

        rho = np.linspace(0, 1, 20)
        self.D = D(rho)
        self.v = v(rho)
        self.influx = influx

    def as_dict(self):
        d = {}
        d['impurity.element'] = self.element
        d['impurity.atomic_weight'] = self.atomic_weight

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


class BackgroundParameters(object):
    pass


class TestProfiles(BackgroundParameters):
    def __init__(self):
        pass

    def as_dict(self):
        d = {}

        rho = np.linspace(0,1,20)
        ne = 1e13 * (1 - rho**2) + 1e11
        Te = 1e3 * (1 - rho**2) + 10

        d['background.atomic_weight'] = 2
        d['background.charge'] = 1
        d['background.decay_length'] = 1
        d['background.electron_temperature'] = Te
        d['background.electron_density'] = ne
        d['background.rho_poloidal'] = rho
        return d

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


class Geometry(object):
    pass


class CircularGeometry(Geometry):
    def __init__(self, radius=0.22):
        self.rho_volume_at_lcfs = radius

    def as_dict(self):
        d = {}
        d['geometry.rho_volume_at_lcfs'] = self.rho_volume_at_lcfs * 100
        d['geometry.major_radius'] = 88
        d['geometry.rho_volume'] = np.linspace(0, 1, 20)
        d['geometry.rho_volume'] *= d['geometry.rho_volume_at_lcfs']
        return d


class STRAHLConfig(object):
    def __init__(self, numerical, impurity, background, geometry):
        assert isinstance(numerical, NumericalParameters)
        assert isinstance(impurity, ImpurityParameters)
        assert isinstance(background, BackgroundParameters)
        assert isinstance(geometry, Geometry)

        self.numerical = numerical
        self.impurity = impurity
        self.background = background
        self.geometry = geometry

        self._render()

    def _render(self):
        d = {}
        d.update(self.numerical.as_dict())
        d.update(self.impurity.as_dict())
        d.update(self.background.as_dict())
        d.update(self.geometry.as_dict())
        self.params = d

        # FIXME add some global parameters
        d['index'] = 0
        d['save_all'] = True
        d['shot'] = 99999


def read_config(filename):
    loc, glob = {}, {}
    execfile(filename, loc, glob)
    return glob['config']


if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=False)
