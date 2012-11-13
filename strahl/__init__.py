__all__ = []

from input import create_input, run
from misc import *
from result import load


def default_params():
    from numpy import linspace, zeros_like
    _r = linspace(0, 1, 20)

    return {
        # geometry
        'geometry.rvol_lcfs': 20.0,
        'geometry.major_radius': 88.0,
        'geometry.rhovol': _r * 0.22,
        'geometry.rhopol': _r,
        'geometry.sol_width': 8.0,
        'geometry.limiter_position': 6.5,

        # plasma background
        'background.rhopol': _r,
        'background.electron_density': 1e16 * (1 - _r ** 2),
        'background.electron_temperature': 1e3 * (1 - _r ** 2),
        'background.atomic_weight': 2,
        'background.charge': 1,
        'background.decay_length': 1,

        # impurity
        'impurity.element': 'Ar',
        'impurity.atomic_weight': 18,
        'impurity.energy_of_neutrals': 1.0,
        'impurity.rhopol': _r,
        'impurity.diffusion_coefficient': _r ** 2 + 0.01,
        'impurity.convection_velocity': zeros_like(_r),
        'impurity.influx': ([1.0], [2.51e21]),
        'impurity.source_position': 1000,

        'impurity.divertor_puff': False,
        'impurity.delta_source': 0,
        'impurity.sol_width': 1.0,
        'impurity.decay_length': .5,
        'impurity.parallel_loss_time': 2.5,

        # recycling
        'recycling.switch': False,
        'recycling.wall_R': 1,
        'recycling.tau_divsol': 1e9,
        'recycling.tau_pump': 1e9,

        # neoclassical transport
        'neoclassical.on': False,
        'neoclassical.type': 3,
        'neoclassical.no_bp_contrib': False,
        'neoclassical.high_limit': 1.,
        'neoclassical.low_limit': 0,

        # numerical parameters
        'numerical.grid.k': 2,
        'numerical.grid.radial_points': 101,
        'numerical.max_internal_steps': 100,
        'numerical.internal_eps': 0.02,
        'numerical.iteration_type': 1,
        'numerical.time.final': 1.0,
        'numerical.time.dt': 1e-3,

        # globals
        'shot': 99999,
        'index': 0,
        'save_all': True
    }


defaultParams = default_params


def sanity_check(params):
    a = set(default_params().keys())
    b = set(params.keys())
    assert a == b, \
        'Parameter set contains an invalid parameter name: %s' %\
        (b - a,)
