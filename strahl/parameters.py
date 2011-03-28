import numpy as np

# Input parameter names as listed in *read_parameter.f*
_description = dict(
    #element
    element='element',
    aimp='atomic weight',
    en0='energy of neutrals (eV)',
    apl='background ion atomic weight',
    zp='background ion atomic charge',
    # grid - file
    shot_grid='shot',
    index_grid='index',
    k='rho = r**K (->K)',
    ir='number of grid points',
    max_internal_steps='max. iterations at fixed time',
    inteps='stop iteration if change below (%%)',
    itmodus='1=iteration, 2=no iteration',
    # start conditions
    startfrom='start new=0/from old impurity distribution=1',
    shot_old='take distr. from shot',
    time_old='taken distr. time',
    # output
    saveall='save all cycles = 1, save final and start distribution=0',
    # timesteps
    n_change='number of changes (start-time+...+stop-time',
    timestep_changes=
        'time dt-at-start increase-of-dt-after-cycle steps per cycle',
    # source
    rl='position (cm)',
    flx_t='constant rate (1/s)',
    flxfromfile='time dependent rate from file (1/0)',
    divbls='divertor puff',
    delta_source='delta_source',
    # edge, recycling
    fall_outsol='decay lenght of impurity outside last grid point (cm)',
    rclswitch='recycling ON=1/OFF=0',
    rcl='Wall-rec.',
    taudiv='tau-div-> SOL9ms)',
    taupump='taup-pump (ms)',
    db='SOL-width (cm)',
    # density, temperature and neutral hydrogen for cx
    shot='take from file with shot',
    index='take from file with index',
    # neoclassical transport
    #TODO
)

_r = np.linspace(0,1,41)

_defaultParams = {
    # geometry
    'geometry.vol_lcfs' : 61.6,
    'geometry.r_maj' : 88.0,
    'geometry.time' : 1.0,
    'geometry.n_grid' : 11,
    'geometry.n_sep' : 9,
    'geometry.n_fourier' : 1,
    'geometry.rho_pol' : np.hstack((np.linspace(0,1,9), [1.02, 1.04])),
    'geometry.rho_vol' : np.hstack((np.linspace(0,1,9), [1.02, 1.04])),
    'geometry.R_lfs' : np.linspace(1, 1.1, 11),
    'geometry.R_hfs' : np.linspace(1, 1.1, 11),

    # plasma background
    'background.rho_poloidal' : _r,
    'background.electron_density' : 1e16 * (1 - _r**2),
    'background.electron_temperature' : 1e3 * (1 - _r**2),
    'background.atomic_weight' : 2,
    'background.charge' : 1,

    # impurity
    'impurity.element' : 'Ar',
    'impurity.atomic_weight' : 18,
    'impurity.energy_of_neutrals' : 1.0,
    'impurity.diffusion_coefficient' : _r**2 + 0.01,
    'impurity.convection_velocity' : np.zeros_like(_r),
    'impurity.influx' : ([1.0], [2.51e21]),
    'impurity.source_position' : 1000,

    'impurity.divertor_puff' : False,
    'impurity.delta_source' : 0,
    'impurity.sol_width' : 5.0,
    'impurity.decay_length' : 2.0,
    'impurity.parallel_loss_time' : 2.5,

    # numerical parameters
    'numerical.grid.k' : 2,
    'numerical.grid.radial_points' : 101,
    'numerical.max_internal_steps' : 100,
    'numerical.internal_eps' : 0.02,
    'numerical.iteration_type' : 1,
    'numerical.time.final' : 1.0,
    'numerical.time.dt' : 1e-3,

    # globals
    'shot' : 99999,
    'index' : 0,
    'save_all': True,
}


def defaultParams():
    return _defaultParams.copy()
