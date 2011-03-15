# Input parameter names as listed in *read_parameter.f*
description = dict(
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


defaults = dict(
    element            = 'C_',
    aimp               = 12,
    en0                = 1.0,
    apl                = 2,
    zp                 = 1,
    shot_grid          = 30073,
    index_grid         = 0,
    k                  = 2.0,
    ir                 = 101,
    max_internal_steps = 100,
    inteps             = 0.02,
    itmodus            = 1,
    startfrom          = 0,
    shot_old           = 30073,
    time_old           = 1.0,
    saveall            = 1,
    n_change           = 1,
    timestep_changes   = '0.0 1.e-3 1.20 10',
    # SOURCE
    rl                 = 1000,
    flx_t              = 2.5e21,
    flxfromfile        = 0,
    divbls             = 0,
    delta_source       = 0,

    # EDGE RECYCLING
    fall_outsol        = 2.0,
    rclswitch          = 0,
    rcl                = 1,
    taudiv             = 55,
    taupump            = 240,
    db                 = 5.0,

    # DENSITY, TEMPERATURE AND NEUTRAL FOR CX
    shot               = 30073,
    index              = 0,

    # NEOCLASSICAL TRANSPORT
    qf                 = 0,
    neo_method         = 2,
    rho_pol_neo        = .95,

    # ANOMALOUS TRANSPORT
    nt_trans           = 1,
    t_trans            = 0.0,
    tau_t              = 2.5,
    fprofil            = "funct",
    vd_nip             = '',
    vd_xip             = '',
    vd_yip             = '',
    n_saw              = 0,
    r_saw              = 35,
    t_saw              = 0.0,
)

