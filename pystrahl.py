
verbose=False
control_file='strahl.control'
create_rad=True
sxr_dat_only=False
diag_only=False
rate_out=False
inp_mode='A'


#import pystrahl
#pystrahl.strahl(verbose,create_rad,inp_mode, sxr_dat_only, diag_only,
#        rate_out,control_file)

strahl_param_file=\
r"""
ELEMENT
=======

cv     element   atomic weight  energy of neutrals(eV)
  '%(element)s'       %(aimp)f             %(en0)f

cv    background ion:  atomic weight    charge
  %(apl)f %(zp)f     

GRID - FILE
===========

cv    shot      index
  %(shot_grid)s      %(index_grid)d 

GRID POINTS AND ITERATION
=========================

cv     rho = r**K (->K)      number of grid points
  %(k)f           %(ir)d                    

max. iterations at fixed time  stop iteration if change below(%%)
      1=iteration,2=no iteration
cv     ion/rec-step
       %(max_internal_steps)d       %(inteps)f %(itmodus)d 

START CONDITIONS
================

cv start new=0/from old impurity distribution=1  take distr. from shot   at    time
    %(startfrom)d       %(shot_old)d	%(time_old)f 

OUTPUT
======

cv    save all cycles = 1, save final and start distribution = 0
	%(saveall)d

TIMESTEPS
=========

cv    number of changes (start-time+... +stop-time) 
  %(n_change)d

cv    time    dt at start    increase of dt after cycle     steps per cycle
  %(timestep_changes)s

SOURCE
======
   
cv    position(cm)    constant rate(1/s)   time dependent rate from file(1/0)
  %(rl)f %(flx_t)e %(flxfromfile)d

cv    divertor puff   delta_source
  %(divbls)f   %(delta_source)f

EDGE, RECYCLING
===============

cv    decay length of impurity outside last grid point(cm)
  %(fall_outsol)f

cv    Rec.:ON=1/OFF=0   wall-rec.  Tau-div->SOL(ms)   Tau-pump(ms) 
  %(rclswitch)d  %(rcl)f %(taudiv)f %(taupump)f

cv    SOL-width(cm)     
  %(db)f

DENSITY, TEMPERATURE  AND NEUTRAL  HYDROGEN  FOR  CX
====================================================
 
cv    take from file with:    shot        index
  %(shot)s %(index)d


NEOCLASSICAL TRANSPORT 
======================

method 0 = off,  >0 = %% of Drift,    1 = approx.
cv  <0 =figure out, but dont use   2/3 = NEOART   neoclassics for rho_pol <
  %(qf)d  %(neo_method)d  %(rho_pol_neo)f

ANOMALOUS  TRANSPORT
====================
 
cv    # of changes for transport 
  %(nt_trans)d
 
cv    time-vector  
  %(t_trans)s
 
cv    parallel loss times (ms) 
  %(tau_t)s              
 
cv      Diffusion  [m^2/s]
  '%(fprofil)s'

cv   	D	D__SOL	Dmax	beta	gamma	r1/r(sol)   
  0.1	1.0     2.5      15.     6.0     0.8

cv    Drift function       Drift Parameter/Velocity
  'const_c'                 'velocity'

cv num of sawteeth   inversion radius (cm)
  %(n_saw)d  %(r_saw)f

cv    times of sawteeth 
  %(t_saw)f
"""

# Input parameter names as listed in *read_parameter.f*
params_descr=dict(
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
    
)


params_defaults=dict(
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


def casedir_init (casedir=None):
    import os, errno

    if casedir:
        os.mkdir(casedir)
    else: 
        casedir = os.getcwd()

    directories=['nete', 'param_files', 'result']

    for d in directories:
        os.mkdir(os.path.join(casedir, d))
# def casedir_init


# Find all the parts which are to be interpolated
#import re
#params_used=re.findall(r'\%[(](.*?)[)][esfd]',strahl_param_file)

#print strahl_param_file%params_defaults
