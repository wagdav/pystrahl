"""
Template strings for generating input files for STRAHL runs.

Test if the the *param_file* string can be interpolated with the default
parameters:

>>> import parameters
>>> config_file =  param_file % parameters.defaults
"""


param_file=\
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

cv  D D_SOL Dmax beta gamma r1/r(sol)
  0.1	1.0     2.5      15.     6.0     0.8

cv    Drift function       Drift Parameter/Velocity
  'const_c'                 'velocity'

cv num of sawteeth   inversion radius (cm)
  %(n_saw)d  %(r_saw)f

cv    times of sawteeth
  %(t_saw)f
"""


pp_datablock=\
r"""
cv time vector
  %(ntimes)d
  %(time_vector)s

cv ne function
  'interp'

cv radial coordinate
  'poloidal rho'

cv # of interpolation points
  %(npoints)d

cv radial grid for interpolation
  %(rgrid)s

cv ne[cm^-3]/ Te[eV] / Ti [eV]
  %(ygrid)s

cv decay length[cm] in rho_volume
  %(decay_length)s
"""


geometry =\
"""
cv rho volume(LCFS)[cm] R_axis [cm] time[s]
%(vol_lcfs)f %(r_maj)f

cv number of grid points, points up to LCFS, Fourier coeffs
%(n_grid)d %(n_sep)d
"""


def used_parameters():
    import re
    return re.findall(r'\%[(](.*?)[)][esfd]',param_file)


if __name__ == '__main__':
    import doctest
    doctest.testmod()


