"""
Template strings for generating input files for STRAHL runs.

Test if the the *param_file* string can be interpolated with the default
parameters:

>>> import defaults
>>> config_file =  param_file % defaults.main
"""


param_file=\
r"""ELEMENT
cv     element   atomic weight  energy of neutrals(eV)
  '%(element)s'       %(aimp)f             %(en0)f

cv    background ion:  atomic weight    charge
  %(apl)f %(zp)f


GRID - FILE
cv    shot      index
  %(shot)s      %(index)d


GRID POINTS AND ITERATION
cv     rho = r**K (->K)      number of grid points
  %(k)f           %(ir)d

   max. iterations     stop iteration       1=iteration,2=no iteration
cv at fixed time      if change below(%%)  ion/rec-step
   %(max_internal_steps)d  %(inteps)f       %(itmodus)d


START CONDITIONS
cv start new=0/from old impurity distribution=1  take distr. from shot   at    time
    %(startfrom)d       %(shot_old)d	%(time_old)f


OUTPUT
cv    save all cycles = 1, save final and start distribution = 0
	%(saveall)d


TIMESTEPS
cv    number of changes (start-time+... +stop-time)
  %(n_change)d

cv    time    dt at start    increase of dt after cycle     steps per cycle
  0   %(dt)1.1e 1.0 10
  %(t_final)1.1f          -1 -1 -1


SOURCE
cv    position(cm)    constant rate(1/s)   time dependent rate from file(1/0)
  %(rl)f 0 1

cv    divertor puff   delta_source
  %(divbls)d   %(delta_source)d


EDGE, RECYCLING
cv    decay length of impurity outside last grid point(cm)
  %(fall_outsol)f

cv    Rec.:ON=1/OFF=0   wall-rec.  Tau-div->SOL(ms)   Tau-pump(ms)
  %(rclswitch)d  %(rcl)f %(taudiv)f %(taupump)f

cv    SOL-width(cm)
  %(db)f

DENSITY, TEMPERATURE  AND NEUTRAL  HYDROGEN  FOR  CX
cv    take from file with:    shot        index
  %(shot)s %(index)d


NEOCLASSICAL TRANSPORT
method 0 = off,  >0 = %% of Drift,    1 = approx.
cv  <0 =figure out, but dont use   2/3 = NEOART   neoclassics for rho_pol <
  %(qf)d  %(neo_method)d  %(rho_pol_neo)f


ANOMALOUS  TRANSPORT
cv    # of changes for transport
  %(nt_trans)d

cv    time-vector
  %(t_trans)s

cv    parallel loss times (ms)
  %(tau_t)s

%(transport_datablock)s

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
"""cv rho volume(LCFS)[cm] R_axis [cm] time[s]
%(vol_lcfs)4.1f %(r_maj)4.1f %(time)4.1f

cv number of grid points, points up to LCFS, Fourier coeffs
%(n_grid)d %(n_sep)d %(n_fourier)d

cv rho poloidal
%(rho_pol)s

cv rho volume / rho_volume (LCFS)
%(rho_vol)s

cv major radius low field side / R_axis
%(R_lfs)s

cv major radius high field side / R_axis
%(R_hfs)s
"""


transport_datablock="""cv function
 'interp' %(drift)s

cv # number of interpolation points
 %(n_points)d

cv rho poloidal grid
 %(x)s

cv D[m**2/s] / v[m/s]
 %(y)s
"""


def used_parameters():
    import re
    return re.findall(r'\%[(](.*?)[)][esfd]',param_file)


if __name__ == '__main__':
    import doctest
    doctest.testmod()


