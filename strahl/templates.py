"""
Template strings for generating input files for STRAHL runs.
"""


param_file=\
r"""ELEMENT
cv     element   atomic weight  energy of neutrals(eV)
 '%(impurity.element)s' %(impurity.atomic_weight)f %(impurity.energy_of_neutrals)f

cv    background ion:  atomic weight    charge
 %(background.atomic_weight)f %(background.charge)f


GRID - FILE
cv    shot      index
  %(shot)s      %(index)d


GRID POINTS AND ITERATION
cv     rho = r**K (->K)      number of grid points
 %(numerical.grid.k)f %(numerical.grid.radial_points)d

   max. iterations     stop iteration       1=iteration,2=no iteration
cv at fixed time      if change below(%%)  ion/rec-step
 %(numerical.max_internal_steps)d  %(numerical.internal_eps)f %(numerical.iteration_type)d


START CONDITIONS
cv start new=0/from old impurity distribution=1  take distr. from shot   at    time
    0    -1   -1


OUTPUT
cv    save all cycles = 1, save final and start distribution = 0
 %(save_all)d


TIMESTEPS
cv    number of changes (start-time+... +stop-time)
 2

cv    time    dt at start    increase of dt after cycle     steps per cycle
 0 %(numerical.time.dt)1.1e 1.0 10
 %(numerical.time.final)1.1f          -1 -1 -1


SOURCE
cv    position(cm)    constant rate(1/s)   time dependent rate from file(1/0)
 %(impurity.source_position)f 0 1

cv    divertor puff   delta_source
 %(impurity.divertor_puff)d   %(impurity.delta_source)d


EDGE, RECYCLING
cv    decay length of impurity outside last grid point(cm)
 %(impurity.decay_length)f

cv    Rec.:ON=1/OFF=0   wall-rec.  Tau-div->SOL(ms)   Tau-pump(ms)
 0 -1 -1 -1

cv    SOL-width(cm)
 %(impurity.sol_width)f

DENSITY, TEMPERATURE  AND NEUTRAL  HYDROGEN  FOR  CX
cv    take from file with:    shot        index
 %(shot)s %(index)d


NEOCLASSICAL TRANSPORT
method 0 = off,  >0 = %% of Drift,    1 = approx.
cv  <0 =figure out, but dont use   2/3 = NEOART   neoclassics for rho_pol <
  0 -1 -1


ANOMALOUS  TRANSPORT
cv    # of changes for transport
 1

cv    time-vector
 0.0

cv    parallel loss times (ms)
  %(impurity.parallel_loss_time)s

%(transport_datablock)s

cv num of sawteeth   inversion radius (cm)
  0 -1

cv    times of sawteeth
  -1
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
%(geometry.rho_volume_at_lcfs)4.1f %(geometry.major_radius)4.1f

cv number of grid points, points up to LCFS, Fourier coeffs
%(n_grid)d %(n_sep)d 1

cv rho poloidal
%(background.rho_poloidal)s

cv rho volume / rho_volume (LCFS)
%(geometry.rho_volume)s

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


