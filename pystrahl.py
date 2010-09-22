
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
    E L E M E N T

cv     element   atomic weight  energy of neutrals(eV)
    %(element)s       %(aimp)f             %(en0)f

cv    background ion:  atomic weight    charge
                      %(apl)f %(zp)f     

			     G R I D - F I L E
cv    shot      index
%(shot_grid)s      %(index_grid)d 

    G R I D   P O I N T S  A N D  I T E R A T I O N

cv     rho = r**K (->K)      number of grid points
             %(k)f           %(ir)d                    

      max. iterations at fixed time  stop iteration if change below(%%)
      1=iteration,2=no iteration
cv     ion/rec-step
       %(max_internal_steps)d       %(inteps)f %(itmodus)d 

    S T A R T   C O N D I T I O N S

cv    start new=0/from old impurity distribution=1     take distr. from shot   at    time
    %(startfrom)d       %(shot_old)d	%(time_old)f 

    O U T P U T   

cv    save all cycles = 1, save final and start distribution = 0
	%(saveall)d

    T I M E S T E P S     

cv    number of changes (start-time+... +stop-time) 
    %(n_change)d

cv    time    dt at start    increase of dt after cycle     steps per cycle
      0.0        1.e-6               1.20                         10
      0.01        2.e-5               1.05                         5
                             
  
			       S O  U R C E  
   
cv    position(cm)    constant rate(1/s)   time dependent rate from file(1/0)
 %(rl)f         %(flx_t)f   %(flxfromfile)d

cv    divertor puff   delta_source
 %(divbls)f   %(delta_source)f

                    E D G E ,  R E C Y C L I N G

cv    decay length of impurity outside last grid point(cm)
    %(fall_outsol)f

cv    Rec.:ON=1/OFF=0   wall-rec.  Tau-div->SOL(ms)   Tau-pump(ms) 
 %(rclswitch)d  %(rcl)f %(taudiv) %(taupump)

cv    SOL-width(cm)     
 %(db)f


  DENSITY, TEMPERATURE  AND NEUTRAL  HYDROGEN  FOR  CX 
 
cv    take from file with:    shot        index
 %(shot)s    %(index)d


  NEOCLASSICAL TRANSPORT 
                   method
     0 = off,  >0 = %% of Drift,    1 = approx.
cv  <0 =figure out, but dont use   2/3 = NEOART   neoclassics for rho_pol <
  %(qf)d  %(neo_method)d  %(rho_pol_neo)f

  ANOMALOUS  TRANSPORT 
 
cv    # of changes for transport 
    %(nt_trans)d
 
cv    time-vector  
    %(t_trans)s
 
cv    parallel loss times (ms) 
    %(tau_t)s              
 
cv      Diffusion  [m^2/s]
    %(fprofil)s

cv   	D	D__SOL	Dmax	beta	gamma	r1/r(sol)   
        0.1	1.0     2.5      15.     6.0     0.8

cv    Drift function       Drift Parameter/Velocity
      'const_c'                 'velocity'

cv    Drift function       Drift Parameter/Velocity
      'interp'                 'velocity'   
 
cv   # of interpolation points
     %(vd_nip)s
 
cv   rho poloidal grid for interpolation
    %(vd_xip)s
 
cv    Drift Parameter/Velocity
    %(vd_yip)s
 
 
cv    # of sawteeth     inversion radius (cm) 
    %(n_saw)f %(r_saw)f
 
 
cv    times of sawteeth 
    %(t_saw)f
"""

params_descr=dict(
    element='element',
    aimp='atomic weight',
    en0='energy of neutrals (eV)',
    apl='background ion atomic weight',
    zp='background ion atomic charge',
    shot_grid='shot',
    index_grid='index',
    )

config_dict=dict(

    element='C_',
    aimp=12.,
    en0=1.0,

    background_ion_atomic_weight=2.0,
    background_ion_charge=1.0,

    shot_index='30073',
    num_grid_points=601,

    num_maxiter=100,
    rtol=0.0002,

#    rho_grid = [0.0, 1.0],

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
 ir                 = 601,
 max_internal_steps = 100,
 inteps             = 0.0002,
 itmodus            = 1,
 startfrom          = 0,
 shot_old           = 30073,
 time_old           = 1.0,
 saveall            = 1,
 n_change           = 2,
 # t_change, dtstart, itz, tinc # time, dt at start, increase of dt, steps
 # cycle
 rl                 = 10,
 flx_t              = 2.5e10,
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
 t_trans            = '',
 tau_t              = '',
 fprofil            = 0.0,
 vd_nip             = '',
 vd_xip             = '',
 vd_yip             = '',
 n_saw              = 0,
 r_saw              = 35,
 t_saw              = 0.0,
)


# Input parameter names as listed in *read_parameter.f*
params_all = [
        't01', 'startfrom', 'k', 'ir', 'therm_cx', 'diag_lin_on', 'rl',
        'n_change', 't_change', 'dtstart', 'saveall', 'aimp', 'apl', 'zp',
        'r', 'rmaj', 'itmodus', 'max_internal_steps', 'ineps', 'itz', 'tinc',
        'db', 'en0', 'v0', 'qf', 'neo_method', 'rho_pol_neo', 'delta_source',
        'taudiv', 'taupump', 'rcl', 'rclswitch', 'r_saw', 'n_saw', 't_saw',
        'element', 'shot', 'shot_grid', 'index_grid', 'shot_old', 'time_old'
        
        ]

# Find all the parts which are to be interpolated
import re
params_used=re.findall(r'\%[(](.*?)[)][sfd]',strahl_param_file)

not_used = set(params_all) - set(params_used)
unknown = set(params_used) - set(params_all)




print strahl_param_file%params_defaults
