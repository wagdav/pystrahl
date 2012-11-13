PyStrahl - Python interface for STRAHL
======================================

Code structure
--------------

strahl/__init__.py
    Module entry point

strahl/config.py
    High-level objects to generate input.

strahl/create_input.py
    Write input files.

strahl/datablocks.py
    Generate blocks of config files using the templates.

strahl/result.py
    Result object to read the output file.

strahl/misc.py
    Predefined profile and pulse shapes.

strahl/templates.py
    Templates of pieces of config files.


In quarantine
-------------
strahl/diagnostics.py
        Line integration routine.
strahl/viz.py
    Quick visualisation routines.


On the radial grids
-------------------

STRAHL uses many different grids of radial coordinates.  In the code the
following convetions are applied:

rvol
    r_vol = sqrt(volume / (2 pi^2 R(axis))

    The radius of a circular torus with a same volume V as the flux surface
    countour.  Note that this is not dimensionless, but it metres [m].

rhovol
    rho_vol = sqrt(volume / volume(LCFS)) = r_vol / r_vol(LCFS)

rhopsi, rhopol
    rho_psi = sqrt((Psi - Psi(axis)) / (Psi(LCFS) - Psi(axis)))
