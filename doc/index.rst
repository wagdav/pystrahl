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
    
