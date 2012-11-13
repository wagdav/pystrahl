import numpy as np

from matplotlib.path import Path
from matplotlib.patches import PathPatch

import crpppy.mds as mds

_nodes = {
    'R.in': 'static("r_v:in")',
    'R.out': 'static("r_v:out")',
    'Z.in': 'static("z_v:in")',
    'Z.out': 'static("z_v:out")',
}

vessel = {}

c = mds.TCV_Connection(-1)
for key, value in _nodes.iteritems():
    vessel[key] = c.tdi(value).data


def vessel_patch():
    vertices_in = [r for r in zip(vessel['R.in'], vessel['Z.in'])]
    vertices_out = [r for r in zip(vessel['R.out'], vessel['Z.out'])][::-1]

    vertices_in.append(vertices_in[0])
    vertices_out.append(vertices_out[0])

    codes_in = [Path.MOVETO] + (len(vertices_in)-1)* [Path.LINETO]
    codes_out = [Path.MOVETO] + (len(vertices_out)-1)* [Path.LINETO]

    vessel_path = Path(vertices_in + vertices_out, codes_in + codes_out)
    vessel_patch = PathPatch(vessel_path, facecolor=(0.6, 0.6, 0.6),
                          edgecolor='black')
    return vessel_patch


def dmv_patch():
    vertices = [(1.16, 0.455),(0.6, 0.455)]
    codes = [Path.MOVETO] + (len(vertices)-1)* [Path.LINETO]
    dmv_path = Path(vertices, codes)
    dmv_patch = PathPatch(dmv_path, edgecolor='black')

    return dmv_patch



