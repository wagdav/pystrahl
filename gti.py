import os

import h5py
import numpy as np

def inverted_data(shot):
    cache_file = temp_filename(shot)
    gti = h5py.File(cache_file, 'r')

    inverted = gti['res/inverted']
    rho = gti['res/rho'][:,0]
    time = gti['res/time'][0,:]

    inverted = np.array(inverted)
    return (rho, time, inverted)


def temp_filename(shot):
    username = os.environ['USER']
    filename = 'gti_%d.h5' % shot

    return os.path.join('/tmp', username, filename)


if __name__ == '__main__':
    shot = 42661
    rho, time, data = inverted_data(shot)
    
    import matplotlib.pyplot as plt
    plt.contour(time, rho, data.T)


