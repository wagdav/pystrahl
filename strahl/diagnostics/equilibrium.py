import numpy as np
import matplotlib.pyplot as plt

major_radius = 0.88
minor_radius = 0.25
z_position = 0.23

elongation = 1.6
triangularity = 0.3

magnetic_axis = (major_radius, z_position)


grid_size = (50, 50)
xx = np.linspace(-minor_radius, +minor_radius, grid_size[0])
yy = np.linspace(-elongation * minor_radius, elongation * minor_radius,
        grid_size[1])
x, y = np.meshgrid(xx, yy)

theta = np.arctan2(y, x)
rhoc = np.hypot(x, y)

R = rhoc * np.cos(theta - triangularity * np.sin(theta))
Z = elongation * rhoc * np.sin(theta)

rho = np.hypot(R, Z) / minor_radius


def plot_flux_function(r, z, f):
    ax = plt.gca()

    levels = np.linspace(0,1,9)
    cb = ax.contourf(r, z, f, levels)
    ax.contour(r, z, f, colors='black', levels=levels)
    plt.colorbar(cb)

    ax.plot([magnetic_axis[0]], [magnetic_axis[1]],'x', ms=5)
    ax.set_aspect('equal')


if __name__ == '__main__':
    plt.figure(1); plt.clf()
    ax = plt.gcf().add_subplot(111)

    plot_flux_function(x + magnetic_axis[0], y + magnetic_axis[1], rho)
    plt.draw()

    plt.show()



