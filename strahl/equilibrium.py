import numpy as np
import matplotlib.pyplot as plt

from geometry.tcv import vessel_patch


class MockUpEquilibrium(object):
    def __init__(self, minor_radius=0.23, magnetic_axis=(0.88, 0.23),
            elongation=1.4, triangularity=0.2, grid_shape=(150,150)):
        self.minor_radius = minor_radius
        self.magnetic_axis = magnetic_axis
        self.elongation = elongation
        self.triangularity = triangularity
        self.grid_shape = grid_shape

    def _rz_contour(self, r, theta):
        Delta = self.triangularity * r / self.minor_radius
        Kappa = self.elongation

        R0, Z0 = self.magnetic_axis
        R = R0 + r * np.cos(theta + Delta * np.sin(theta))
        Z = Z0 + r * Kappa * np.sin(theta)

        return R, Z

    def get_rho(self):
        r, theta = self._polar_grid()
        r *= self.minor_radius
        surface_r, surface_z = self._rz_contour(r, theta)

        return surface_r, surface_z, r

    def _polar_grid(self):
        xx = np.linspace(-1.0, 1.0, self.grid_shape[0])
        yy = np.linspace(-1.0, 1.0, self.grid_shape[1])
        x, y = np.meshgrid(xx, yy)

        theta = np.arctan2(y, x) + np.pi
        rho = np.hypot(x, y)

        return rho, theta

    def get_volume(self, n_points):
        r = np.linspace(0, self.minor_radius, n_points)
        theta = np.linspace(0, 2 * np.pi, 30)

        R0, Z0 = self.magnetic_axis
        V = np.zeros_like(r)
        for i in xrange(len(r)):
            R, Z = self._rz_contour(r[i], theta)
            V[i] = 0.5 * np.trapz((R-R0)**2 + (Z-Z0)**2, theta)

        V *= 2 * np.pi * R0

        return r, V


def plot_equilibrium(eq):
    """
    >>> eq = MockUpEquilibrium()
    >>> plot_equilibrium(eq)
    """
    ax = plt.gca()
    ax.add_patch(vessel_patch())

    r, z, rho = eq.get_rho()

    levels = np.linspace(0,eq.minor_radius,9)
    ax.contour(r, z, rho, colors='black', levels=levels)
    ax.contour(r, z, rho, colors='black', levels=[eq.minor_radius], linewidths=2)

    major_radius, z_position = eq.magnetic_axis
    ax.plot([major_radius], [z_position],'kx', ms=5)
    ax.set_aspect('equal')
    ax.figure.canvas.draw()


if __name__ == '__main__':
    import doctest
    doctest.testmod()
