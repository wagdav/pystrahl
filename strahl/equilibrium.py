import numpy as np
import matplotlib.pyplot as plt

from geometry.tcv import vessel_patch

class Equilibrium(object):
    """
    Abstract class to define an uniform equilibrium interface.
    """
    def get_volume(self, radius=None):
        pass

    @property
    def major_radius(self):
        return self.magnetic_axis[0]

    @property
    def z_position(self):
        return self.magnetic_axis[1]

class MockUpEquilibrium(Equilibrium):
    """
    Example:
    --------

    Circular plasma with TCV-like size:
    >>> eq = MockUpEquilibrium(magnetic_axis=(0.88, 0.23))

    The volume of a this can be easily calculated analytically as:
    >>> from math import pi
    >>> volume_an = 2 * pi**2 * 0.88 * 0.23**2

    We compare the two values:
    >>> volume_eq, = eq.get_volume()
    >>> round(volume_an, 4) == round(volume_eq, 4)
    True
    """
    def __init__(self, minor_radius=0.23, magnetic_axis=(0.88, 0.23),
            elongation=1.0, triangularity=0.0, grid_shape=(150,150)):
        self.minor_radius = minor_radius
        self.magnetic_axis = magnetic_axis
        self.elongation = elongation
        self.triangularity = triangularity
        self.grid_shape = grid_shape

    def get_volume(self, radius=None):
        """
        Give the plasma volume on the grid *radius*. If radius is None the
        total volume is given, i. e. radius = minor_radius is used.
        """
        if radius == None: radius = [self.minor_radius]

        radius = np.asarray(radius)
        theta = np.linspace(0, 2 * np.pi, 30)

        R0, Z0 = self.magnetic_axis
        area = []
        for r in radius:
            R, Z = self._rz_contour(r, theta)
            area.append(0.5 * np.trapz((R-R0)**2 + (Z-Z0)**2, theta))

        area = np.array(area)
        volume = 2 * np.pi * R0 * area

        return volume

    def get_radius_contours(self):
        r, theta = self._polar_grid()
        r *= self.minor_radius
        surface_r, surface_z = self._rz_contour(r, theta)

        return surface_r, surface_z, r

    def get_rhovol(self, r):
        """
        Return rho_vol = sqrt (volume / (2 * pi**2 * R_0)).

        >>> circular = MockUpEquilibrium(magnetic_axis=(0.88, 0.23))
        >>> rho_vol, = circular.get_rhovol([0.23])
        >>> print round(rho_vol,2)
        0.23
        """
        volume = self.get_volume(r)
        return np.sqrt(volume/(2 * np.pi**2 * self.major_radius))

    def _polar_grid(self):
        xx = np.linspace(-1.0, 1.0, self.grid_shape[0])
        yy = np.linspace(-1.0, 1.0, self.grid_shape[1])
        x, y = np.meshgrid(xx, yy)

        theta = np.arctan2(y, x) + np.pi
        rho = np.hypot(x, y)

        return rho, theta

    def _rz_contour(self, r, theta):
        Delta = self.triangularity * r / self.minor_radius
        Kappa = self.elongation

        R0, Z0 = self.magnetic_axis
        R = R0 + r * np.cos(theta + Delta * np.sin(theta))
        Z = Z0 + r * Kappa * np.sin(theta)

        return R, Z

    def plot(self, ax=None):
        """
        >>> eq = MockUpEquilibrium()
        >>> eq.plot()
        """
        ax = ax or plt.gca()
        ax.add_patch(vessel_patch())

        r, z, rho = self.get_radius_contours()

        levels = np.linspace(0, self.minor_radius,9)
        ax.contour(r, z, rho, colors='black', levels=levels)
        ax.contour(r, z, rho, colors='black', levels=[self.minor_radius], linewidths=2)

        ax.plot([self.major_radius], [self.z_position],'kx', ms=5)
        ax.set_aspect('equal')
        ax.figure.canvas.draw()


if __name__ == '__main__':
    import doctest
    doctest.testmod()
