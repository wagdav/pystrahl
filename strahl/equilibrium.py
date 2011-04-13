import numpy as np
import matplotlib.pyplot as plt

import crpppy.mds as mds

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

    def plot(self, ax=None):
        ax = ax or plt.gca()

        r, z, psi = self.get_psi_contours()
        levels = np.linspace(0, 1, 9)

        ax.contour(r, z, psi, colors='black', levels=levels)
        ax.contour(r, z, psi, colors='black', levels=[1], linewidths=2)
        ax.plot([self.major_radius], [self.z_position],'kx', ms=5)
        ax.set_aspect('equal')
        ax.figure.canvas.draw()


class MockUpEquilibrium(Equilibrium):
    """
    Example:
    --------

    Circular plasma with TCV-like size:
    >>> eq = MockUpEquilibrium(magnetic_axis=(0.88, 0.23))
    >>> eq.plot()

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

    def get_volume(self, rho=None):
        """
        Give the plasma volume on the grid *radius*. If radius is None the
        total volume is given, i. e. radius = minor_radius is used.
        """
        if rho == None: rho = [1]

        radius = np.asarray(rho)
        theta = np.linspace(0, 2 * np.pi, 30)

        R0, Z0 = self.magnetic_axis
        area = []
        for r in rho:
            R, Z = self._rz_contour(r * self.minor_radius, theta)
            area.append(0.5 * np.trapz((R-R0)**2 + (Z-Z0)**2, theta))

        area = np.array(area)
        volume = 2 * np.pi * R0 * area

        return volume

    def get_psi_contours(self):
        r, theta = self._polar_grid()
        surface_r, surface_z = self._rz_contour(r * self.minor_radius, theta)

        return surface_r, surface_z, r**2

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


class LiuqeEquilibrium(Equilibrium):
    """
    >>> eq = LiuqeEquilibrium(42661,1)
    >>> eq.plot()
    """
    def __init__(self, shot, time):
        conn = mds.TCV_Connection(shot)
        rho = conn.get(r'dim_of(\results::psitbx:vol)')
        vol = conn.get(r'\results::psitbx:vol[*,$1]', time)

        r = conn.get(r'dim_of(\results::psi,0)')
        z = conn.get(r'dim_of(\results::psi,1)')
        psi = conn.get(r'\results::psi[*,*,$1]/\results::psi_axis[$1]',
                time)
        psi = 1 - psi

        grho = conn.get(r'\results::psitbx:grho[*,$1]',time)

        self.rho = np.array(rho)
        self.vol = np.array(vol)
        self.r = np.array(r)
        self.z = np.array(z)
        self.psi = np.array(psi)
        self.grho = np.array(grho)

        r0 = conn.get(r'\results::r_axis[$1]', time)
        z0 = conn.get(r'\results::z_axis[$1]', time)
        self.magnetic_axis = float(r0), float(z0)
        conn.close()

    def get_volume(self, rho_out):
        vol_out = np.interp(rho_out, self.rho, self.vol)
        return vol_out

    def get_grho(self, rho_out):
        grho_out = np.interp(rho_out, self.rho, self.grho)
        return grho_out

    def get_psi_contours(self):
        rr, zz = np.meshgrid(self.r, self.z)
        return rr, zz, self.psi


if __name__ == '__main__':
    import doctest
    doctest.testmod()
