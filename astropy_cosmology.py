#!/usr/bin/env python

from __future__ import division, print_function
import numpy as np
import scipy

import astropy.cosmology
import astropy.constants as const
import astropy.units as u

cosmo = astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)


def run_comoving_distance(n=10000):
    z = np.arange(n) + 0.001
    return cosmo.comoving_distance(z)


def run_comoving_distance_elliptic(n=10000):
    z = np.arange(n) + 0.001
    return comoving_distance_elliptic(z, cosmo.H0, cosmo.Om0)


def comoving_distance_elliptic(z, H0, Om0):
    """Calculate comoving distance using incomplete elliptic integration
    """
    print("z, H0, Om0: ", z, H0, Om0)
    s = ((1-Om0)/Om0) ** (-1/3)
    prefactor = const.c/(H0*np.sqrt(s*Om0))
    prefactor = prefactor.to(u.Mpc)
    print("s, prefactor: ", s, prefactor)
    print("T(s): ", T_legendre(s))
    print("T(s/(1+z)): ", T_legendre(s/(1+z)))
    return prefactor * (T_legendre(s) - T_legendre(s/(1+z)))


def T_legendre(x):
    """Compute T(x) using Legendre elliptical integrals of the first kind.

    T(x) = 3^{-\frac{1}{4}} F\left(arccos\left(\frac{1+(1-\sqrt{3}x}{1+(1+\sqrt{3})x}\right), \cos\frac{\pi}{12}\right)
    """
    F = scipy.special.ellipkinc
    phi = np.arccos((1+(1-np.sqrt(3))*x)/(1+(1+np.sqrt(3))*x))
    k = np.cos(np.pi/12)
    print("phi, k: ", phi, k)
    return 3**(1./4) * F(phi, k)


def T_carlson(z):
    """Compute T(x) using Carlson elliptical integrals of the first kind.

    R_F(x, y, z)  Carlson 1977

    Not yet implemented in SciPy
    """
    pass


if __name__ == "__main__":
    dc = run_comoving_distance(n=100)
    dc_elliptic = run_comoving_distance_elliptic(n=100)

    print("DC - DC_elliptic: ", dc-dc_elliptic)
