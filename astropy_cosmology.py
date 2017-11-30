#!/usr/bin/env python

from __future__ import division, print_function

import math

import numpy as np
import scipy.special

import astropy.cosmology
import astropy.constants as const
import astropy.units as u

cosmo = astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)


def run_comoving_distance(n=10000, z_max=10.0):
    z = (z_max/n) * np.arange(n) + 0.001
    return z, cosmo.comoving_distance(z)


def run_comoving_distance_elliptic(n=10000, z_max=10.0):
    z = (z_max/n) * np.arange(n) + 0.001
    return z, comoving_distance_elliptic(z, cosmo.H0, cosmo.Om0)


def comoving_distance_elliptic(z, H0, Om0):
    """Calculate comoving distance using incomplete elliptic integration
    """
    s = ((1-Om0)/Om0) ** (1/3)
    prefactor = (const.c/H0)*(1/np.sqrt(s*Om0))
    prefactor = prefactor.to(u.Mpc)  # normalize from km/s/Mpc
    return prefactor * (T_legendre(s) - T_legendre(s/(1+z)))


def T_legendre(x):
    """Compute T(x) using Legendre elliptical integrals of the first kind.

    T(x) = 3^{-\frac{1}{4}}
           F\left(arccos\left(\frac{1+(1-\sqrt{3}x}{1+(1+\sqrt{3})x}\right),
                  \cos\frac{\pi}{12}\right)
    where
    F(phi, m) = int_0^phi {1/sqrt(1 - m \sin^2{t})} dt
    """
    F = scipy.special.ellipkinc
    # math.sqrt is several times faster than np.sqrt for scalars
    phi = np.arccos((1 + (1-math.sqrt(3))*x) /
                    (1 + (1+math.sqrt(3))*x))
    k = np.cos(np.pi/12)
    m = k*k  # np.cos(np.pi/12)**2 == 1/2 + math.sqrt(3)/4
    # ellipkinc expects m=k*k as its second argument.
    return 3**(-1./4) * F(phi, m)


def T_carlson(z):
    """Compute T(x) using Carlson elliptical integrals of the first kind.

    R_F(x, y, z)  Carlson 1977

    Not yet implemented in SciPy
    """
    pass


def mpc_to_modulus(dist):
    modulus = 5*np.log10(dist/u.Mpc) + 25
    return modulus


def plot_z_dc_elliptic(z, dc, elliptic, plotname='z_dc_elliptic.pdf',
                       debug=False):
    import matplotlib.pyplot as plt
    if debug:
        print("DC - DC_elliptic: ", dc-dc_elliptic)
    fig = plt.figure(figsize=(6,8))
    plt.subplot(2,1,1)
    plt.plot(z, dc, label='dc')
    plt.plot(z, dc_elliptic, label='dc_elliptic')
    plt.plot(z, dc - dc_elliptic, label='dc - dc_elliptic')
    plt.plot(z, 10000*(dc/dc_elliptic), linestyle='--',
             label='10,000 * (dc / dc_elliptic)')
    plt.ylabel('Mpc')
    plt.xlabel('z')
    plt.ylim(10000*np.array([-1, +1]))
    plt.legend()

    dc_mag = mpc_to_modulus(dc)
    dc_elliptic_mag = mpc_to_modulus(dc_elliptic)
    plt.subplot(2,1,2)
    plt.plot(z, dc_mag, label='dc')
    plt.plot(z, dc_elliptic_mag, label='dc_elliptic')
    plt.plot(z, dc_mag - dc_elliptic_mag, label='dc - dc_elliptic')
    plt.ylabel('Distance Modulus')
    plt.xlabel('z')

    plt.savefig(plotname)


if __name__ == "__main__":
    z, dc = run_comoving_distance(n=100)
    z, dc_elliptic = run_comoving_distance_elliptic(n=100)

    plot = False
    if plot:
        plot_z_dc_elliptic(z, dc, dc_elliptic)
