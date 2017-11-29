#!/usr/bin/env python

import numpy as np
import scipy

import astropy.cosmology

cosmo = astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)


def run_comoving_distance(n=10000):
    z = np.arange(n) + 0.001
    cosmo.comoving_distance(z)


def run_elliptic_comoving_distance(n=10000):
    z = np.arange(n) + 0.001


def elliptic_integral_distance(z):
    scipy.special.ellipseinc()



if __name__ == "__main__":
    run_comoving_distance(n=100000)
