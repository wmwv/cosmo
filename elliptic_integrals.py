#!/usr/bin/env python

from itertools import cycle

import matplotlib.pyplot as plt
import numpy as np
import scipy.special

"""
Explore Elliptic Integrals
"""


def elliptic_integral_first_integrand(theta, m):
    """Integrand of elliptic integral of the first kind

    m=k*k, where k is the elliptic modulus
    """
    return 1/np.sqrt(1-m*np.sin(theta))


def plot_elliptic_integrand(n=1000):
    """Plot the integrand of the elliptic function."""
    theta_max = np.pi
    theta = (theta_max/(n-1)) * np.arange(n)
    # m = k*k
    m_values = [0.01, 0.1, 0.2, 0.5, 0.8, 0.9, 0.99]
    # Let's go through these from high to low
    #   to match the height they will appear on the plot
    m_values = m_values[::-1]
    colors = cycle(['black', 'red', 'green', 'blue', 'yellow', 'orange', 'purple'])
    for m, c in zip(m_values, colors):
        plt.plot(theta, elliptic_integral_first_integrand(theta, m),
                 label="m = %.2f" % m, linestyle='--', color=c)
        plt.plot(theta, scipy.special.ellipkinc(theta, m), linestyle='-',
                 color=c)
    plt.xlabel(r'$\theta$')
    plt.legend()


if __name__ == "__main__":
    plot_elliptic_integrand()
    plt.show()
