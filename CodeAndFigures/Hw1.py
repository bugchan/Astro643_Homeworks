# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 20:13:17 2020

@author: sandy

1. summary
2. extended summary (optional)
3. see also (optional)
4. references (optional)
5. examples (optional)

See: https://numpydoc.readthedocs.io/en/latest/format.html

"""

import numpy as np
import matplotlib.pyplot as plt
import SetupPlots as SP

# % Definitions


def calcVinf(G, M, R, L, calcGamma, constDict):
    """Calculate v_inf."""
    Lsun = constDict['Lsun']
    c = constDict['c']
    G_CGS = constDict['G']
    Msun = constDict['Msun']
    if calcGamma:
        kappa = 0.3  # cm^2 g^-1
        Gamma = (kappa*L/(4*np.pi*c*G_CGS*M))*(Lsun/Msun)
        print(Gamma)
        v_inf = np.sqrt(2*(1 - Gamma)*G*M/R)/2.6  # km/s
    else:
        v_inf = np.sqrt(2*G*M/R)/2.6  # km/s
    return v_inf


def calcMdot(G, M, R, L, calcGamma, constDict):
    r"""Return the mass loss rate.

    Several sentences providing an extended description. Refer to
    variables using back-ticks, e.g. `var`.

    Parameters
    ----------
    Gamma : array_like
        Array_like means all those objects -- lists, nested lists, etc. --
        that can be converted to an array.  We can also refer to
        variables like `var1`.
    G : float
        The gravitational constant in units of km^2 R_sun / s^2 M_sun.
    M : float / array
        Mass in solar units.
    R : float
        Radius in solar units.
    L : Array
        Luminosity in solar units.

    Returns
    -------
    type
        Explanation of anonymous return value of type ``type``.
    describe : type
        Explanation of return value named `describe`.
    out : type
        Explanation of `out`.
    type_without_description

    Other Parameters
    ----------------
    only_seldom_used_keywords : type
        Explanation
    common_parameters_listed_above : type
        Explanation

    Raises
    ------
    BadException
        Because you shouldn't have done that.

    Notes
    -----
    Notes about the implementation algorithm (if needed).

    This can have multiple paragraphs.

    You may include some math:

    .. math:: X(e^{j\omega } ) = x(n)e^{ - j\omega n}

    And even use a Greek symbol like :math:`\omega` inline.

    References
    ----------
    Cite the relevant literature, e.g. [1]_.  You may also cite these
    references in the notes section above.

    .. [1] O. McNoleg, "The integration of GIS, remote sensing,
       expert systems and adaptive co-kriging for environmental habitat
       modelling of the Highland Haggis using object-oriented, fuzzy-logic
       and neural-network techniques," Computers & Geosciences, vol. 22,
       pp. 585-588, 1996.
    """
    Lsun = constDict['Lsun']
    c = constDict['c']
    G_CGS = constDict['G']
    Msun = constDict['Msun']
    Mdot = []
    for l in L:
        v_inf = calcVinf(G, M, R, l, calcGamma, constDict)
        temp1 = (-1.37 + 2.07*np.log10(l/(10**6)) - np.log10(v_inf*R**(1/2)))
        Mdot.append(temp1)
    return Mdot

# %% Constants


constCGS = {'G': 6.67408e-8,  # cm^3 g^-1 s^-2
            'Msun': 1.9891e33,  # g
            'Lsun': 3.847e33,  # erg/s or g cm**2 s**-2/s
            'Rsun': 6.96e10,  # cm
            'c': 2.997e10,  # cm s^-1
            }
constSI = {'G': 6.67408e-11,  # m^3 kg^-1 s^-2
           'Msun': 1.9891e30,  # kg
           'Rsun': 6.96e8,  # m
           'c': 2.997e8,  # m s^-1
           }
constants = (constCGS, constSI)

G = constSI['G']*constSI['Msun']/(constSI['Rsun']*10**6)
M = np.linspace(20, 100, 100)  # Msun
L = np.array((1e5, 3e5, 1e6, 2e6))  # Lsun
R = 1

# %% Operations

Mdot1 = calcMdot(G, M, R, L, calcGamma=False, constDict=constCGS)
Mdot2 = calcMdot(G, M, R, L, calcGamma=True, constDict=constCGS)

LM3 = M**3
v_inf3 = np.sqrt(2*G*M/R)/2.6  # km/s
Mdot3 = (-1.37 + 2.07*np.log10(LM3/(10**6)) - np.log10(v_inf3*R**(1/2)))


# %% Plot
fig = plt.figure()
gs = SP.setupPlot(fig, figsize=(1, 2.5), grid=(1, 1))
ax = fig.add_subplot(gs[0])

for i in range(len(L)):
    ax.plot(M, Mdot1[i], '-', color='C%i' % (i + 1), label=r'$L/L_\odot={:.0e}$'.format(L[i]))
    ax.plot(M, Mdot2[i], '--', color='C%i' % (i + 1))

ax.plot(M, Mdot3, color='C0', label=r'$L\propto M^3$')

# ax.set_ylim((10**-7, 10**-3.3))
ax.set_ylabel(r'$\log\dot{M}\quad[M_\odot yr^{-1}]$')
ax.set_xlim((20, 100))
ax.set_xlabel(r'$M\quad[M_\odot]$')
ax.legend(loc='lower right')
ax.grid()

fig.tight_layout()
fig.savefig('hw1problem3fig1.pdf')
