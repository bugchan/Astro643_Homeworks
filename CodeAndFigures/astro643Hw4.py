# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 22:31:18 2020

@author: sandy

1. summary
2. extended summary (optional)
3. see also (optional)
4. references (optional)
5. examples (optional)

See: https://numpydoc.readthedocs.io/en/latest/format.html

"""

import numpy as np
import constants as cons
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import SetupPlots as SP
from matplotlib.backends.backend_pgf import FigureCanvasPgf
matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)

width, height = SP.setupPlot(False)

# matplotlib.use('Qt5Agg', warn=False)
# matplotlib.rcParams.update({
# #     "pgf.texsystem": "pdflatex",
# 'backend': 'module://ipykernel.pylab.backend_inline',
# #     'font.family': 'serif',
#     'text.usetex': False,
# #     'pgf.rcfonts': False,
# })

# %% Problem 1

M = 1.4*cons.mSun
R_WD = cons.rEarth
R_NS = 10e3  # meter
mCarbon = 12/cons.NA  # g
Qcarbon = 4.62e6*cons.eV2erg
Qsophia = 2.49e9*cons.eV2erg
Qbook = 8.25e-5  # ergs
n = 3
coeff = 3/(5-n)

EgravWD = coeff*cons.G*M**2/R_WD
EgravNS = coeff*cons.G*M**2/R_NS

Enuclear = Qbook*M/mCarbon
# Enuclear = Qsophia*M/58.6

print('Problem 1')
print('  EgravWD = {:.2e} ergs'.format(EgravWD))
print('  EgravNS = {:.2e} ergs'.format(EgravNS))
print('  Enuclear = {:.2e} ergs'.format(Enuclear))
print('')

# %% Problem 3

MvSun = 5.072
mvStar = 14

logd = (mvStar-MvSun+5)/5
d = 10**(logd)

print('Problem 3')
print('  log(d) = {:.3f}'.format(logd))
print('  d = {:.3f}'.format(d))

P = 25*cons.day2sec  # seconds
vobs = 20e5  # cm/s

sini = vobs*(2*P/(np.pi*cons.G*cons.mSun))**(1/3)
inclination = np.arcsin(sini)
print('  sin(i) = {:.3f}'.format(sini))
print('  i = {:.3f} rad = {:.3f} deg'.format(inclination,
                                             np.rad2deg(inclination)))

a = (P**2*cons.G*2*cons.mSun/(4*np.pi**2))**(1/3)  # cm
print('  a = {:.3e} AU'.format(a/cons.AU2cm))
eclipseAngle = np.rad2deg(np.arccos((2*cons.rSun)/a))
print('  Eclipsing Angle = {:.3f} deg'.format(eclipseAngle))
hubbleRes = 551e-9/2.3
hubbleResSec = np.rad2deg(hubbleRes)*60*60
print('  Hubble\'s resolution = {:.3e}rad'.format(hubbleRes))
print('  Hubble\'s resolution = {:.3e}arcsec'.format(hubbleResSec))

# %% Problem 4

MainSequence = pd.read_csv('ASTRO643_HW4P4_Tables - MainSequence.csv')
Giants = pd.read_csv('ASTRO643_HW4P4_Tables - Giants.csv')
SuperGiants = pd.read_csv('ASTRO643_HW4P4_Tables - SuperGiants.csv')
extinction = pd.read_csv('ASTRO643_HW4P4_Tables - Extinction.csv')


# %%
MS_JH = MainSequence['(V-H)0']-MainSequence['(V-J)0']
MS_HK = MainSequence['(V-K)0']-MainSequence['(V-H)0']

G_JH = Giants['(V-H)0']-Giants['(V-J)0']
G_HK = Giants['(V-K)0']-Giants['(V-H)0']

SG_JH = SuperGiants['(V-H)0']-SuperGiants['(V-J)0']
SG_HK = SuperGiants['(V-K)0']-SuperGiants['(V-H)0']


# # %% Reddening Vector

A_jv = extinction[extinction.Filter == 'J']['A_lambda/A_V'].iloc[0]
A_hv = extinction[extinction.Filter == 'H']['A_lambda/A_V'].iloc[0]
A_kv = extinction[extinction.Filter == 'K']['A_lambda/A_V'].iloc[0]
# color excess in terms of Av
E_JHv = A_jv - A_hv
E_HKv = A_hv - A_kv
Rv = 3.1
E_BV = 1
# calculation of Av
Av = Rv*E_BV
# real color excess
E_JH = E_JHv*Av
E_HK = E_HKv*Av
# origin for arrow
x0 = 0.2
y0 = 0.3

# # %% NIR-Object
J_mag = 12.66
H_mag = 11.77
K_mag = 11.51

NIR_JH = J_mag - H_mag
NIR_HK = H_mag - K_mag

# # %% Plot
fig1 = plt.figure(figsize=(1.5*width, 1.5*height))
ax1 = fig1.add_subplot()
ax1.plot(MS_HK, MS_JH, '-', color='C0', label='Main Sequence')
ax1.plot(G_HK, G_JH, '-', color='C1', label='Giants')
ax1.plot(SG_HK, SG_JH, '-', color='C2', label='SuperGiants')

ax1.arrow(x0, y0, E_HK, E_JH, head_width=.03, fc='black',
          linewidth=1)
ax1.text(x0+.1, y0, r'$A_v=3.1$')

ax1.errorbar(NIR_HK, NIR_JH, xerr=0.06, yerr=0.09, fmt='*', color='C4',
             markersize=5, label='NIR Object')

ax1.set_xlabel(r'$(H-K)_0$')
ax1.set_ylabel(r'$(J-H)_0$')
xlim1, xlim2 = -1.1, 0.6
ylim1, ylim2 = -0.25, 1.4
ax1.set_xlim(xlim1, xlim2)
ax1.set_ylim(ylim1, ylim2)

# extinction for Near Infrared Object
Av_NIR2MS = -.4
Av_NIR2G = -1.5
Av_NIR2SG = -1.9


ax1.errorbar(NIR_HK+E_HKv*Av_NIR2MS, NIR_JH+E_JHv*Av_NIR2MS,
             xerr=0.06, yerr=0.09, fmt='*', color='C0',
             markersize=3, label='Av = {:.2f}'.format(Av_NIR2MS))
ax1.errorbar(NIR_HK+E_HKv*Av_NIR2G, NIR_JH+E_JHv*Av_NIR2G,
             xerr=0.06, yerr=0.09, fmt='*', color='C1',
             markersize=3, label='Av = {:.2f}'.format(Av_NIR2G))
ax1.errorbar(NIR_HK+E_HKv*Av_NIR2SG, NIR_JH+E_JHv*Av_NIR2SG,
             xerr=0.06, yerr=0.09, fmt='*', color='C2',
             markersize=3, label='Av = {:.2f}'.format(Av_NIR2SG))

ax1.grid()
ax1.legend()

fig1.savefig('ASTRO643_HW4P4plot.pdf')


# %% Problem 5

def IMFunc(m, index, coeff=1):
    return coeff*m**(-index)


# constants
expa = 1.3      # index for .1 < M < .5mSun
expb = 2.35     # index for >= .5mSun
minMass = .1    # mSun
medMass = .5    # mSun
maxMass = 100   # mSun
n = 6           # mSun yr^-1

# Math
M = np.logspace(np.log10(minMass), np.log10(maxMass), 100)

# C in terms of D
AB = IMFunc(medMass, expb, 1) / IMFunc(medMass, expa, 1)  # B

# integration
a = (-expa + 1.)
b = (-expb + 1.)
AxM = (1/a)*(medMass**a - minMass**a)  # xA
BxM = (1/b)*(maxMass**b - medMass**b)  # xB
nB = AxM*AB + BxM
B = n/nB
A = AB*B

print('Problem 5')
print('  IMF: {:.2f}*M^{:.2f} for .1 < M < .5mSun'.format(A, expa))
print('  IMF: {:.2f}*M^{:.2f} for M >= .5mSun'.format(B, expb))

# IMF
IMF_a = IMFunc(M[M < medMass], expa, A)
IMF_b = IMFunc(M[M >= medMass], expb, B)
IMF = np.concatenate((IMF_a, IMF_b))

# Supernovae rate
SNRyear_100 = (B/b)*(100**b-8**b)
SNRyear_inf = (B/b)*(-8**b)
print('  SNR from 8-100 Msun per century: {:.2f}'.format(SNRyear_100*100))
print('  SNR from 8-inf Msun per century: {:.2f}'.format(SNRyear_inf*100))


fig2 = plt.figure(figsize=(1.5*width, 1.5*height))
ax4 = fig2.add_subplot(1, 1, 1)
ax4.loglog(M, IMF)
ax4.text(.3, 10, r'$M^{-1.3}$', fontsize=6)
ax4.text(10, .01, r'$M^{-2.35}$', fontsize=6)
ax4.set_xlabel(r'm $[M_\odot]$')
ax4.set_ylabel(r'IMF $ [M_\odot] yr^{-1} $')
ax4.grid(True, which='both')

fig2.savefig('ASTRO643_HW4P5plot.pdf')
