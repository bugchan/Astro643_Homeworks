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

matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})

#%% Problem 4

MainSequence = pd.read_csv('ASTRO643_HW4P4_Tables - MainSequence.csv')
Giants = pd.read_csv('ASTRO643_HW4P4_Tables - Giants.csv')
SuperGiants = pd.read_csv('ASTRO643_HW4P4_Tables - SuperGiants.csv')
extinction = pd.read_csv('ASTRO643_HW4P4_Tables - Extinction.csv')

MS_JH = MainSequence['(V-H)0']-MainSequence['(V-J)0']
MS_HK = MainSequence['(V-K)0']-MainSequence['(V-H)0']

G_JH = Giants['(V-H)0']-Giants['(V-J)0']
G_HK = Giants['(V-K)0']-Giants['(V-H)0']

SG_JH = SuperGiants['(V-H)0']-SuperGiants['(V-J)0']
SG_HK = SuperGiants['(V-K)0']-SuperGiants['(V-H)0']

A_jv = extinction[extinction.Filter=='J']['A_lambda/A_V'].iloc[0]
A_hv = extinction[extinction.Filter=='H']['A_lambda/A_V'].iloc[0]
A_kv = extinction[extinction.Filter=='K']['A_lambda/A_V'].iloc[0]
E_JHv = A_jv - A_hv
E_HKv = A_hv - A_kv


#%% Reddening Vector

Rv = 3.1
E_BV = 1
Av = Rv*E_BV
E_JH = E_JHv*Av
E_HK = E_HKv*Av
x0 = -.4
y0 = 0

#%% NIR-Object

J_mag = 12.66 
H_mag = 11.77
K_mag = 11.51

NIR_JH = J_mag - H_mag
NIR_HK = H_mag - K_mag

# %% Plot

fig1 = plt.figure(1, figsize = (5, 7))
ax1 = fig1.add_subplot(3, 1, 1)
ax1.scatter(MS_HK, MS_JH, color='C0', label = 'Main Sequence')
ax1.arrow(x0,y0,E_HK,E_JH, head_width=.05, fc='black',
          label = 'Reddening E(B-V)')
ax1.plot(NIR_HK, NIR_JH,'*',color='C4', label='NIR Object')
ax1.set_ylabel(r'$(J-H)_0$')
# ax1.set_xlabel(r'$(H-K)_0$')
ax1.set_xlim(-1.2,.8)
ax1.set_ylim(-.5,1.5)
ax1.grid()
ax1.legend()


ax2 = fig1.add_subplot(3, 1, 2, sharey=ax1, sharex = ax1)
ax2.scatter(G_HK, G_JH, color='C1', label = 'Giants')
ax2.arrow(x0,y0,E_HK,E_JH, head_width=.05, fc='black',
          label = 'Reddening E(B-V)')
ax2.plot(NIR_HK, NIR_JH,'*',color='C4', label='NIR Object')
ax2.set_ylabel(r'$(J-H)_0$')
# ax2.set_xlabel(r'$(H-K)_0$')
ax2.grid()
ax2.legend()


ax3 = fig1.add_subplot(3, 1, 3, sharey=ax1, sharex = ax1)
ax3.scatter(SG_HK, SG_JH, color='C2', label = 'SuperGiants')
ax3.arrow(x0,y0,E_HK,E_JH, head_width=.05, fc='black',
          label = 'Reddening E(B-V)')
ax3.plot(NIR_HK, NIR_JH,'*',color='C4', label='NIR Object')
ax3.set_ylabel(r'$(J-H)_0$')
# ax3.set_xlabel(r'$(H-K)_0$')
ax3.grid()
ax3.legend()

fig1.tight_layout()

fig1.savefig('ASTRO643_HW4P4.pgf')

# plt.annotate('Reddening E(B-V)=1',(1.8*x0,y0+.2))













