# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 15:16:05 2020

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

# % Definitions


# % Constants


NA = 6.022e23
e = 4.8032068e-10
c = 2.99792458e10
B = 9.739e5
me = 9.1093898e-28
a = 7.56577e-15
R_sun = 6.96e10
k_B = 1.380650e-16
h = 6.6260688e-27


# Numeric values:
L = (1/(NA*4))*((e**2)/(me*(c**2)*((1./B)**(2./3.))))**(-2.)
K = 4*a*c/3

# print(L)

# %% White Dwarf

rho = 10**6
mu_I = 12
mu_e = 2
Z_c = 6
R = 0.01*R_sun

T = 10**7

# L_whitedwarf = L*(rho/mu_e)**(4/3)*(mu_I/rho)*(1/Z_c**(2))
L_cons = (me**2*c**4) / (NA*np.pi*e**4*B**(4/3))
L2_2 = L_cons * (rho**(1/3)*mu_I)/(Z_c**2*mu_e**(4/3))
L2 = (mu_I*me**2*c**4*rho**(4/3)) / \
     (rho*NA*4*np.pi*Z_c**2*e**4*mu_e**(4/3)*B**(4/3))
v_e = c*(rho/(mu_e*B))**(1/3)
# t = R**2/(L_whitedwarf*v_e)
t2_cons = B**(1/3)/(L*c)
t2 = t2_cons * (R**2*Z_c**2*rho/mu_I)*(mu_e/rho)**(5/3)
# cv = (8*(np.pi**3)*(me**2)*(k_B**2)*T*v_e)/(3*h**3)
cv2 = 8*np.pi**3*me**2*c*k_B**2*T*rho**(1/3)/(3*h**3*mu_e**(1/3)*B**(1/3))
# De = cv*v_e*L_whitedwarf/3
# De2 = (8*np.pi**3*me**2*c*k**2*T*rho**(1/3)/(3*h**3*mu_e**(1/3)*B**(1/3)))*\
#       (c*rho**(1/3)/(mu_e**(1/3)*B**(1/3)))*L2/3
De2 = (2*np.pi**2*me**4*c**6*k_B**2*T*mu_I*rho) / \
      (9*h**3*NA*Z_c**2*e**4*mu_e**2*B**2)
# k_cond = (K*T**3) / (De*rho)
k_cons = (3*h**3*NA*e**4*B**2*a) / (2*np.pi**2*me**4*c**5*k_B**2)
k_cond2 = k_cons*mu_e**2*Z_c**2*T**2 / (rho**2*mu_I)

# k_cond = 4*a*c*T**3/(3*De*rho)
# print('White Dwarf Lambda: %.3e' % (L_whitedwarf))
print('White Dwarf Lambda: %.3e' % (L2*np.pi))
print('White Dwarf v_e: %.3e' % (v_e))
# print('White Dwarf t: %.3e' % (t))
print('White Dwarf t2_cons: %.3e' % (t2_cons))
print('White Dwarf t2: %.3e' % (t2))
# print('White Dwarf c_V: {:.3e}'.format(cv))
print('White Dwarf c_V: {:.3e}'.format(cv2))
# print('White Dwarf D_e: {:.3e}'.format(De))
print('White Dwarf D_e: {:.3e}'.format(De2*np.pi))
# print('White Dwarf k: %.3e' % (k_cond))
print('White Dwarf k: {:.3e}'.format(k_cond2/np.pi))
print('White Dwarf k: {:.3e}'.format(k_cond2/4))
print('')

# %% Sun

# rho = 1.41
# rho = 0.2  # NASA density at tachokline
rho = 1
mu_I = 1.3
mu_e = 1.2
# Z_sun = 1.3
Z_sun = 1
Gamma = 5/3
# A = 1.9
A = 1
T = 10**7
# kappa_t = 1
kappa_t = 0.2*(1+.7)  # from Daniel notes on opacity
mu = 1/((1/mu_I)+(1/mu_e))

L_sun = L*(rho/mu_e)**(4/3)*(mu_I/rho)*(1/Z_sun**(2))
v_e = c*(rho/(mu_e*B))**(1/3)
cv = (8*(np.pi**3)*(me**2)*(k_B**2)*T*v_e)/(3*h**3)
De = cv*v_e*L_sun/3
k_cond = (K*T**3) / (De*rho)
nu_T = v_e*L_sun/(3*rho*Gamma)
nu_T2 = (8*a*c*T**3*mu) / (15*kappa_t*rho**2*NA*k_B)
ne = rho*NA/mu_e
nu = (2e-15*T**(5/2)*A**(1/2))/(rho*Z_sun**4*np.log(10**4*T**(3/2)*ne**(-1/2)))
PR2 = nu / nu_T2

print('Sun Lambda: %.3e' % (L_sun))
print('Sun v_e: %.3e' % (v_e))
# print('Sun nu_T: %e' % (nu_T))
print('Sun nu_T2: %.3e' % (nu_T2))
print('Sun nu: %.3e' % (nu))
print('Sun k_Cond: %.3e' % (k_cond))
# print('Sun PR: %e' % (nu / nu_T))
print('Sun PR2: %.3e' % (nu / nu_T2))

# %% problem 4

l = np.linspace(0, 10, 10)
Nsq = np.logspace(-3, 3, 10)

Ra = -Nsq*l**4


# %% problem 5
eMg = -8.26
nMg = 24
Z_Mg = 12
A_Mg = 24
eHe = -7.074
nHe = 4
A_He = 4
Z_He = 2
eSi = -8.447
nSi = 28
mu = (A_He**(-1)+A_Mg**(-1))**(-1)

Q = nMg*eMg + nHe*eHe - nSi*eSi
print('Q-value: %e' % (Q))

T6 = ((Q*10**3)/(1.22*(Z_He**2*Z_Mg**2*mu)**(1/3)))**(3/2)
print('Temp:    %e' % (T6*10**6))
