# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 15:39:04 2020

@author: Sandra Bustamante

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import constants as cons

#%% Problem 2

def RK4(dvdt, a, b, h, IV, dim=2, *args):
    """Integrate 1st Order Diferrential Equations using RK4. We are assuming
    that v is the first derivative of x and dvdt is the second derivative of 
    x.

    Parameters
    ----------
    dvdt : definition
        Function of the differential equation to solve.
    a : float
        Initial value of time.
    b : float
        Last value of time.
    h : float
        Step size.
    IV : array
        An array of the initial values of [x0,v0] at time 0.
    dim : int, optional
        Dimensions of the solutions. If x is space then 
        2 dimension means there is an x and y directions. The default is 2.

    Returns
    -------
    t : Array
        Time array.
    x : Array
        Solution of the position array.
    v : Array
        Solution of the velocity array.

    """
    t = np.arange(a, b+h, h)          # create time
    x = np.zeros((len(t), dim))      # initialize x
    v = np.zeros((len(t), dim))      # initialize x
    x[0], v[0] = IV                 # set initial values
    # apply Fourth Order Runge-Kutta Method
    for i in np.arange(1, len(t)):
        if x[i-1,0] >= 0: # stops after the first negative number of x
            k1 = h*dvdt(t[i-1], x[i-1], v[i-1], *args)
            j1 = h*v[i-1]
            k2 = h*dvdt(t[i-1]+h/2.0, x[i-1]+j1/2.0,
                        v[i-1] + k1/2.0, *args)
            j2 = h*(v[i-1]+k1/2.0)
            k3 = h*dvdt(t[i-1]+h/2.0, x[i-1]+j2/2.0,
                        v[i-1] + k2/2.0, *args)
            j3 = h*(v[i-1]+k2/2.0)
            k4 = h*dvdt(t[i], x[i-1] + j3, v[i-1] + k3, *args)
            j4 = h*(v[i-1] + k3)
            v[i] = v[i-1] + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
            x[i] = x[i-1] + (j1 + 2.0*j2 + 2.0*j3 + j4)/6.0
        else:
            break
    if dim == 1:
        x = x.reshape(len(t))
        v = v.reshape(len(t))
    return t, x, v

def dzdx(x, y, z, n):
    """
    Uses the nomenclature as in the textbook where x=xi, y=theta and z=dydx 
    """    
    return  -y**n-2*z/x

nArray = np.array((0, 1, 1.5, 2, 3, 4))
x_a = 0
x_b = 15
h = .01 
y0 = 1
z0 = 0
IV = [y0,z0]
val = 0

fig = plt.figure()
ax = fig.add_subplot()
data = []


for i in range(len(nArray)):
    n = nArray[i]
    x, y, z = RK4(dzdx,h, x_b, h, IV, 1, n)
    index = np.argmin(np.abs(y-val))
    # interpolate to find the value of x for which y=0
    x1 = np.interp(0., y[index-1:index], x[index-1:index])
    z1 = np.interp(0., y[index-1:index], z[index-1:index])
    y=np.trim_zeros(y,'b')
    ax.plot(x[:len(y)], y, color='C{:}'.format(i), label = 'n = {:.1f}'.format(n))
    ax.plot(x1, 0, 'o', color='C{:}'.format(i))
    #save the value in data
    data.append(np.array((n, x1, -z1, (-1/3)*(x1/z1))))

# Stylize plot    
ax.grid()
ax.legend()
ax.set_ylim(-.1, 1)
ax.set_xlim(x_a, x_b)
ax.set_xlabel(r'$\xi$')
ax.set_ylabel(r'$\theta_n$')
fig.tight_layout()
fig.savefig('Astro643Hw3P2Plot.png')


# create dataframe
df = pd.DataFrame(data, columns=('n', 'x1', '-z1','ratio'))
# values of table 7.1 HKT
x1True = np.array((np.sqrt(6),np.pi, 3.6538, 4.3529, 6.8969, 14.972))
z1True = np.array((np.sqrt(6)/3, 1/np.pi, 0.20330, 0.12725, 0.04243, 0.00802))
ratioTrue = np.array((1, np.pi**2/3, 5.9907, 11.402, 54.183, 622.41))
df.insert(loc=2, column='xTrue', value=x1True)
df.insert(loc=4, column='zTrue', value=z1True)
df.insert(loc=6, column='ratioTrue', value=ratioTrue)
print(df)

# Save table in tex file
df.rename(columns = {'n':r'n',
                     'x1':r'$\xi_1$',
                     'xTrue':r'$\xi_1$ True', 
                     '-z1':r'$-\theta_n^\prime(\xi_1)$',
                     'zTrue':r'$-\theta_n^\prime(\xi_1)$ True',
                     'ratio':r'$\rho_c/\langle\rho\rangle$',
                     'ratioTrue':r'$\rho_c/\langle\rho\rangle$ True',
                     }, inplace = True)

with open('Astro643Hw3P2Table.tex', 'w') as tf:
    tf.write(df.to_latex(float_format='%.3f',
                         index=False,
                         escape=False))
    
# %% Problem 3

print('Problem 3:')
n3 = 3
i3 = np.where(nArray==n3)[0][0]
    
K = (1/(4*np.pi))*np.sqrt(3/2)*(cons.h*cons.c/(cons.G*cons.mA**(4/3)))**(3/2)
M = K*x1True[i3]**2*z1True[i3]/4
M_solarUnits = M/cons.mSun

print('  M = {:.3e}(2/mu_e)**2'.format(M))
print('  M/Msun = {:.3f}(2/mu_e)**2'.format(M_solarUnits))

# %% Problem 4
print('Problem 4')

n4 = 3/2
i4 = np.where(nArray==n4)[0][0]
#convective with solar composition
mu = 0.61
mu_e = 1.17


rhoCentral = (x1True[i4]/z1True[i4])*(1/(4*np.pi))
rhoCentral_solarUnits = rhoCentral*cons.mSun*cons.rSun**(-3)
print('  rhoCentral = {:.3f} M R^-3'.format(rhoCentral))
print('  rhoCentral = {:.3f} (M/Msun) (R/Rsun)^-3'.format(rhoCentral_solarUnits))

Sophia = 2/(5*z1True[i4]*x1True[i4])*mu
Tcentral = 2*cons.G/(5*z1True[i4]*x1True[i4]*cons.NA*cons.kB)
Tcentral_solarUnits = Tcentral*cons.mSun*cons.rSun**-1
print('  TempCentral = {:.3e} mu M R^-1'.format(Tcentral))
print('  TempCentral = {:.3e} M R^-1'.format(Tcentral*mu))
print('  TempCentral = {:.3e} mu (M/Msun) (R/Rsun)^-1'.format(Tcentral_solarUnits))
print('  TempCentral = {:.3e} (M/Msun) (R/Rsun)^-1'.format(Tcentral_solarUnits*mu))

Tbpart1 = 1.004e13*(25/cons.G**2)*(cons.NA*cons.kB)/(4**(5/3)*np.pi**(2/3))
Tbpart2 = x1True[i4]**(8/3)*z1True[i4]**(4/3)
Tcentralb = (Tbpart1*Tbpart2)**(-1)
Tcentralb_solarUnits = Tcentralb*cons.mSun**(4/3)
print('  When Pgas=Pe, Tcentral = {:.3e} mu mu_e^5/3 M^4/3'.format(Tcentralb))
print('       or {:.2e} M'.format(Tcentralb*mu*mu_e**(5/3)))
print('       or {:.3e} mu mu_e^5/3 (M/Msun)^4/3'.format(Tcentralb_solarUnits))
print('       or {:.2e} (M/Msun)^4/3'.format(Tcentralb_solarUnits*mu*mu_e**(5/3)))

Tcentralc = 4e6 #K
mu_i = (1/mu -1/mu_e)**(-1)
print('  mu_i = {:.3e}'.format(mu_i))
Mpart1 = (1.004e13*cons.NA*cons.kB*25*cons.G**(-2))**(3/4)*4**(-5/4)
Mpart2 = np.pi**(-1/2)*x1True[i4]**(2)*z1True[i4]*Tcentralc**(3/4)*mu_e**(-5/4)
M = Mpart1*Mpart2*mu**(-3/4)
# M = Mpart1*Mpart2*mu_i**(-3/4)
# M2 = (Tcentralc*(Tcentralb*mu*mu_e**(5/3))**(-1))**(3/4)
print('  M = {:.3e} g'.format(M))
print('  M = {:.3f} Msun'.format(M/cons.mSun))

# %% Problem 5

print('Problem 5')
n5 = 3/2
i5 = np.where(nArray==n5)[0][0]

PpCons = (8*np.pi*cons.G*cons.sigma/(3*2.5e-31))**(2/3)*(cons.NA*cons.kB)**(1/3)
print('  Pp = {:.3e} (Z/0.02)^-2/3 (M/L)^2/3 mu^-1/3 T_eff^-3'.format(PpCons))

KprimePart1 = (cons.NA*cons.kB*x1True[i5])**(5/2)*z1True[i5]**(1/2)
KprimePart2 = (5/(2*cons.G))**(3/2)*(4*np.pi*cons.sigma)**(3/4)*(4*np.pi)**(-1)
KprimeCons = KprimePart1*KprimePart2
print('  Kprime = {:.3e} mu^-5/2 M^-1/2 L^-3/4 T_eff^3'.format(KprimeCons))

TeffCons = ((2**(2/3)*PpCons/KprimeCons)**(2/5)*((8/5)**(-2/9)))**(5/17)
TeffCons = 2**(4/51)*PpCons**(10/85)*KprimeCons**(-10/85)
print('  TeffCons = {:.3e} (Z/0.02)^-4/51 mu^13/51 M^7/51 L^1/102'.format(TeffCons))