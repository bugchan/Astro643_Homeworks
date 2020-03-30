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
    """Integrate 1st Order Diferrential Equations using RK4.

    Parameters
    ----------
    dvdt : definition
        DESCRIPTION.
    a : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    h : TYPE
        DESCRIPTION.
    IV : TYPE
        DESCRIPTION.
    dim : TYPE, optional
        DESCRIPTION. The default is 2.

    Returns
    -------
    t : TYPE
        DESCRIPTION.
    x : TYPE
        DESCRIPTION.
    v : TYPE
        DESCRIPTION.

    """
    t = np.arange(a, b+h, h)          # create time
    x = np.zeros((len(t), dim))      # initialize x
    v = np.zeros((len(t), dim))      # initialize x
    x[0], v[0] = IV                 # set initial values
    # apply Fourth Order Runge-Kutta Method
    for i in np.arange(1, len(t)):
        if x[i-1,0] >= 0:
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
n3 = 3
index3 = np.where(nArray==n3)[0][0]
    
K = (1/(4*np.pi))*np.sqrt(3/2)*(cons.h*cons.c/(cons.G*cons.mA**(4/3)))**(3/2)
M = K*x1True[index3]**2*z1True[index3]/(4*cons.mSun)

print('M/Msun = {:.3f}(2/mu_e)**2'.format(M))

# %% Problem 4

n4 = 3/2
index4 = np.where(nArray==n4)[0][0]
#convective with solar composition
mu = 0.61
mu_e = 1.17


rhoCentral = (x1True[index4]/z1True[index4])*(1/(4*np.pi))
rhoCentral_solarUnits = rhoCentral*cons.mSun*cons.rSun**(-3)
print('rhoCentral = {:.3f} M R^-3'.format(rhoCentral))
print('rhoCentral = {:.3f} (M/Msun) (R/Rsun)^-3'.format(rhoCentral_solarUnits))

Tcentral = 2*cons.G/(5*z1True[index4]*x1True[index4]*cons.NA*cons.kB)
Tcentral_solarUnits = Tcentral*cons.mSun*cons.rSun**-1
print('TempCentral = {:.3e} mu M R^-1'.format(Tcentral))
print('TempCentral = {:.3e} M R^-1'.format(Tcentral*mu))
print('TempCentral = {:.3e} mu (M/Msun) (R/Rsun)^-1'.format(Tcentral_solarUnits))

Tcentralb = (1.004e13/(cons.NA*cons.kB))*(x1True[index4]/(4*np.pi*z1True[index4]))**(2/3)
Tcentralb_solarUnits = Tcentralb*cons.mSun**(2/3)*cons.rSun**(-2)
print('When Pgas=Pe, Tcentral = {:.3e} (mu/mu_e^5/3) M^2/3 R^-2'.format(Tcentralb))
print('             or {:.2e} M R^-2'.format(Tcentralb*mu/mu_e**(5/3)))
print('             or {:.3e} (mu/mu_e^5/3) (M/Msun)^2/3 (R/Rsun)^-2'.format(Tcentralb_solarUnits))

Tcentralc = 4e6 #K
Mpart1 = 1.004e13*(cons.NA*cons.kB/mu)*(25/(4*cons.G**2))*(4*np.pi)**(-2/3)
Mpart2 = x1True[index4]**(8/3)*z1True[index4]**(4/3)*Tcentralc*mu_e**(-5/3)
M = (Mpart1*Mpart2)**(3/4)
print('M = {:.3e} g'.format(M))
print('M = {:.3f} Msun'.format(M/cons.mSun))

# %% Problem 5

n5 = 3/2
i5 = np.where(nArray==n5)[0][0]

PpCons = ((8*np.pi*cons.G*cons.sigma)/(3*2.5e-31))**(2/5)*(cons.NA*cons.kB)**(3/5)
print('Pp = {:.3e} (Z/0.02)^-2/5 (M/L)^3/5 mu^-3/5 T_eff^-7/5'.format(PpCons))

KprimePart1 = (cons.NA*cons.kB*x1True[i5])**(5/2)*z1True[i5]**(1/2)
KprimePart2 = (5/(2*cons.G))**3/2*(4*np.pi*cons.sigma)**(3/4)*(4*np.pi)**(-1)
KprimeCons = KprimePart1*KprimePart2
print('Kprime = {:.3e} mu^-5/2 M^-1/2 L^-3/4 T_eff^3/4'.format(KprimeCons))

TeffCons = ((2**(2/3)*PpCons/KprimeCons)**(2/5)*((8/5)**(-2/9)))**(20/32)
print('TeffCons = {:.3e} (Z/0.02)^-50/93 mu^38/93 M^18/93 L^7/93'.format(TeffCons))