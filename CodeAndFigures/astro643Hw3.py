# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 15:39:04 2020

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
import pandas as pd


#%% Problem 2

def RK4(dvdt, a, b, h, IV, dim=2, *args):
    """Integrate 1st Order Diferrential Equations using RK4.

    Parameters
    ----------
    dvdt : TYPE
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

nArray = (0, 1, 1.5, 2, 3, 4)
x_a = 0
x_b = 15
h = .001 
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
    x1 = np.interp(0, y[index-1:index], x[index-1:index])
    z1 = np.interp(0, y[index-1:index], z[index-1:index])
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

#create dataframe
df = pd.DataFrame(data, columns=('n', 'x1', '-z1','ratio'))
x1True = np.array((np.sqrt(6),np.pi, 3.6538, 4.3529, 6.8969, 14.972))
z1True = np.array((np.sqrt(6)/3, 1/np.pi, 0.20330, 0.12725, 0.04243, 0.00802))
df.insert(loc=2, column='xTrue', value=x1True)
df.insert(loc=4, column='zTrue', value=z1True)
print(df)

# Create the file.tex
df.rename(columns = {'n':r'n',
                     'x1':r'$\xi_1$',
                     'xTrue':r'$\xi_1$ True', 
                     '-z1':r'$-\theta_n^\prime(\xi_1)$',
                     'zTrue':r'$-\theta_n^\prime(\xi_1)$ True',
                     'ratio':r'$\frac{1}{3}\left(\frac{\xi}{-\theta_n^\prime}\right)_\xi_1'
                     }, inplace = True)

with open('Astro643Hw3Problem2.tex', 'w') as tf:
    tf.write(df.to_latex(float_format='%.3f',
                         index=False,
                         escape=False))





