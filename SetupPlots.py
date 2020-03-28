#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 17:24:55 2019

@author: Sandra Bustamante
"""
import numpy as np
import matplotlib.pyplot as plt

def setupPlot(fig, figsize = (1, 1), grid=(1,1)):
    '''Configure the size of the plots.'''

    nwidth, nheight = figsize
    width = 3.39*nwidth
    height = nheight*width*(np.sqrt(5.)-1.)/2.
    fontsize = 8
    linewidth = 0.4
    markersize = 1

    params = {'axes.labelsize': fontsize,
              'axes.titlesize': fontsize,
              'font.size': fontsize,
              'legend.fontsize': fontsize-2,
              'xtick.labelsize': fontsize,
              'ytick.labelsize': fontsize,
              'lines.linewidth': linewidth,
              'grid.linewidth': linewidth*.7,
              'axes.axisbelow': True,
              'pgf.rcfonts': False,
              'lines.markersize': markersize,
              }
    plt.rcParams.update(params)
    fig.set_size_inches((width, height))
    gs = fig.add_gridspec(2, 2, hspace=0, width_ratios=[10., .5])

    return gs
