#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 21:54:22 2020

@author: abduhu

Complex unitary matrix plotting
********


"""
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import pi, cos, exp, atan
import numpy as np

def cplot(cmatrix, title=''):
    """
    Plots complex matrix using chromatic values.
    """
    dims = np.array(cmatrix).shape
    img = []
    for r, row in enumerate(cmatrix):
        img.append([])
        for c in row:
            y = c.imag
            x = c.real
            
            if x == 0:
                if y > 0:
                    theta = pi/2
                else:
                    theta = -pi/2
            else:
                theta = atan(y/x)
            if x < 0:
                theta += pi
            rad = 1- exp(-(x**2 + y**2)/0.5)
            img[-1].append([cos(theta/2)**2,
                            cos(theta/2 + 2*pi/3)**2,
                            cos(theta/2 - 2*pi/3)**2,
                            rad])
#    plt.imshow(img) #, extent=(1, dims[0], 1, dims[1]))
#    plt.xticks([i for i in range(dims[0])])
#    plt.yticks([dims[1]-1-i for i in range(dims[1])])
#    plt.title(title)
#    plt.axis()
    mpl.rcParams['figure.figsize'] = (10, 10)
    fig, (pl, sc) = plt.subplots(nrows=1, ncols=2, sharex=False)
                                 #figsize=[8, 25])
    sc.imshow(img)
    pl.imshow(scale(), extent=(-1, 1, -1, 1))
    pl.set_xlabel('Re')
    pl.set_ylabel('Img')
    pl.grid(True)
    plt.title(title)
    plt.show()
    return

def scale():
    """
    Plot the scaling spectrum of the complex plane [-1, 1, -i, i]
    """
    img = []
    sc = 50
    for r in range(sc, -sc, -1):
        img.append([])
        for c  in range(-sc, sc, 1):
            y = (r/sc)
            x = (c/sc)
            if x == 0:
                if y > 0:
                    theta = pi/2
                else:
                    theta = -pi/2
            else:
                theta = atan(y/x)
            if x < 0:
                theta += pi
            rad = 1- exp(-(x**2 + y**2)/0.5)
            img[-1].append([cos(theta/2)**2,
                            cos(theta/2 + 2*pi/3)**2,
                            cos(theta/2 - 2*pi/3)**2,
                            rad])
#    plt.imshow(img, extent=(-1, 1, -1, 1))
#    plt.xlabel('Re')
#    plt.ylabel('Img')
#    #plt.title('Scale')
#    plt.show()
    return img
    
    