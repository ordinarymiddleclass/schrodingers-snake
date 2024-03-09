"""
    hydrogen_radial.py
    ~~~~~~~~~~~~~~~~~~~~~
    This script is used to visualize the radial part of hydrogen atom wavefunctions.
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.special import genlaguerre
from scipy.special import sph_harm
import mpl_toolkits.mplot3d.axes3d as axes3d
import matplotlib.colors as mcolors

def radial_hydrogen_wavefunction(n,l,r):
    """
    This function calculates the radial part of the hydrogen atom wavefunctions.
    """    
    first_part = np.sqrt((2/n)**3 * math.factorial(n-l-1)/(2*n*math.factorial(n+l)))
    second_part = np.exp(-r/n)
    third_part = (2*r/n)**l
    #calculate laguerre polynomial
    fourth_part = genlaguerre(n-l-1, 2*l+1)(2*r/n)
    print(first_part)
    print(second_part)
    print(third_part)
    print(fourth_part)
    return np.sqrt((2/n)**3 * math.factorial(n-l-1)/(2*n*math.factorial(n+l))) * np.exp(-r/n) * (2*r/n)**l * genlaguerre(n-l-1, 2*l+1)(2*r/n)

def radial_hydrogen_wavefunction_plot(n,l,r):
    """
    This function plots the radial part of the hydrogen atom wavefunctions.
    """
    #calculate the radial part of the hydrogen atom wavefunctions
    rsh = radial_hydrogen_wavefunction(n,l,r)
    #plot the radial part of the hydrogen atom wavefunctions
    plt.plot(r,rsh,label=r'$n=$'+str(n)+r', $l=$'+str(l))
    plt.title(r'Radial part of the hydrogen atom wavefunctions, $R_{nl}(r)$')
    plt.xlabel(r'$r$')
    plt.ylabel(r'$R_{nl}(r)$')
    plt.legend()
    plt.show()


n = 2
l = 1 
r = np.linspace(0,30,100)
#plot radial part of the hydrogen atom wavefunctions with multiple n and l on the same plot
r21 = radial_hydrogen_wavefunction(n,l,r)
r31 = radial_hydrogen_wavefunction(3,1,r)
r32 = radial_hydrogen_wavefunction(3,2,r)
#create a new figure
plt.figure()
#plot the radial part of the hydrogen atom wavefunctions
plt.plot(r,r21,label=r'$n=2, l=1$')
plt.plot(r,r31,label=r'$n=3, l=1$')
plt.plot(r,r32,label=r'$n=3, l=2$')
plt.title(r'Radial part of the hydrogen atom wavefunctions, $R_{nl}(r)$')
plt.xlabel(r'$r$')
plt.ylabel(r'$R_{nl}(r)$')
plt.legend()
plt.show()
