import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.special import lpmv
from scipy.special import sph_harm
import mpl_toolkits.mplot3d.axes3d as axes3d
import matplotlib.colors as mcolors

ls = [0,1,2,3]
x=np.linspace(-1,1,100)

plt.figure()

for l in ls:
    plt.plot(x,lpmv(0,l,x),label=r'$l=$'+str(l))
plt.title(r'Legendre Polynomials, $P_l(x)$')
plt.xlabel(r'$x$')
plt.ylabel(r'$P_l(x)$')
plt.legend()
plt.show()

thetas = np.linspace(0,2*np.pi,200)

plt.figure(figsize=(8,8))
for i in range(0, len(ls)):
    l=ls[i]
    r = lpmv(0,l,np.cos(thetas))
    plt.polar(thetas, abs(r),label=r'$l=$'+str(l))
plt.title(r'Associated Legendre Polynomials, $||P_l(x)||$')
plt.ylabel(r'$||P_l(x)||$')
plt.legend()
plt.show()

ls = [0,1,1,1,2]
ms = [0,-1,0,1,1]

x=np.linspace(-1,1,100)

plt.figure()

for i in range(0, len(ls)):
    l=ls[i]
    m=ms[i]
    plt.plot(x,lpmv(m,l,x),label=r'$l=$'+str(l))
plt.title(r'Associated Legendre Polynomials, $P_l^m(x)$')
plt.xlabel(r'$x$')
plt.ylabel(r'$P_l^m(x)$')
plt.legend()
plt.show()

thetas = np.linspace(0,2*np.pi,200)

plt.figure(figsize=(8,8))
for i in range(0, len(ls)):
    l=ls[i]
    m=ms[i]
    r = lpmv(m,l,np.cos(thetas))
    plt.polar(thetas, abs(r),label=r'$l=$'+str(l)+r', $m=$'+str(m))
plt.title(r'Associated Legendre Polynomials, $||P_l^m(x)||$')
plt.ylabel(r'$||P_l^m(x)||$')
plt.legend()
plt.show()

l=2
m=1

thetas = np.linspace(0, np.pi, 20)
phis = np.linspace(0, 2*np.pi, 20)

(Theta,Phi)=np.meshgrid(thetas,phis) 
s_harm=sph_harm(m, l, Phi, Theta)
   
R = np.real(s_harm)
X = R * np.sin(Theta) * np.cos(Phi)
Y = R * np.sin(Theta) * np.sin(Phi)
Z = R * np.cos(Theta)

cmap = plt.get_cmap('jet')
norm = mcolors.Normalize(vmin=Z.min(), vmax=Z.max())

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1,1,1, projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=plt.get_cmap('jet'),facecolors=cmap(norm(R)),
    linewidth=0, antialiased=False, alpha=0.5)
plt.title(r'Spherical Harmonics, $Y_l^m(\theta,\phi)$'+r', $l=$'+str(l)+r', $m=$'+str(m))
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.ylabel(r'$z$')
plt.legend()
plt.show()