import numpy
import math
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from scipy.special import sph_harm
import scipy.special
from scipy.special import assoc_laguerre
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import ListedColormap
from skimage import measure

def hydrogen_wf(n,l,m,X,Y,Z):
    R = numpy.sqrt(X**2+Y**2+Z**2)
    Theta = numpy.arccos(Z/R)
    #Phi = numpy.arctan(Y/X)
    Phi = numpy.arctan2(Y,X)
    
    rho = 2.*R/n
    s_harm=sph_harm(m, l, Phi, Theta)
    #s_harm=sph_harm(m, l, Theta, Phi)
    #l_poly = assoc_laguerre(rho, n-l-1, 2*l+1)
    l_poly = scipy.special.genlaguerre(n-l-1,2*l+1)(rho)
    #l_poly = scipy.special.genlaguerre(n+l,2*l+1)(rho)
    
    prefactor = numpy.sqrt((2./n)**3*math.factorial(n-l-1)/(2.*n*math.factorial(n+l)))
    #prefactor=1
    wf = prefactor*numpy.exp(-rho/2.)*rho**l*s_harm*l_poly
    wf = numpy.nan_to_num(wf)
    return wf

dz=0.05
zmin=-2
zmax=2
x = numpy.arange(zmin,zmax,dz)
y = numpy.arange(zmin,zmax,dz)
z = numpy.arange(zmin,zmax,dz)
X,Y,Z = numpy.meshgrid(x,y,z)

#Change these to change which orbital to plot
n=2
l=0
m=0

data = hydrogen_wf(n,l,m,X,Y,Z)
data = abs(data)**2
#data = data.real
print(sum(sum(sum(data)))*dz**3)
R = numpy.sqrt(X**2+Y**2+Z**2)
#print("Theta")
#print(numpy.arccos(Z/R))
#print("Phi")
#print(numpy.arctan2(Y,X))

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.15, bottom=0.15)
#im = plt.imshow(data[:,:,0], vmin=0, vmax=0.05,extent=[zmin,zmax,zmin,zmax])
im = plt.imshow(data[:,:,int((0-zmin)/dz)], vmin=0, vmax = numpy.max(data), extent=[zmin,zmax,zmin,zmax])
#im = plt.contour(data[:,:,0])
plt.colorbar()
sli = Slider(plt.axes([0.25, 0.025, 0.65, 0.03]), "Z", z[0], z[len(z)-1], valinit=0)

def update(val):
    index = int((sli.val-zmin) / dz)
    im.set_data(data[:,:,index])
       
sli.on_changed(update)
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim([0,len(x)])
ax.set_ylim([0,len(y)])
ax.set_zlim([0,len(z)])
max_val = numpy.max(data)
#ax.plot_surface(X,Y,zval+numpy.zeros_like(X), rstride=1, cstride=1, facecolors = my_cmap(colors))
verts, faces, nors, vas = measure.marching_cubes(data, max_val/2, spacing = (1,1,1))
result=ax.plot_trisurf(verts[:,0], verts[:,1], faces, verts[:,2], cmap ='Spectral', lw=0)
ax.set_title("Hydrogen Orbital Isosurface - n="+str(n)+", l="+str(l)+", m="+str(m))

sli = Slider(plt.axes([0.25, 0.025, 0.65, 0.03]), "iso", 0, max_val, valinit=max_val/2)

def update(val):
    ax.clear()
    verts, faces = measure.marching_cubes(data, sli.val, spacing = (1,1,1))
    result = ax.plot_trisurf(verts[:,0], verts[:,1], faces, verts[:,2], cmap ='Spectral', lw=0)
    ax.set_xlim([0,len(x)])
    ax.set_ylim([0,len(y)])
    ax.set_zlim([0,len(z)])
    ax.set_title("Hydrogen Orbital Isosurface - n="+str(n)+", l="+str(l)+", m="+str(m))
           
sli.on_changed(update)
plt.show()