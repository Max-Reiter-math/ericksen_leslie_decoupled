# my plots
from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
import os


def sliced_quiver(mesh, f,  id, t, path, z=0.0, tol= 1e-10, scale = "auto", abs=True):
    """
    scale is either a string "auto" or "unit" or a tuple
    """
    plt.close("all")
    os.makedirs(path+"/plots",exist_ok=True)
    #-
    vertices = mesh.coordinates()
    slice_vertices_x = []
    slice_vertices_y = []
    X = []
    Y = []
    Z = []
    #-
    if abs: transf = lambda o : np.abs(o)
    else: transf = lambda o : o
    for x in vertices:
        if (x[2]-z)**2 <= tol:
            slice_vertices_x.append(x[0])
            slice_vertices_y.append(x[1])
            tmp = f(x)
            X.append(tmp[0])
            Y.append(tmp[1])
            Z.append(transf(tmp[2]))
    plt.clf()
    if Z == []:
        # this means there are no vertices in that plane
        raise Exception(" Error: no vertices in the plane z="+str(z))
    else:
        if scale == "auto":
            cnorm = mpl.colors.Normalize(vmin=0.0,vmax=np.max(Z))
        elif scale == "unit":
            cnorm = mpl.colors.Normalize(vmin=0.0,vmax=1.0)
        else:
            cnorm = mpl.colors.Normalize(vmin=scale[0],vmax=scale[1])
        cmap = mpl.cm.winter
        sm = plt.cm.ScalarMappable(cmap=cmap,norm=cnorm)
        plt.quiver(slice_vertices_x, slice_vertices_y, X,Y, color= cmap(cnorm(Z)))
        plt.colorbar(sm)
        plt.savefig(path+"/plots/"+id+"_"+str(t).replace(".","-")+".png", dpi=400)

def quiver_2d(mesh, f, id, t, path, scale="auto"):
    plt.close("all")
    os.makedirs(path+"/plots",exist_ok=True)
    #-
    vertices = mesh.coordinates()
    slice_vertices_x = []
    slice_vertices_y = []
    X = []
    Y = []
    magnitudes = []
    for x in vertices:
        slice_vertices_x.append(x[0])
        slice_vertices_y.append(x[1])
        tmp = f(x)
        X.append(tmp[0])
        Y.append(tmp[1])  
        magnitudes.append(np.sqrt(tmp[0]**2 + tmp[1]**2))
    plt.clf()
    if scale == "auto":
        cnorm = mpl.colors.Normalize(vmin=0.0,vmax=np.max(magnitudes))
    elif scale == "unit":
        cnorm = mpl.colors.Normalize(vmin=0.0,vmax=1.0)
    else:
        cnorm = mpl.colors.Normalize(vmin=scale[0],vmax=scale[1])
    cmap = mpl.cm.winter_r
    sm = plt.cm.ScalarMappable(cmap=cmap,norm=cnorm)
    
    plt.quiver(slice_vertices_x, slice_vertices_y, X,Y, color= cmap(cnorm(magnitudes)))
    plt.colorbar(sm)
    plt.savefig(path+"/plots/"+id+"_"+str(t).replace(".","-")+".png", dpi=400)

def energy_plot(x_data, Energy_data, labels, path):
    """
    takes the energy data as list of lists, where there is a list of all energy types at every time step
    """
    plt.close("all")
    os.makedirs(path+"/plots",exist_ok=True)
    plt.clf()
    plt.ylabel("Energy")
    plt.xlabel("Time")
    time_data = x_data
    cmap = mpl.cm.get_cmap('Spectral')
    Energy_data = list(zip(*Energy_data))
    Energy_data = [list(a) for a in Energy_data]
    n = len(Energy_data)
    for i in range(len(Energy_data)):
        plt.plot(time_data, Energy_data[i], label = labels[i], color = cmap(i/n))
    plt.legend(loc = "upper right")
    plt.savefig(path+"/plots/energy_plot.png", dpi=400)
