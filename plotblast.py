#! /usr/bin/env python

"""
Plots results of spiralarm models. Assumes cartesian grid.

Run "plot_spiralarm.py -h" to see description of inputs.

See documentation on athena_read.athdf() for important notes about reading files
with mesh refinement.

The -c <colormap> option is highly recommended to change the default. Consider
"RdBu_r" or "gist_heat" for example.

Users are encouraged to make their own versions of this script for improved
results by adjusting figure size, spacings, tick locations, axes labels, etc.
The script must also be modified to plot any functions of the quantities in the
file, including combinations of multiple quantities.

Requires scipy if making streamplot.
"""

# Python standard modules
import argparse
import warnings

# Python plotting modules
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D

# Other Python modules
import h5py
import numpy as np

# Athena++ modules
import athena_read

# Streamline function
from streamlines import plot_streamlines

# Main function
def map(path,prob,run,out,d,**kwargs):

    data_file  = (("%s%s/%s/%s.out%1i.%05i.athdf") % (path,prob,run,prob,out,d))
    output_file= (("%s%s/plots/map_%s_%s_a%1i_%05i.png") % (path,prob,run,kwargs['quantity'],kwargs['axis'],d))
    print('[plotraw]: reading file %s' % (data_file))

    inputs     = athena_read.athinput("%s%s/%s/athinput.%s" % (path,prob,run,prob))
    gamma      = (inputs['hydro'])['gamma']

    # Determine refinement level to use
    if kwargs['level'] is not None:
        level = kwargs['level']
    else:
        level = None

    vector = (kwargs['vector'] is not None) 

    if vector:
        if kwargs['vector'] == 'vel':
            vecstr = 'mom'
        else:
           vecstr = kwargs['vector']

    # Extract basic coordinate information
    with h5py.File(data_file, 'r') as f:
        coordinates = f.attrs['Coordinates'].decode("utf-8")
        Time        = f.attrs['Time']
    if (coordinates != "cartesian"):
        raise Exception("[plotraw]: Coordinates required: cartesian. Found: %s" % (coordinates))
    data = athena_read.athdf(data_file, face_func_1=None)
    dim  = np.array(data['dens'].shape)
    x     = data['x1v']
    y     = data['x2v']
    z     = data['x3v']
    x1f   = data['x1f']
    x2f   = data['x2f']
    x3f   = data['x3f']    
    
    if 'Eint' in data:
        prss  = data['Eint']*(gamma-1)
    else:
        prss  = data['Etot']-0.5*(data['mom1']**2+data['mom2']**2+data['mom3']**2)/data['dens']
        if 'Bcc1' in data:
            prss -= 0.5*(data['Bcc1']**2+data['Bcc2']**2+data['Bcc3']**2)
        prss *= (gamma-1)

    if kwargs['quantity'] == 'temp':
        dq = prss/data['dens']
    elif (kwargs['quantity'])[0:3] == 'vel':
        dq = data[("mom%1s" % ((kwargs['quantity'])[3]))]/data['dens']
    elif kwargs['quantity'] == 'prss':
        dq = prss 
    elif kwargs['quantity'] == 'mach':
        cs   = np.sqrt(gamma*prss/data['dens'])
        dq   = np.sqrt((data['mom1']/data['dens'])**2 + (data['mom2']/data['dens'])**2 + (data['mom3']/data['dens'])**2)/cs
    elif kwargs['quantity'] == 'c0':
        dq   = data['s0']/data['dens']
    elif kwargs['quantity'] == 'vabs':
        dq   = np.sqrt((data['mom1']/data['dens'])**2 + (data['mom2']/data['dens'])**2 + (data['mom3']/data['dens'])**2)
    elif kwargs['quantity'] == 'grdp':
        if dim[0] == 1:
            g    = np.gradient(np.squeeze(prss))  
            X,Y  = np.meshgrid(x,y)
            dq   = (g[0]*Y+g[1]*X)/np.sqrt(X**2+Y**2)
            dq   = np.abs(dq)
        else:
            g     = np.gradient(np.squeeze(prss))
            X,Y,Z = np.meshgrid(x,y,z)
            dq    = (g[0]*Z+g[1]*Y+g[2]*X)/np.sqrt(X**2+Y**2+Z**2)
            dq    = np.abs(dq)
    elif kwargs['quantity'] == 'pmag':
        if 'Bcc1' in data:
            dq   = 0.5*(data['Bcc1']**2+data['Bcc2']**2+data['Bcc3']**2)
        else:
            raise Exception('[plotblast]: no magnetic fields in data set for quantity %s' % (kwargs['quantity']))
    elif kwargs['quantity'] == 'grdm':
        if 'Bcc1' in data:
            pmag   = 0.5*(data['Bcc1']**2+data['Bcc2']**2+data['Bcc3']**2)
            if dim[0] == 1:
                g    = np.gradient(np.squeeze(pmag))
                X,Y  = np.meshgrid(x,y)
                dq   = (g[0]*Y+g[1]*X)/np.sqrt(X**2+Y**2)
                dq   = np.abs(dq)
            else:
                g     = np.gradient(np.squeeze(pmag))
                X,Y,Z = np.meshgrid(x,y,z)
                dq    = (g[0]*Z+g[1]*Y+g[2]*X)/np.sqrt(X**2+Y**2+Z**2)
                dq    = np.abs(dq)
        else:
            raise Exception('[plotblast]: no magnetic fields in data set for quantity %s' % (kwargs['quantity']))
    elif kwargs['quantity'] == 'divb':
        if 'Bcc1' in data:
            if dim[0] == 1:
                gx   = np.gradient(np.squeeze(data['Bcc1']),axis=1)
                gy   = np.gradient(np.squeeze(data['Bcc2']),axis=0)
                dq   = gx+gy
            else:
                gx   = np.gradient(np.squeeze(data['Bcc1']),axis=2)
                gy   = np.gradient(np.squeeze(data['Bcc2']),axis=1)
                gz   = np.gradient(np.squeeze(data['Bcc3']),axis=0)
                dq   = gx+gy+gz
        else:
            raise Exception('[plotblast]: no magnetic fields in data set for quantity %s' % (kwargs['quantity']))
    elif kwargs['quantity'] == 'ptot':
        if 'Bcc1' in data:
            dq   = prss + 0.5*(data['Bcc1']**2+data['Bcc2']**2+data['Bcc3']**2)
        else:
            dq   = prss
    elif kwargs['quantity'] == 'beta':
        if 'Bcc1' in data:
            dq = prss/(0.5*(data['Bcc1']**2+data['Bcc2']**2+data['Bcc3']**2))
        else:
            raise Exception('[plotraw]: no magnetic fields in data set for quantity %s' % (kwargs['quantity']))
    else:
        dq    = data[kwargs['quantity']]

    # dimensions
    if dim[0]==1:
        dq = np.squeeze(dq)
        if kwargs['center']:
            dq = dq[:,len(x)//4:3*len(x)//4]
            x  = x[len(x)//4:3*len(x)//4]
        xplot = x
        yplot = y
        lplot = z
        xstr  = 'x'
        ystr  = 'y'
        if vector:
            vecx = np.squeeze(data[vecstr+'1'])
            vecy = np.squeeze(data[vecstr+'2'])
    else: # 3d data. Transpose 3D volume such that last axis can be summed over.
        if kwargs['axis'] == 0:
            xplot   = y
            yplot   = z
            lplot   = x
            nplot   = np.array([dim[1],dim[0]])
            xstr    = 'y'
            ystr    = 'z'
            dd      = data['dens']
            dplot   = np.array(dim[[0,1,2]]) # y,z,x
            if vector:
                vecx = data[vecstr+'2']
                vecy = data[vecstr+'3']
        elif kwargs['axis'] == 1:
            xplot   = x
            yplot   = z
            lplot   = y
            nplot   = np.array([dim[2],dim[0]])
            xstr    = 'x'
            ystr    = 'z'
            dq      = np.transpose(dq,(0,2,1))
            dd      = np.transpose(data['dens'],(0,2,1))
            dplot   = np.array(dim[[0,2,1]])
            if vector:
                vecx = np.transpose(data[vecstr+'1'],(0,2,1))
                vecy = np.transpose(data[vecstr+'3'],(0,2,1))
        elif kwargs['axis'] == 2:
            xplot   = x
            yplot   = y
            lplot   = z
            nplot   = np.array([dim[2],dim[1]])
            xstr    = 'x'
            ystr    = 'y'
            dq      = np.transpose(dq,(1,2,0))
            dd      = np.transpose(data['dens'],(1,2,0))
            dplot   = np.array(dim[[1,2,0]])
            if vector:
                vecx = np.transpose(data[vecstr+'1'],(1,2,0))
                vecy = np.transpose(data[vecstr+'2'],(1,2,0))
        else: 
            print('[plotblast]: invalid axis value %i' % (kwargs['axis']))
            quit()

        if vector:
            if vecstr == 'mom':
                vecx /= dd
                vecy /= dd

        if kwargs['midplane'] and kwargs['surface']:
            print('[plotblast]: keywords midplane and surface are exclusive')
            quit()
        if kwargs['midplane']:
            if dplot[2] % 2 == 1:
                dq      = dq[:,:,dplot[2]//2]
                if vector:
                    vecx = vecx[:,:,dplot[2]//2]
                    vecy = vecy[:,:,dplot[2]//2]
            else:
                dq      = 0.5*(dq[:,:,dplot[2]//2-1]+dq[:,:,dplot[2]//2])
                if vector:
                    vecx = 0.5*(vecx[:,:,dplot[2]//2-1]+vecx[:,:,dplot[2]//2])
                    vecy = 0.5*(vecy[:,:,dplot[2]//2-1]+vecy[:,:,dplot[2]//2])
        if kwargs['surface']:
            if (kwargs['quantity'])[0:3] == 'vel':
                dq = np.sum(dq*dd,axis=2)/np.sum(dd,axis=2)
            elif kwargs['quantity'] == 'temp':
                dq = np.sum(dq*dd,axis=2)/np.sum(dd,axis=2)
            elif kwargs['quantity'] == 'c0':
                dq = np.sum(dq*dd,axis=2)/np.sum(dd,axis=2)
            else:
                dq = np.mean(dq,axis=2)#*(lplot[1]-lplot[0])
            if vector:
                vecx = np.sum(vecx*dd,axis=2)/np.sum(dd,axis=2)
                vecy = np.sum(vecy*dd,axis=2)/np.sum(dd,axis=2)


    if kwargs['units']:
        Time     = Time * 9.52517888e+01
        xplot    = xplot    * 8.87652670
        yplot    = yplot    * 8.87652670
        lplot    = lplot    * 8.87652670
        ulengstr = ' [pc]'
        utimestr = ' [Myr]' 
        if kwargs['quantity'] == 'dens':
            quantlabel = r'n [cm$^{-3}$]'
        elif (kwargs['quantity'])[0:3] == 'vel':
            quantlabel = r'v [km s$^{-1}$]'
            dq = dq * 9.11836593e-02
        elif kwargs['quantity'] == 'temp':
            quantlabel = 'T [K]'
        elif kwargs['quantity'] == 'c0':
            quantlabel = 'color'
        else:
            quantlabel = 'unknown'
    else:
        ulengstr = ''
        utimestr = ''
        quantlabel = kwargs['quantity']
    nx1 = len(x)
    nx2 = len(y)
    nx3 = len(z)
    if kwargs['units']:
        print('[plotraw]: time        = %13.5e%s' % (Time,utimestr))
    else: 
        print('[plotraw]: time        = %13.5e' % (Time))
    print('[plotraw]: dimensions  = %5i %5i %5i minmax(x)=%13.5e %13.5e minmax(y)=%13.5e %13.5e minmax(z)=%13.5e %13.5e' 
          % (nx1,nx2,nx3,np.min(x1f),np.max(x1f),np.min(x2f),np.max(x2f),np.min(x3f),np.max(x3f)))
    print('[plotraw]: min/max(dq) = %13.5e %13.5e' % (np.min(dq),np.max(dq)))
  
    x_min = np.min(xplot)
    x_max = np.max(xplot)
    y_min = np.min(yplot)
    y_max = np.max(yplot)
    z_min = np.min(lplot)
    z_max = np.max(lplot)

    r_max = x_max

    x_grid,y_grid = np.meshgrid(xplot,yplot)

    vals = np.squeeze(dq)
    if kwargs['logc']:
        quantlabel = 'log '+quantlabel
        vals = np.log10(vals)
    print(vals[0,0],vals[0,1],vals[1,0],vals[1,1])
  
    # Determine colormapping properties
    cmap = plt.get_cmap(kwargs['colormap'])
    vmin = kwargs['vmin']
    vmax = kwargs['vmax']
    #if kwargs['logc']:
    #  norm = colors.LogNorm()
    #else:
    #  norm = colors.Normalize()
  
    # Make plot
    xsize = 7.5
    ysize = 6
    xlo   = 0.1
    ylo   = 0.1
    yup   = 0.9
    xup   = xlo+(yup-ylo)*ysize/xsize
    fig = plt.figure(num=0,figsize=(xsize,ysize),facecolor='white',dpi=100)
    ax  = fig.add_axes([xlo,ylo,(xup-xlo),(yup-ylo)])
    pos = ax.get_position()
    if kwargs['contour'] > 0:
        im  = ax.contour(x_grid,y_grid,vals,vmin=vmin,vmax=vmax,levels=kwargs['contour'],cmap=cmap,zorder=1)
    else:
        im  = ax.pcolormesh(x_grid, y_grid, vals, cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')
    if vector or stream:
        skip = kwargs['samples']
        #vim = ax.quiver(x_grid[skip//2::skip,skip//2::skip],y_grid[skip//2::skip,skip//2::skip],
        #                vecx[skip//2::skip,skip//2::skip],vecy[skip//2::skip,skip//2::skip],color='white',pivot='middle')  
        print(x_grid.shape)
        sp  = np.zeros((2*len(xplot)//skip,2))
        off = len(xplot)//skip 
        sp[:off,0] = xplot[skip//2]
        sp[:off,1] = yplot[skip//2::skip]
        sp[off:,0] = xplot[skip//2::skip]
        sp[off:,1] = yplot[skip//2] 
        if kwargs['vector'] == 'Bcc':
            vim = ax.streamplot(x_grid,y_grid,vecx,vecy,color=kwargs['cvmap'],start_points=sp,density=30)
        elif kwargs['vector'] == 'vel':
            vim = ax.quiver(x_grid[skip//2::skip,skip//2::skip],y_grid[skip//2::skip,skip//2::skip],
                            vecx[skip//2::skip,skip//2::skip],vecy[skip//2::skip,skip//2::skip],color='white',pivot='middle')  

    #fig.gca().set_aspect('equal')
    ax.set_xlim((x_min, x_max))
    ax.set_ylim((y_min, y_max))
    x_string = xstr
    y_string = ystr
    if kwargs['grid']:
        ax.grid()
    ax.set_xlabel('$'+x_string+'$'+ulengstr)
    ax.set_ylabel('$'+y_string+'$'+ulengstr)
    if kwargs['title']:
        if kwargs['units']:
            ax.set_title('Time = %4.0f%s' % (Time,utimestr))
        else:
            ax.set_title('Time = %13.5e' % Time)
    ax.set_position(pos)
    if 'circ' in kwargs:
        if kwargs['circ'] > 0.0:
            circle1 = plt.Circle((0.0,0.0), kwargs['circ'], color='red',fill=False,zorder=2)
            ax.add_patch(circle1)
    if kwargs['label'] != 'None':
        ax.text(0.05,0.95,kwargs['label'],transform=ax.transAxes,fontsize=12)
    #cb = plt.colorbar(im,ax=ax,orientation="vertical")
    #cb.set_label(quantlabel)

    cax = fig.add_axes([0.8, ylo, 0.05, yup-ylo])
    cb = fig.colorbar(im, orientation='vertical', cax=cax)
    cb.set_label(quantlabel)



    #fig.tight_layout()
    if kwargs['save']:
        print('[plotblast]: Saving file %s' % (output_file))
        fig.savefig(output_file,format='png',dpi=200)
        plt.close() 
    else:
        plt.show()










#3d voxel projection

def voxel(path,prob,run,out,d,**kwargs):
    data_file  = (("%s%s/%s/%s.out%1i.%05i.athdf") % (path,prob,run,prob,out,d))
    output_file= (("%s%s/plots/map_%s_%s_a%1i_%05i.png") % (path,prob,run,kwargs['quantity'],kwargs['axis'],d))
    print('[plotraw]: reading file %s' % (data_file))

    inputs     = athena_read.athinput("%s%s/%s/athinput.%s" % (path,prob,run,prob))
    gamma      = (inputs['hydro'])['gamma']

    # Determine refinement level to use
    if kwargs['level'] is not None:
        level = kwargs['level']
    else:
        level = None

    vector = (kwargs['vector'] is not None) 

    if vector:
        if kwargs['vector'] == 'vel':
            vecstr = 'mom'
        else:
           vecstr = kwargs['vector']

    # Extract basic coordinate information
    with h5py.File(data_file, 'r') as f:
        coordinates = f.attrs['Coordinates'].decode("utf-8")
        Time        = f.attrs['Time']
    if (coordinates != "cartesian"):
        raise Exception("[plotraw]: Coordinates required: cartesian. Found: %s" % (coordinates))
    data = athena_read.athdf(data_file, face_func_1=None)
    dim  = np.array(data['dens'].shape)
    x     = data['x1v']
    y     = data['x2v']
    z     = data['x3v']
    x1f   = data['x1f']
    x2f   = data['x2f']
    x3f   = data['x3f']    
    
    if 'Eint' in data:
        prss  = data['Eint']*(gamma-1)
    else:
        prss  = data['Etot']-0.5*(data['mom1']**2+data['mom2']**2+data['mom3']**2)/data['dens']
        if 'Bcc1' in data:
            prss -= 0.5*(data['Bcc1']**2+data['Bcc2']**2+data['Bcc3']**2)
        prss *= (gamma-1)

    if kwargs['quantity'] == 'temp':
        dq = prss/data['dens']
    elif (kwargs['quantity'])[0:3] == 'vel':
        dq = data[("mom%1s" % ((kwargs['quantity'])[3]))]/data['dens']
    elif kwargs['quantity'] == 'prss':
        dq = prss 
    elif kwargs['quantity'] == 'mach':
        cs   = np.sqrt(gamma*prss/data['dens'])
        dq   = np.sqrt((data['mom1']/data['dens'])**2 + (data['mom2']/data['dens'])**2 + (data['mom3']/data['dens'])**2)/cs
    elif kwargs['quantity'] == 'c0':
        dq   = data['s0']/data['dens']
    elif kwargs['quantity'] == 'vabs':
        dq   = np.sqrt((data['mom1']/data['dens'])**2 + (data['mom2']/data['dens'])**2 + (data['mom3']/data['dens'])**2)
    elif kwargs['quantity'] == 'grdp':
        if dim[0] == 1:
            g    = np.gradient(np.squeeze(prss))  
            X,Y  = np.meshgrid(x,y)
            dq   = (g[0]*Y+g[1]*X)/np.sqrt(X**2+Y**2)
            dq   = np.abs(dq)
        else:
            g     = np.gradient(np.squeeze(prss))
            X,Y,Z = np.meshgrid(x,y,z)
            dq    = (g[0]*Z+g[1]*Y+g[2]*X)/np.sqrt(X**2+Y**2+Z**2)
            dq    = np.abs(dq)
    elif kwargs['quantity'] == 'pmag':
        if 'Bcc1' in data:
            dq   = 0.5*(data['Bcc1']**2+data['Bcc2']**2+data['Bcc3']**2)
        else:
            raise Exception('[plotblast]: no magnetic fields in data set for quantity %s' % (kwargs['quantity']))
    elif kwargs['quantity'] == 'grdm':
        if 'Bcc1' in data:
            pmag   = 0.5*(data['Bcc1']**2+data['Bcc2']**2+data['Bcc3']**2)
            if dim[0] == 1:
                g    = np.gradient(np.squeeze(pmag))
                X,Y  = np.meshgrid(x,y)
                dq   = (g[0]*Y+g[1]*X)/np.sqrt(X**2+Y**2)
                dq   = np.abs(dq)
            else:
                g     = np.gradient(np.squeeze(pmag))
                X,Y,Z = np.meshgrid(x,y,z)
                dq    = (g[0]*Z+g[1]*Y+g[2]*X)/np.sqrt(X**2+Y**2+Z**2)
                dq    = np.abs(dq)
        else:
            raise Exception('[plotblast]: no magnetic fields in data set for quantity %s' % (kwargs['quantity']))
    elif kwargs['quantity'] == 'divb':
        if 'Bcc1' in data:
            if dim[0] == 1:
                gx   = np.gradient(np.squeeze(data['Bcc1']),axis=1)
                gy   = np.gradient(np.squeeze(data['Bcc2']),axis=0)
                dq   = gx+gy
            else:
                gx   = np.gradient(np.squeeze(data['Bcc1']),axis=2)
                gy   = np.gradient(np.squeeze(data['Bcc2']),axis=1)
                gz   = np.gradient(np.squeeze(data['Bcc3']),axis=0)
                dq   = gx+gy+gz
        else:
            raise Exception('[plotblast]: no magnetic fields in data set for quantity %s' % (kwargs['quantity']))
    elif kwargs['quantity'] == 'ptot':
        if 'Bcc1' in data:
            dq   = prss + 0.5*(data['Bcc1']**2+data['Bcc2']**2+data['Bcc3']**2)
        else:
            dq   = prss
    elif kwargs['quantity'] == 'beta':
        if 'Bcc1' in data:
            dq = prss/(0.5*(data['Bcc1']**2+data['Bcc2']**2+data['Bcc3']**2))
        else:
            raise Exception('[plotraw]: no magnetic fields in data set for quantity %s' % (kwargs['quantity']))
    else:
        dq    = data[kwargs['quantity']]

    # Condense arrays

    def shrink(a):
        a = .125*(       a[ ::2, ::2, ::2]
                       + a[1::2, ::2, ::2]
                       + a[ ::2, 1::2, ::2]
                       + a[1::2, 1::2, ::2]
                       + a[ ::2, ::2, 1::2]
                       + a[1::2, ::2, 1::2]
                       + a[ ::2, 1::2, 1::2]
                       + a[1::2, 1::2, 1::2])
        return a
    
    dq_new = shrink(dq)
    dq_new = shrink(dq_new)
    dq_new = dq_new[16:, 16:, 16:]

    x1f_new = x1f[::2]
    x1f_new = x1f_new[::2]
    x1f_new = x1f_new[16:]
    x2f_new = x2f[::2]
    x2f_new = x2f_new[::2]
    x2f_new = x2f_new[16:]
    x3f_new = x3f[::2]
    x3f_new = x3f_new[::2]
    x3f_new = x3f_new[16:]

    x_new = x1f[1::2]
    x_new = x_new[1::2]
    x_new = x_new[16:]
    y_new = x2f[1::2]
    y_new = y_new[1::2]
    y_new = y_new[16:]
    z_new = x3f[1::2]
    z_new = z_new[1::2]
    z_new = z_new[16:]

    # Create rgb values
    def midpoints(x):
        sl = ()
        for _ in range(x.ndim):
            x = (x[sl + np.index_exp[:-1]] + x[sl + np.index_exp[1:]]) / 2.0
            sl += np.index_exp[:]
        return x

    xr = (x1f_new - x1f_new[0]) / (x1f_new[-1] - x1f_new[0])
    yg = (x2f_new - x2f_new[0]) / (x2f_new[-1] - x2f_new[0])
    zb = (x3f_new - x3f_new[0]) / (x3f_new[-1] - x3f_new[0])

    r, g, b = np.meshgrid(xr, yg, zb)
    
    
    #r, g, b = np.meshgrid((x1f_new - x1f_new[0]) / (x1f_new[-1] - x1f_new[0]),
                        #(x2f_new - x2f_new[0]) / (x2f_new[-1] - x2f_new[0]),
                        #(x3f_new - x3f_new[0]) / (x3f_new[-1] - x3f_new[0]))
    
    rc = midpoints(r)
    gc = midpoints(g)
    bc = midpoints(b)

    # define a sphere
    sphere = dq_new > 1

    # combine the color components
    colors = np.zeros(sphere.shape + (4,))
    colors[..., 0] = rc
    colors[..., 1] = gc
    colors[..., 2] = bc
    colors[..., 3] = 1*((dq_new - np.min(dq_new)) / (np.max(dq_new) - np.min(dq_new)))

    # Plot
    ax = plt.figure().add_subplot(projection='3d')
    ax.voxels(r, g, b, sphere, facecolors=colors, edgecolors=np.clip(2*colors - 0.5, 0, 1), linewidth=0.5)
    ax.set(xlabel='x', ylabel='y', zlabel='z')
    ax.set_aspect('equal')

    if ("Bcc1" in data):
        u = data["Bcc1"]
        u = shrink(u)
        u = shrink(u)
        u = u[16:, 16:, 16:]
        u = u*sphere
        
        v = data["Bcc2"]
        v = shrink(v)
        v = shrink(v)
        v = v[16:, 16:, 16:]
        v = v*sphere
        
        w = data["Bcc3"]
        w = shrink(w)
        w = shrink(w)
        w = w[16:, 16:, 16:]
        w = w*sphere
        
        stride = 4
        rc_stride = rc[::stride, ::stride, ::stride]
        gc_stride = gc[::stride, ::stride, ::stride]
        bc_stride = bc[::stride, ::stride, ::stride]
        u_stride = u[::stride, ::stride, ::stride]
        v_stride = v[::stride, ::stride, ::stride]
        w_stride = w[::stride, ::stride, ::stride]

        c = (0, 1, 0, .5)
        colors = [c, c, c]
        
        #ax.quiver(rc_stride, gc_stride, bc_stride, u_stride, v_stride, w_stride, color= colors, length = .1, normalize=False, pivot = "middle")

    u = data["Bcc1"]
    v = data["Bcc2"]
    w = data["Bcc3"]

    #random starting pts, inside vs outside sampling, colors, 64 cubed figure
    """
    s = 0
    for i in range(4):

        q = 0
        r = 0
        for i in range(5):
            plot_streamlines(q, r, s, x, y, z, u, v, w, 0)
            q += .175

        q = 0
        r = 0
        for i in range(5):
            plot_streamlines(q, r, s, x, y, z, u, v, w, 0)
            r += .175
            
        s += .2

    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.set_zlim([0, .75])
    """
    plt.show()
    
    #stop streamlines once slope 0

    
