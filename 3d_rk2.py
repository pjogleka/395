# Implementation of 2nd order Runge-Kutta integration in 3 dimensions.

import matplotlib.pyplot as plt
import numpy as np

def find_nearest(x, y, z):

    x0 = np.floor(((x-np.min(a))/(np.max(a)-np.min(a))*(len(a)-1))).astype(int)
    x1 = x0 + 1

    y0 = np.floor(((y-np.min(b))/(np.max(b)-np.min(b))*(len(b)-1))).astype(int)
    y1 = y0 + 1

    z0 = np.floor(((z-np.min(c))/(np.max(c)-np.min(c))*(len(c)-1))).astype(int)
    z1 = z0 + 1
    #print(x0, y0, z0)
    return x0, x1, y0, y1, z0, z1


def interpolate(x, x0, x1, y, y0, y1, z, z0, z1, c000, c100, c001, c101, c010, c110, c011, c111):
    #print(x0, x1, y0, y1, z0, z1)
    xd = (x-x0)/(x1-x0)
    yd = (y-y0)/(y1-y0)
    zd = (z-y0)/(z1-z0)

    c00 = c000*(1-xd)+c100*xd
    c01 = c001*(1-xd)+c101*xd
    c10 = c010*(1-xd)+c110*xd
    c11 = c011*(1-xd)+c111*xd

    c0 = c00*(1-yd)+c10*yd
    c1 = c01*(1-yd)+c11*yd

    c = c0*(1-zd)+c1*zd

    return c


def plot_streamlines(first_x, first_y, first_z):

    xs = [first_x]
    ys = [first_y]
    zs = [first_z]
    k2x = first_x
    k2y = first_y
    k2z = first_z
    iters = 0

    while k2x < 1 and k2x > -1 and k2y < 1 and k2y > -1 and k2z < 1 and k2z > -1 and iters < 100000:

        dt = .01
        
        x0, x1, y0, y1, z0, z1 = find_nearest(xs[-1], ys[-1], zs[-1])
        
        u_int1 = interpolate(xs[-1], a[x0], a[x1], ys[-1], b[y0], b[y1], zs[-1], c[z0], c[z1], u[x0, y0, z0], u[x1, y0, z0], u[x0, y0, z1], u[x1, y0, z1], u[x0, y1, z0], u[x1, y1, z0], u[x0, y1, z1], u[x1, y1, z1])
        v_int1 = interpolate(xs[-1], a[x0], a[x1], ys[-1], b[y0], b[y1], zs[-1], c[z0], c[z1], v[x0, y0, z0], v[x1, y0, z0], v[x0, y0, z1], v[x1, y0, z1], v[x0, y1, z0], v[x1, y1, z0], v[x0, y1, z1], v[x1, y1, z1])
        w_int1 = interpolate(xs[-1], a[x0], a[x1], ys[-1], b[y0], b[y1], zs[-1], c[z0], c[z1], w[x0, y0, z0], w[x1, y0, z0], w[x0, y0, z1], w[x1, y0, z1], w[x0, y1, z0], w[x1, y1, z0], w[x0, y1, z1], w[x1, y1, z1])

        k1x = xs[-1] + .5*dt*u_int1
        k1y = ys[-1] + .5*dt*v_int1
        k1z = zs[-1] + .5*dt*w_int1

        x0k, x1k, y0k, y1k, z0k, z1k = find_nearest(k1x, k1y, k1z)

        u_int2 = interpolate(k1x, a[x0k], a[x1k], k1y, b[y0k], b[y1k], k1z, c[z0k], c[z1k], u[x0k, y0k, z0k], u[x1k, y0k, z0k], u[x0k, y0k, z1k], u[x1k, y0k, z1k], u[x0k, y1k, z0k], u[x1k, y1k, z0k], u[x0k, y1k, z1k], u[x1k, y1k, z1k])
        v_int2 = interpolate(k1x, a[x0k], a[x1k], k1y, b[y0k], b[y1k], k1z, c[z0k], c[z1k], v[x0k, y0k, z0k], v[x1k, y0k, z0k], v[x0k, y0k, z1k], v[x1k, y0k, z1k], v[x0k, y1k, z0k], v[x1k, y1k, z0k], v[x0k, y1k, z1k], v[x1k, y1k, z1k])
        w_int2 = interpolate(k1x, a[x0k], a[x1k], k1y, b[y0k], b[y1k], k1z, c[z0k], c[z1k], w[x0k, y0k, z0k], w[x1k, y0k, z0k], w[x0k, y0k, z1k], w[x1k, y0k, z1k], w[x0k, y1k, z0k], w[x1k, y1k, z0k], w[x0k, y1k, z1k], w[x1k, y1k, z1k])

        k2x = k1x + .5*dt*u_int2
        k2y = k1y + .5*dt*v_int2
        k2z = k1z + .5*dt*w_int2

        xs.append(k2x)
        ys.append(k2y)
        zs.append(k2z)

        iters += 1
        
    return xs, ys, zs


# Make the grid
a = np.linspace(-1, 1, 7)
b = np.linspace(-1, 1, 7)
c = np.linspace(-1, 1, 7)
A, B, C = np.meshgrid(a, b, c, indexing='ij')

# Make the direction data for the arrows
u = -B
v = A
w = .01 * np.ones(C.shape)

ax = plt.figure().add_subplot(projection='3d')
ax.quiver(A, B, C, u, v, w, length=0.1, normalize=True)
xs, ys, zs = plot_streamlines(.5, .5, -.9)
ax.plot(xs, ys, zs, color = 'black')
plt.show()
