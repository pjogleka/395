import matplotlib.pyplot as plt
import numpy as np

def find_nearest(x, y, z, a, b, c):
    
    err = 0
    
    x0 = np.floor(((x-np.min(a))/(np.max(a)-np.min(a))*(len(a)-1))).astype(int)
    if x0 >= len(a) - 1:
        err = 1
    x1 = x0 + 1

    y0 = np.floor(((y-np.min(b))/(np.max(b)-np.min(b))*(len(b)-1))).astype(int)
    if y0 >= len(b) - 1:
        err = 2
    y1 = y0 + 1

    z0 = np.floor(((z-np.min(c))/(np.max(c)-np.min(c))*(len(c)-1))).astype(int)
    if z0 >= len(c) - 1:
        err = 3
    z1 = z0 + 1

    return x0, x1, y0, y1, z0, z1, err


def interpolate(x, x0, x1, y, y0, y1, z, z0, z1, c000, c100, c001, c101, c010, c110, c011, c111):

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


def plot_streamlines(first_x, first_y, first_z, a, b, c, u, v, w, print):

    xs = [first_x]
    ys = [first_y]
    zs = [first_z]
    x = first_x
    y = first_y
    z = first_z
    iters = 0

    while iters < 100000:

        dt = .01

        x = xs[-1]
        y = ys[-1]
        z = zs[-1]
        x0, x1, y0, y1, z0, z1, err = find_nearest(x, y, z, a, b, c)
        if err:
            break
        unnn = interpolate(x, a[x0], a[x1], y, b[y0], b[y1], z, c[z0], c[z1], u[x0, y0, z0], u[x1, y0, z0], u[x0, y0, z1], u[x1, y0, z1], u[x0, y1, z0], u[x1, y1, z0], u[x0, y1, z1], u[x1, y1, z1])
        xn13 = x + (1/3)*dt*unnn

        x = xn13
        y = ys[-1]
        z = zs[-1]
        x0, x1, y0, y1, z0, z1, err = find_nearest(x, y, z, a, b, c)
        if err:
            break
        v1nn = interpolate(x, a[x0], a[x1], y, b[y0], b[y1], z, c[z0], c[z1], v[x0, y0, z0], v[x1, y0, z0], v[x0, y0, z1], v[x1, y0, z1], v[x0, y1, z0], v[x1, y1, z0], v[x0, y1, z1], v[x1, y1, z1])
        yn13 = y + (1/3)*dt*v1nn

        x = xn13
        y = yn13
        z = zs[-1]
        x0, x1, y0, y1, z0, z1, err = find_nearest(x, y, z, a, b, c)
        if err:
            break
        w11n = interpolate(x, a[x0], a[x1], y, b[y0], b[y1], z, c[z0], c[z1], w[x0, y0, z0], w[x1, y0, z0], w[x0, y0, z1], w[x1, y0, z1], w[x0, y1, z0], w[x1, y1, z0], w[x0, y1, z1], w[x1, y1, z1])
        zn13 = z + (1/3)*dt*w11n

        x = xn13
        y = yn13
        z = zn13
        x0, x1, y0, y1, z0, z1, err = find_nearest(x, y, z, a, b, c)
        if err:
            break
        v111 = interpolate(x, a[x0], a[x1], y, b[y0], b[y1], z, c[z0], c[z1], v[x0, y0, z0], v[x1, y0, z0], v[x0, y0, z1], v[x1, y0, z1], v[x0, y1, z0], v[x1, y1, z0], v[x0, y1, z1], v[x1, y1, z1])
        yn23 = y + (1/3)*dt*v111

        x = xn13
        y = yn23
        z = zn13
        x0, x1, y0, y1, z0, z1, err = find_nearest(x, y, z, a, b, c)
        if err:
            break
        u121 = interpolate(x, a[x0], a[x1], y, b[y0], b[y1], z, c[z0], c[z1], u[x0, y0, z0], u[x1, y0, z0], u[x0, y0, z1], u[x1, y0, z1], u[x0, y1, z0], u[x1, y1, z0], u[x0, y1, z1], u[x1, y1, z1])
        xn23 = x + (1/3)*dt*u121

        x = xn23
        y = yn23
        z = zn13
        x0, x1, y0, y1, z0, z1, err = find_nearest(x, y, z, a, b, c)
        if err:
            break
        w221 = interpolate(x, a[x0], a[x1], y, b[y0], b[y1], z, c[z0], c[z1], w[x0, y0, z0], w[x1, y0, z0], w[x0, y0, z1], w[x1, y0, z1], w[x0, y1, z0], w[x1, y1, z0], w[x0, y1, z1], w[x1, y1, z1])
        zn23 = z + (1/3)*dt*w221

        x = xn23
        y = yn23
        z = zn23
        x0, x1, y0, y1, z0, z1, err = find_nearest(x, y, z, a, b, c)
        if err:
            break
        w222 = interpolate(x, a[x0], a[x1], y, b[y0], b[y1], z, c[z0], c[z1], w[x0, y0, z0], w[x1, y0, z0], w[x0, y0, z1], w[x1, y0, z1], w[x0, y1, z0], w[x1, y1, z0], w[x0, y1, z1], w[x1, y1, z1])
        zn_1 = z + (1/3)*dt*w222

        x = xn23
        y = yn23
        z = zn_1
        x0, x1, y0, y1, z0, z1, err = find_nearest(x, y, z, a, b, c)
        if err:
            break
        v22_ = interpolate(x, a[x0], a[x1], y, b[y0], b[y1], z, c[z0], c[z1], v[x0, y0, z0], v[x1, y0, z0], v[x0, y0, z1], v[x1, y0, z1], v[x0, y1, z0], v[x1, y1, z0], v[x0, y1, z1], v[x1, y1, z1])
        yn_1 = y + (1/3)*dt*v22_

        x = xn23
        y = yn_1
        z = zn_1
        x0, x1, y0, y1, z0, z1, err = find_nearest(x, y, z, a, b, c)
        if err:
            break
        u2__ = interpolate(x, a[x0], a[x1], y, b[y0], b[y1], z, c[z0], c[z1], u[x0, y0, z0], u[x1, y0, z0], u[x0, y0, z1], u[x1, y0, z1], u[x0, y1, z0], u[x1, y1, z0], u[x0, y1, z1], u[x1, y1, z1])
        xn_1 = x + (1/3)*dt*u2__

        xs.append(xn_1)
        ys.append(yn_1)
        zs.append(zn_1)

        iters += 1

    if err and print:
        print("Domain reached along axis %1i" % (err - 1))
    plt.plot(xs, ys, zs, color = 'black')

"""
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
plot_streamlines(.5, .5, -.9, a, b, c, u, v, w, False)
plt.show()
"""

