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
    next_x = first_x
    next_y = first_y
    next_z = first_z
    iters = 0
    while next_x < 1 and next_x > -1 and next_y < 1 and next_y > -1 and next_z < 1 and next_z > -1 and iters < 100000:
        x0, x1, y0, y1, z0, z1 = find_nearest(xs[-1], ys[-1], zs[-1])
        dt = .01
        u_int = interpolate(xs[-1], a[x0], a[x1], ys[-1], b[y0], b[y1], zs[-1], c[z0], c[z1], u[x0, y0, z0], u[x1, y0, z0], u[x0, y0, z1], u[x1, y0, z1], u[x0, y1, z0], u[x1, y1, z0], u[x0, y1, z1], u[x1, y1, z1])
        v_int = interpolate(xs[-1], a[x0], a[x1], ys[-1], b[y0], b[y1], zs[-1], c[z0], c[z1], v[x0, y0, z0], v[x1, y0, z0], v[x0, y0, z1], v[x1, y0, z1], v[x0, y1, z0], v[x1, y1, z0], v[x0, y1, z1], v[x1, y1, z1])
        w_int = interpolate(xs[-1], a[x0], a[x1], ys[-1], b[y0], b[y1], zs[-1], c[z0], c[z1], w[x0, y0, z0], w[x1, y0, z0], w[x0, y0, z1], w[x1, y0, z1], w[x0, y1, z0], w[x1, y1, z0], w[x0, y1, z1], w[x1, y1, z1])
        next_x = xs[-1] + dt*u_int
        next_y = ys[-1] + dt*v_int
        next_z = zs[-1] + dt*w_int
        xs.append(next_x)
        ys.append(next_y)
        zs.append(next_z)
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

#higher order (RK2)
