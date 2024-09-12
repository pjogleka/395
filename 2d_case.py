import matplotlib.pyplot as plt
import numpy as np

def interpolate(x, y, x1, y1, x2, y2, q11, q12, q21, q22):
    R1 = ((x2-x)/(x2-x1))*q11 + ((x-x1)/(x2-x1))*q21
    R2 = ((x2-x)/(x2-x1))*q12 + ((x-x1)/(x2-x1))*q22
    S = ((y2-y)/(y2-y1))*R1 + ((y-y1)/(y2-y1))*R2
    return S

def find_nearest(x, y):
    #returns indices of surrounding x and y values

    x_low = np.floor(((x-np.min(a))/(np.max(a)-np.min(a))*(len(a)-1))).astype(int)
    x_high = x_low + 1

    y_low = np.floor(((y-np.min(b))/(np.max(b)-np.min(b))*(len(b)-1))).astype(int)
    y_high = y_low + 1

    """
    x_less = np.where(x < a)
    x_low = np.max(x_less)

    x_more = np.where(x > a)
    x_high = np.min(x_more)

    y_less = np.where(y < b)
    y_low = np.max(y_less)

    y_more = np.where(y > b)
    y_high = np.min(y_more)
    """

    return x_low, x_high, y_low, y_high

def plot_streamlines(first_x, first_y):
    xs = [first_x]
    ys = [first_y]
    next_x = first_x
    next_y = first_y
    iters = 0
    while next_x < 4 and next_x > -4 and next_y < 4 and next_y > -4 and iters < 100000:
    #for i in range(4):
        x_low, x_high, y_low, y_high = find_nearest(xs[-1], ys[-1])
        dt = .1
        u_int = interpolate(xs[-1], ys[-1], a[x_low], b[y_low], a[x_high], b[y_high], U[x_low, y_low], U[x_low, y_high], U[x_high, y_low], U[x_high, y_high])
        v_int = interpolate(xs[-1], ys[-1], a[x_low], b[y_low], a[x_high], b[y_high], V[x_low, y_low], V[x_low, y_high], V[x_high, y_low], V[x_high, y_high])
        next_x = xs[-1] + dt*u_int
        next_y = ys[-1] + dt*v_int
        xs.append(next_x)
        ys.append(next_y)
        #print(f"xl= {x_low}, xh={x_high}, yl={y_low}, yh={y_high}, u={u_int},v= {v_int}, x={xs[-2]}, y={ys[-2]}")
        #print(f"v11={V[x_low, y_low]}, v12={V[x_low, y_high]}, v21={V[x_high, y_low]}, v22={V[x_high, y_high]}")
        iters += 1
    plt.plot(xs, ys, 'pink')

# make data
a = np.linspace(-4, 4, 9)
b = np.linspace(-4, 4, 9)
A, B = np.meshgrid(a, b, indexing = 'ij')
U = A + B
V = B - A

# plot
fig, ax = plt.subplots()
ax.quiver(A, B, U, V, angles='xy')
ax.set(xlim=(-5, 5), ylim=(-5, 5))

plot_streamlines(.3,.5)
plt.show()

# make data
a = np.linspace(-4, 4, 9)
b = np.linspace(-4, 4, 9)
A, B = np.meshgrid(a, b, indexing= 'ij')
U = -B
V = A/2

# plot
fig, ax = plt.subplots()
ax.quiver(A, B, U, V, angles='xy')
ax.set(xlim=(-5, 5), ylim=(-5, 5))

plot_streamlines(.1, .5)
plt.show()
