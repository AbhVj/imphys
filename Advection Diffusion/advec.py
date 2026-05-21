import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ============================================================
# PARAMETERS
# ============================================================

Nx, Ny = 200, 100
Lx, Ly = 2.0, 1.0

dx = Lx / Nx
dy = Ly / Ny

dt = 0.0002
Nt = 3000

D = 0.01

u = 1.0
v = 0.0

# ============================================================
# GRID
# ============================================================

x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)

X, Y = np.meshgrid(x, y)

# ============================================================
# INITIAL DYE BLOB
# ============================================================

# Blob near left boundary

c = np.exp(
    -(
        (X - 0.2)**2 +
        (Y - 0.5)**2
    ) / 0.005
)

# ============================================================
# PLOT SETUP
# ============================================================

fig, ax = plt.subplots(figsize=(10,4))

im = ax.imshow(
    c,
    origin='lower',
    extent=[0, Lx, 0, Ly],
    cmap='viridis',
    vmin=0,
    vmax=0.2,
    aspect='auto'
)

plt.colorbar(im, ax=ax, label='Concentration')

ax.set_title("Advection-Diffusion of Dye Blob")
ax.set_xlabel("x")
ax.set_ylabel("y")

# ============================================================
# UPDATE FUNCTION
# ============================================================

def update(frame):

    global c

    c_new = c.copy()

    # --------------------------------------------------------
    # DERIVATIVES
    # --------------------------------------------------------

    c_x = (
        c[1:-1, 2:] - c[1:-1, :-2]
    ) / (2 * dx)

    c_y = (
        c[2:, 1:-1] - c[:-2, 1:-1]
    ) / (2 * dy)

    c_xx = (
        c[1:-1, 2:]
        - 2 * c[1:-1, 1:-1]
        + c[1:-1, :-2]
    ) / dx**2

    c_yy = (
        c[2:, 1:-1]
        - 2 * c[1:-1, 1:-1]
        + c[:-2, 1:-1]
    ) / dy**2

    # --------------------------------------------------------
    # PHYSICS TERMS
    # --------------------------------------------------------

    advection = u * c_x + v * c_y

    diffusion = D * (c_xx + c_yy)

    # --------------------------------------------------------
    # UPDATE INTERIOR
    # --------------------------------------------------------

    c_new[1:-1, 1:-1] = (
        c[1:-1, 1:-1]
        + dt * (-advection + diffusion)
    )

    # --------------------------------------------------------
    # BOUNDARY CONDITIONS
    # --------------------------------------------------------

    # Top/bottom no-flux

    c_new[0, :] = c_new[1, :]
    c_new[-1, :] = c_new[-2, :]

    # Left boundary:
    # allow clean inflow/outflow behaviour

    c_new[:, 0] = c_new[:, 1]

    # Right boundary:
    # zero-gradient outflow

    c_new[:, -1] = c_new[:, -2]

    c = c_new

    im.set_array(c)

    return [im]

# ============================================================
# ANIMATION
# ============================================================

anim = FuncAnimation(
    fig,
    update,
    frames=Nt,
    interval=5,
    blit=True
)

plt.show()
