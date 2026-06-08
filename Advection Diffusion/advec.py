import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ============================================================
# DOMAIN PARAMETERS
# ============================================================

Nx, Ny = 200, 200

Lx, Ly = 2.0, 2.0

dx = Lx / Nx
dy = Ly / Ny

# ============================================================
# TIME PARAMETERS
# ============================================================

dt = 0.0001
Nt = 4000

# physics steps per rendered frame
steps_per_frame = 50

# ============================================================
# PHYSICAL PARAMETERS
# ============================================================

# very weak physical diffusion
D = 1e-6

# vortex angular strength
strength = 8.0

# ============================================================
# CREATE GRID
# ============================================================

x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)

X, Y = np.meshgrid(x, y)

# ============================================================
# VORTEX VELOCITY FIELD
# ============================================================

x0 = Lx / 2
y0 = Ly / 2

# solid-body rotation

u = -strength * (Y - y0)
v =  strength * (X - x0)

# ============================================================
# INITIAL DYE BLOB
# ============================================================

c = np.exp(
    -(
        (X - 0.7)**2 +
        (Y - 1.0)**2
    ) / 0.005
)

# ============================================================
# PLOT SETUP
# ============================================================

fig, ax = plt.subplots(figsize=(7,7))

im = ax.imshow(
    c,
    origin='lower',
    extent=[0, Lx, 0, Ly],
    cmap='plasma',
    vmin=0,
    vmax=1
)

plt.colorbar(im, ax=ax, label='Concentration')

# streamlines of velocity field

ax.streamplot(
    x,
    y,
    u,
    v,
    color='white',
    density=1.5,
    linewidth=0.7
)

ax.set_title("2D Advection-Diffusion with Upwind Scheme")
ax.set_xlabel("x")
ax.set_ylabel("y")

# ============================================================
# UPDATE FUNCTION
# ============================================================

def update(frame):

    global c

    # --------------------------------------------------------
    # ADVANCE MANY PHYSICS STEPS PER FRAME
    # --------------------------------------------------------

    for _ in range(steps_per_frame):

        c_new = c.copy()

        # ----------------------------------------------------
        # INTERIOR VELOCITY FIELD
        # ----------------------------------------------------

        u_inner = u[1:-1, 1:-1]
        v_inner = v[1:-1, 1:-1]

        # ====================================================
        # UPWIND FIRST DERIVATIVES
        # ====================================================

        # ----------------------------------------------------
        # x derivative
        # ----------------------------------------------------

        c_x = np.where(

            u_inner > 0,

            # backward difference
            (
                c[1:-1, 1:-1]
                - c[1:-1, :-2]
            ) / dx,

            # forward difference
            (
                c[1:-1, 2:]
                - c[1:-1, 1:-1]
            ) / dx
        )

        # ----------------------------------------------------
        # y derivative
        # ----------------------------------------------------

        c_y = np.where(

            v_inner > 0,

            # backward difference
            (
                c[1:-1, 1:-1]
                - c[:-2, 1:-1]
            ) / dy,

            # forward difference
            (
                c[2:, 1:-1]
                - c[1:-1, 1:-1]
            ) / dy
        )

        # ====================================================
        # SECOND DERIVATIVES (DIFFUSION)
        # ====================================================

        c_xx = (
            c[1:-1, 2:]
            - 2*c[1:-1, 1:-1]
            + c[1:-1, :-2]
        ) / dx**2

        c_yy = (
            c[2:, 1:-1]
            - 2*c[1:-1, 1:-1]
            + c[:-2, 1:-1]
        ) / dy**2

        # ====================================================
        # PHYSICS TERMS
        # ====================================================

        advection = (
            u_inner * c_x
            + v_inner * c_y
        )

        diffusion = D * (c_xx + c_yy)

        # ====================================================
        # EXPLICIT EULER UPDATE
        # ====================================================

        c_new[1:-1, 1:-1] = (

            c[1:-1, 1:-1]

            + dt * (

                -advection
                + diffusion

            )
        )

        # ====================================================
        # BOUNDARY CONDITIONS
        # ====================================================

        # no-flux top/bottom

        c_new[0, :] = c_new[1, :]
        c_new[-1, :] = c_new[-2, :]

        # no-flux left/right

        c_new[:, 0] = c_new[:, 1]
        c_new[:, -1] = c_new[:, -2]

        # update concentration field

        c = c_new

    # ========================================================
    # UPDATE VISUALISATION
    # ========================================================

    im.set_array(c)

    # dynamic colour scaling

    im.set_clim(
        vmin=0,
        vmax=np.max(c)
    )

    return [im]

# ============================================================
# ANIMATION
# ============================================================

anim = FuncAnimation(
    fig,
    update,
    frames=Nt,
    interval=5,
    blit=False
)

plt.tight_layout()

plt.show()