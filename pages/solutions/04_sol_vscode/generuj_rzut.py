"""Generate figures for the oblique projectile motion report.

Run from pages/solutions/04_sol_vscode with:
    python generuj_rzut.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

OUT = Path("files/rzut_ukosny")
OUT.mkdir(parents=True, exist_ok=True)

V0 = 26.0
ALPHA_DEG = 97.0
ALPHA = np.deg2rad(ALPHA_DEG)
G = 9.81
M = 0.20
C_LINEAR = 0.12
C_QUADRATIC = 0.006
GAMMA = C_LINEAR / M
BETA = C_QUADRATIC / M
VX0 = V0 * np.cos(ALPHA)
VY0 = V0 * np.sin(ALPHA)

BLUE, ORANGE, GREEN = "#2463A6", "#D47722", "#168578"
PURPLE, INK, GRAY = "#8455A4", "#233044", "#64748B"
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.labelcolor": INK,
    "text.color": INK,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.2,
    "lines.linewidth": 2.2,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=180)
    plt.close(fig)


def no_drag(t):
    x = VX0 * t
    y = VY0 * t - 0.5 * G * t**2
    vx = np.full_like(t, VX0)
    vy = VY0 - G * t
    ax = np.zeros_like(t)
    ay = np.full_like(t, -G)
    return x, y, vx, vy, ax, ay


def linear_drag(t):
    decay = np.exp(-GAMMA * t)
    vx = VX0 * decay
    vy = (VY0 + G / GAMMA) * decay - G / GAMMA
    x = VX0 * (1 - decay) / GAMMA
    y = (VY0 + G / GAMMA) * (1 - decay) / GAMMA - G * t / GAMMA
    ax = -GAMMA * vx
    ay = -G - GAMMA * vy
    return x, y, vx, vy, ax, ay


def quadratic_rhs(_t, state):
    x, y, vx, vy = state
    speed = np.hypot(vx, vy)
    return [vx, vy, -BETA * speed * vx, -G - BETA * speed * vy]


def ground_event(_t, state):
    return state[1]


ground_event.terminal = True
ground_event.direction = -1


def quadratic_solution():
    result = solve_ivp(
        quadratic_rhs,
        (0.0, 60.0),
        [0.0, 0.0, VX0, VY0],
        events=ground_event,
        dense_output=True,
        rtol=1e-10,
        atol=1e-12,
        max_step=0.04,
    )
    if not result.success or len(result.t_events[0]) != 1:
        raise RuntimeError("The quadratic-drag trajectory did not reach the ground.")
    return result, float(result.t_events[0][0])


T0 = 2 * VY0 / G
TL = brentq(lambda t: linear_drag(np.array([t]))[1][0], 1e-9, 60.0)
Q_SOL, TQ = quadratic_solution()
TIMES = {
    "Bez oporu": np.linspace(0, T0, 700),
    "Opór liniowy": np.linspace(0, TL, 700),
    "Opór kwadratowy": np.linspace(0, TQ, 700),
}


def states():
    outputs = {}
    t = TIMES["Bez oporu"]
    outputs["Bez oporu"] = (t, no_drag(t))
    t = TIMES["Opór liniowy"]
    outputs["Opór liniowy"] = (t, linear_drag(t))
    t = TIMES["Opór kwadratowy"]
    values = Q_SOL.sol(t)
    speed = np.hypot(values[2], values[3])
    ax = -BETA * speed * values[2]
    ay = -G - BETA * speed * values[3]
    outputs["Opór kwadratowy"] = (t, (*values, ax, ay))
    return outputs


DATA = states()
COLORS = {
    "Bez oporu": BLUE,
    "Opór liniowy": ORANGE,
    "Opór kwadratowy": GREEN,
}

# Figure 1: a schematic defines the launch angle and coordinate axes.
fig, ax = plt.subplots(figsize=(8.5, 4.7), layout="constrained")
schematic_t = np.linspace(0, 0.25, 200)
sx, sy, *_ = no_drag(schematic_t)
ax.plot(sx, sy, color=BLUE, label="Tor bez oporu (fragment)")
ax.scatter([0], [0], color=INK, zorder=4)
ax.annotate("", xy=(2.0, 0), xytext=(0, 0),
            arrowprops={"arrowstyle": "->", "color": GRAY, "lw": 1.5})
ax.annotate("", xy=(0, 5.2), xytext=(0, 0),
            arrowprops={"arrowstyle": "->", "color": GRAY, "lw": 1.5})
ax.annotate("", xy=(VX0 * 0.14, VY0 * 0.14), xytext=(0, 0),
            arrowprops={"arrowstyle": "->", "color": ORANGE, "lw": 2.5})
ax.text(VX0 * 0.14 - 0.3, VY0 * 0.14 + 0.4, r"$\vec{v}_0$", color=ORANGE)
ax.text(1.0, 0.35, r"$x$", color=GRAY)
ax.text(0.25, 4.6, r"$y$", color=GRAY)
ax.text(-0.75, 0.25, r"$\alpha=97^\circ$", color=INK)
ax.set(xlabel="Położenie x [m]", ylabel="Wysokość y [m]",
       xlim=(-2.2, 1.5), ylim=(-0.8, 8.0))
ax.set_title("Rzut ukośny: układ współrzędnych i kierunek prędkości początkowej", loc="left")
ax.legend(loc="upper right", fontsize=9)
save(fig, "schemat")

# Figure 2: trajectories and impact points for all models.
fig, ax = plt.subplots(figsize=(9.5, 5.4), layout="constrained")
for name, (t, (x, y, *_)) in DATA.items():
    ax.plot(x, y, color=COLORS[name], label=name)
    ax.scatter([x[-1]], [0], color=COLORS[name], s=30, zorder=4)
ax.axhline(0, color=INK, lw=1)
ax.set(xlabel="Położenie poziome x [m]", ylabel="Wysokość y [m]")
ax.set_title("Tory ruchu dla trzech modeli oporu", loc="left")
ax.legend()
save(fig, "tory_ruchu")

# Figure 3: position and velocity components as functions of time.
fig, axes = plt.subplots(2, 1, figsize=(9.5, 7.0), sharex=True,
                         layout="constrained")
for name, (t, (x, y, vx, vy, *_)) in DATA.items():
    axes[0].plot(t, x, color=COLORS[name], label=f"x, {name}")
    axes[0].plot(t, y, color=COLORS[name], ls="--", label=f"y, {name}")
    axes[1].plot(t, vx, color=COLORS[name], label=f"vₓ, {name}")
    axes[1].plot(t, vy, color=COLORS[name], ls="--", label=f"vᵧ, {name}")
axes[0].set_ylabel("Położenie [m]")
axes[1].set(xlabel="Czas t [s]", ylabel="Składowa prędkości [m/s]")
axes[0].set_title("Położenie i składowe prędkości w czasie", loc="left")
for axis in axes:
    axis.axhline(0, color=INK, lw=0.7, alpha=0.5)
    axis.legend(ncol=2, fontsize=8)
save(fig, "polozenie_predkosc")

# Figure 4: speed and acceleration magnitude.
fig, axes = plt.subplots(1, 2, figsize=(10, 4.3), layout="constrained")
for name, (t, (x, y, vx, vy, ax_values, ay_values)) in DATA.items():
    axes[0].plot(t, np.hypot(vx, vy), color=COLORS[name], label=name)
    axes[1].plot(t, np.hypot(ax_values, ay_values), color=COLORS[name], label=name)
axes[0].set(xlabel="Czas t [s]", ylabel="Szybkość |v| [m/s]", title="Szybkość")
axes[1].set(xlabel="Czas t [s]", ylabel="Wartość przyspieszenia |a| [m/s²]",
            title="Wartość przyspieszenia")
axes[0].legend(fontsize=8)
axes[1].legend(fontsize=8)
save(fig, "szybkosc_przyspieszenie")

# Figure 5: mechanical energy and dissipated energy for the drag models.
fig, ax = plt.subplots(figsize=(9.5, 4.8), layout="constrained")
energy0 = 0.5 * M * V0**2
for name in ("Opór liniowy", "Opór kwadratowy"):
    t, (x, y, vx, vy, *_rest) = DATA[name]
    kinetic = 0.5 * M * (vx**2 + vy**2)
    potential = M * G * y
    ax.plot(t, kinetic + potential, color=COLORS[name], label=f"E mechaniczna: {name}")
    ax.plot(t, energy0 - kinetic - potential, color=COLORS[name], ls="--",
            label=f"Energia rozproszona: {name}")
ax.axhline(energy0, color=INK, ls=":", label="E₀")
ax.set(xlabel="Czas t [s]", ylabel="Energia [J]", ylim(bottom=0))
ax.set_title("Energia mechaniczna maleje wskutek oporu", loc="left")
ax.legend(ncol=2, fontsize=8)
save(fig, "energia")

# Figure 6: sensitivity of signed horizontal displacement to launch angle.
angles = np.linspace(5, 175, 500)
angles_rad = np.deg2rad(angles)
range_no_drag = V0**2 * np.sin(2 * angles_rad) / G
fig, ax = plt.subplots(figsize=(9.5, 4.4), layout="constrained")
ax.plot(angles, range_no_drag, color=PURPLE)
ax.axvline(ALPHA_DEG, color=ORANGE, ls="--", label=r"$\alpha=97^\circ$")
ax.axhline(0, color=INK, lw=0.8)
ax.scatter([ALPHA_DEG], [V0**2 * np.sin(2 * ALPHA) / G], color=ORANGE, zorder=4)
ax.set(xlabel="Kąt początkowy α [°]", ylabel="Przemieszczenie poziome do lądowania [m]")
ax.set_title("Znak zasięgu zależy od strony pionu, w którą rzucono ciało", loc="left")
ax.legend()
save(fig, "zasieg_kat")

# Numerical checks: ground level, initial speed, and monotonic energy loss.
for name, (t, (x, y, vx, vy, ax_values, ay_values)) in DATA.items():
    assert abs(y[0]) < 1e-12
    assert abs(y[-1]) < 1e-7
    assert abs(np.hypot(vx[0], vy[0]) - V0) < 1e-10
for name in ("Opór liniowy", "Opór kwadratowy"):
    t, (x, y, vx, vy, *_rest) = DATA[name]
    mechanical = 0.5 * M * (vx**2 + vy**2) + M * G * y
    assert np.max(np.diff(mechanical)) < 1e-7

print(f"v0 = {V0:.2f} m/s; alpha = {ALPHA_DEG:.1f} deg; g = {G:.2f} m/s^2")
print(f"Bez oporu: T = {T0:.4f} s; x_ladowania = {VX0*T0:.4f} m")
print(f"Opór liniowy: c = {C_LINEAR:.3f} kg/s; T = {TL:.4f} s; x_ladowania = {linear_drag(np.array([TL]))[0][0]:.4f} m")
print(f"Opór kwadratowy: c_d = {C_QUADRATIC:.3f} kg/m; T = {TQ:.4f} s; x_ladowania = {Q_SOL.y[0, -1]:.4f} m")
print(f"Zapisano 6 wykresów w {OUT.resolve()}")