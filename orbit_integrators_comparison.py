"""Compare numerical integrators for a satellite orbit (units: GM = 1).

A circular orbit at radius 1 has speed 1 and period 2*pi.
Methods: explicit Euler, semi-implicit Euler, RK4, leapfrog.
Checks: relative energy error, angular-momentum error, and the position
difference from a high-accuracy SciPy reference (DOP853).

Run:  python orbit_integrators_comparison.py
Needs: numpy, scipy, matplotlib
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")  # remove this line to show plots in a window
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

GM = 1.0


def accel(r):
    return -GM * r / np.linalg.norm(r) ** 3


def energy(r, v):
    return 0.5 * np.dot(v, v) - GM / np.linalg.norm(r)


def ang_mom(r, v):
    return r[0] * v[1] - r[1] * v[0]


def explicit_euler(r, v, dt):
    return r + v * dt, v + accel(r) * dt


def semi_implicit_euler(r, v, dt):
    v_new = v + accel(r) * dt
    return r + v_new * dt, v_new


def rk4(r, v, dt):
    k1r, k1v = v, accel(r)
    k2r, k2v = v + 0.5 * dt * k1v, accel(r + 0.5 * dt * k1r)
    k3r, k3v = v + 0.5 * dt * k2v, accel(r + 0.5 * dt * k2r)
    k4r, k4v = v + dt * k3v, accel(r + dt * k3r)
    r_new = r + dt / 6 * (k1r + 2 * k2r + 2 * k3r + k4r)
    v_new = v + dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
    return r_new, v_new


def leapfrog(r, v, dt):
    v_half = v + 0.5 * dt * accel(r)
    r_new = r + dt * v_half
    return r_new, v_half + 0.5 * dt * accel(r_new)


METHODS = {
    "Explicit Euler": explicit_euler,
    "Semi-implicit Euler": semi_implicit_euler,
    "RK4": rk4,
    "Leapfrog": leapfrog,
}


def run(step, dt, n_orbits=10, v0=1.0):
    """Return positions, relative energy error, angular-momentum error."""
    n = int(n_orbits * 2 * np.pi / dt)
    r, v = np.array([1.0, 0.0]), np.array([0.0, v0])
    E0, L0 = energy(r, v), ang_mom(r, v)
    pos, e_err, l_err = [], [], []
    for _ in range(n):
        r, v = step(r, v, dt)
        pos.append(r.copy())
        e_err.append(abs((energy(r, v) - E0) / E0))
        l_err.append(abs(ang_mom(r, v) - L0))
    return np.array(pos), np.array(e_err), np.array(l_err)


def reference(dt, n_orbits=10, v0=1.0):
    n = int(n_orbits * 2 * np.pi / dt)
    times = dt * np.arange(1, n + 1)
    rhs = lambda t, y: np.concatenate([y[2:], accel(y[:2])])
    sol = solve_ivp(rhs, [0, times[-1]], [1, 0, 0, v0], method="DOP853",
                    rtol=1e-12, atol=1e-12, t_eval=times)
    return sol.y[:2].T


def main():
    dt = 0.01
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for ax, v0, title in [(axes[0], 1.0, "Circular orbit (v = 1.0)"),
                          (axes[1], 0.8, "Elliptical orbit (v = 0.8)")]:
        ref = reference(dt, v0=v0)
        print(f"\n{title}, dt={dt}, 10 orbits")
        for name, f in METHODS.items():
            pos, e_err, l_err = run(f, dt, v0=v0)
            t = dt * np.arange(1, len(e_err) + 1)
            ax.semilogy(t, np.maximum(e_err, 1e-18), label=name, lw=1)
            d = np.linalg.norm(pos - ref, axis=1).max()
            print(f"{name:20s} final E err {e_err[-1]:.2e}  max E err {e_err.max():.2e}"
                  f"  final L err {l_err[-1]:.2e}  max pos diff vs SciPy {d:.2e}")
        ax.set_title(title)
        ax.set_xlabel("time (orbit period = 2*pi)")
    axes[0].set_ylabel("relative energy error")
    axes[0].legend()
    plt.tight_layout()
    plt.savefig("energy_error_comparison.png", dpi=120)


if __name__ == "__main__":
    main()
