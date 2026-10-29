"""Convergence test: how does position error shrink as the step size dt shrinks?

Circular orbit (GM = 1), integrated to t = 5 (a bit under one period of 2*pi), compared
with a high-accuracy SciPy reference. (Ending exactly at one full period makes semi-implicit
Euler look second-order because errors cancel, so a non-period end time is used.) On a log-log plot the slope is the observed order.

Run:  python convergence.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")  # remove this line to show plots in a window
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from orbit_integrators_comparison import METHODS, accel

T = 5.0
DTS = np.array([0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625])


def final_state(step, dt):
    n = int(round(T / dt))
    dt = T / n  # adjust dt slightly so the run ends exactly at t = T
    r, v = np.array([1.0, 0.0]), np.array([0.0, 1.0])
    for _ in range(n):
        r, v = step(r, v, dt)
    return r


def reference():
    rhs = lambda t, y: np.concatenate([y[2:], accel(y[:2])])
    sol = solve_ivp(rhs, [0, T], [1, 0, 0, 1], method="DOP853", rtol=1e-13, atol=1e-13)
    return sol.y[:2, -1]


def main():
    ref = reference()
    plt.figure(figsize=(6, 4.5))
    for name, f in METHODS.items():
        err = np.array([np.linalg.norm(final_state(f, dt) - ref) for dt in DTS])
        slope = np.polyfit(np.log(DTS), np.log(err), 1)[0]
        print(f"{name:20s} observed order ~ {slope:.2f}   errors: " +
              " ".join(f"{e:.1e}" for e in err))
        plt.loglog(DTS, err, "o-", label=f"{name} (slope {slope:.1f})")
    plt.xlabel("step size dt")
    plt.ylabel("position error at t = 5")
    plt.title("Convergence: error vs step size")
    plt.legend()
    plt.tight_layout()
    plt.savefig("convergence.png", dpi=120)


if __name__ == "__main__":
    main()
