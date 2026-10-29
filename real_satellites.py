"""Real satellite orbits from public TLE data (CelesTrak) propagated with SGP4.

1. Downloads current TLEs for a few satellites (cached in tle_cache.txt for 2 hours).
2. Propagates each for 3 orbits with SGP4 and plots the orbits and altitudes.
3. Compares the ISS SGP4 path with a plain two-body RK4 integration started from the
   same state, to show how much the real orbit differs from the simple model
   (SGP4 includes Earth's oblateness and drag; two-body does not).

Altitude here is |r| - 6378.137 km (spherical Earth), so it is approximate.

Run:  python real_satellites.py        Needs: numpy, matplotlib, sgp4
"""
import os
import time
import urllib.request
import numpy as np
import matplotlib
matplotlib.use("Agg")  # remove this line to show plots in a window
import matplotlib.pyplot as plt
from sgp4.api import Satrec, jday


SATS = {"ISS": 25544, "Hubble": 20580, "Tiangong (CSS)": 48274,
        "Terra": 25994, "NOAA 20": 43013}
CACHE = "tle_cache.txt"
URL = "https://celestrak.org/NORAD/elements/gp.php?CATNR={}&FORMAT=tle"
R_EARTH = 6378.137          # km
MU = 398600.4418            # km^3/s^2


def tle_checksum_ok(line):
    total = sum(int(c) if c.isdigit() else (1 if c == "-" else 0) for c in line[:68])
    return total % 10 == int(line[68])


def load_tles():
    """Return {name: (line1, line2)}; download if the cache is missing or old."""
    if not os.path.exists(CACHE) or time.time() - os.path.getmtime(CACHE) > 2 * 3600:
        lines = []
        for name, cat in SATS.items():
            with urllib.request.urlopen(URL.format(cat), timeout=30) as f:
                txt = [l.strip() for l in f.read().decode().splitlines() if l.strip()]
            if len(txt) < 2 or not txt[-2].startswith("1 "):
                print(f"Skipping {name}: no TLE returned (is it still in orbit?)")
                continue
            lines += [name, txt[-2], txt[-1]]
        with open(CACHE, "w") as f:
            f.write("\n".join(lines) + "\n")
    with open(CACHE) as f:
        L = [l.rstrip("\n") for l in f if l.strip()]
    tles = {L[i]: (L[i + 1], L[i + 2]) for i in range(0, len(L), 3)}
    for name, (l1, l2) in tles.items():
        assert tle_checksum_ok(l1) and tle_checksum_ok(l2), f"bad TLE checksum for {name}"
    return tles


def propagate(sat, n_orbits=3, step_s=30.0):
    """SGP4 positions (km, TEME) over n_orbits starting at the TLE epoch."""
    period_s = 2 * np.pi / sat.no_kozai * 60.0   # no_kozai is rad/min
    t_s = np.arange(0, n_orbits * period_s, step_s)
    jd = sat.jdsatepoch + t_s / 86400.0
    fr = sat.jdsatepochF * np.ones_like(t_s)
    err, r, v = sat.sgp4_array(jd, fr)
    if np.any(err != 0):
        raise RuntimeError("SGP4 error code(s): %s" % set(err[err != 0]))
    return t_s, r, v


def two_body_rk4(r0, v0, t_s, dt=10.0):
    """Plain two-body RK4 (km, km/s) started from the SGP4 state; sampled at t_s."""
    acc = lambda r: -MU * r / np.linalg.norm(r) ** 3

    def step(r, v, h):
        k1r, k1v = v, acc(r)
        k2r, k2v = v + 0.5 * h * k1v, acc(r + 0.5 * h * k1r)
        k3r, k3v = v + 0.5 * h * k2v, acc(r + 0.5 * h * k2r)
        k4r, k4v = v + h * k3v, acc(r + h * k3r)
        return (r + h / 6 * (k1r + 2 * k2r + 2 * k3r + k4r),
                v + h / 6 * (k1v + 2 * k2v + 2 * k3v + k4v))

    r, v, t, out = r0.copy(), v0.copy(), 0.0, []
    for target in t_s:
        while t + dt <= target + 1e-9:
            r, v = step(r, v, dt)
            t += dt
        out.append(r.copy())
    return np.array(out)


def main():
    tles = load_tles()
    fig = plt.figure(figsize=(12, 5))
    ax3 = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2)
    u, w = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
    ax3.plot_wireframe(R_EARTH * np.cos(u) * np.sin(w), R_EARTH * np.sin(u) * np.sin(w),
                       R_EARTH * np.cos(w), color="lightgray", lw=0.3)
    results = {}
    for name, (l1, l2) in tles.items():
        sat = Satrec.twoline2rv(l1, l2)
        t_s, r, v = propagate(sat)
        results[name] = (sat, t_s, r, v)
        ax3.plot(r[:, 0], r[:, 1], r[:, 2], lw=1, label=name)
        alt = np.linalg.norm(r, axis=1) - R_EARTH
        ax2.plot(t_s / 3600, alt, lw=1, label=name)
        print(f"{name:16s} mean altitude {alt.mean():7.1f} km   range {alt.min():7.1f}-{alt.max():7.1f} km")
    ax3.set_title("Real satellite orbits (SGP4, 3 orbits each)")
    ax3.set_box_aspect((1, 1, 1))
    ax3.legend(fontsize=7)
    ax2.set_xlabel("hours since TLE epoch")
    ax2.set_ylabel("approx. altitude (km)")
    ax2.set_title("Altitude over time")
    plt.tight_layout()
    plt.savefig("real_satellites.png", dpi=120)

    # Two-body vs SGP4 for the ISS
    sat, t_s, r, v = results["ISS"]
    r_tb = two_body_rk4(r[0], v[0], t_s)
    diff = np.linalg.norm(r_tb - r, axis=1)
    plt.figure(figsize=(6, 4.5))
    plt.plot(t_s / 3600, diff)
    plt.xlabel("hours since TLE epoch")
    plt.ylabel("distance between two-body RK4 and SGP4 (km)")
    plt.title("ISS: simple two-body model vs SGP4")
    plt.tight_layout()
    plt.savefig("iss_twobody_vs_sgp4.png", dpi=120)
    print(f"ISS two-body vs SGP4 after 3 orbits: {diff[-1]:.1f} km")


if __name__ == "__main__":
    main()
