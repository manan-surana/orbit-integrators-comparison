"""Closest approach between two satellites, using TLE (general perturbation) data and SGP4.

SGP4 is used only to generate the orbits. For a pair of satellites this script:
  1. propagates both over a time window and computes their separation vs time,
  2. finds the time of closest approach (TCA) and the distance at that time (DCA) numerically,
  3. predicts TCPA/DCPA from the relative position and velocity at an earlier moment using
     straight-line (constant relative velocity) motion, and compares with the numerical answer,
  4. plots separation vs time, and a zoom around the closest approach.

Straight-line prediction:   TCPA = -(r . v) / (v . v),    DCPA = |r + v * TCPA|
where r and v are the relative position and velocity (satellite 2 minus satellite 1).
It is only a good approximation for short lead times, because orbits curve.

Run:  python closest_approach.py                      (ISS and Tiangong, 24 hours)
      python closest_approach.py "Hubble" "Terra" 12   (names must match real_satellites.py; hours)
Needs: numpy, scipy, matplotlib, sgp4. Positions are in km (TEME frame).
"""
import sys
from datetime import datetime, timedelta, timezone
import numpy as np
import matplotlib
matplotlib.use("Agg")  # remove this line to show plots in a window
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from sgp4.api import Satrec

STEP_S = 10.0


def epoch_jd(sat):
    return sat.jdsatepoch + sat.jdsatepochF


def state_at(sat, jd0, t_s):
    """Position (km) and velocity (km/s) at t_s seconds after jd0 (array or scalar)."""
    t_s = np.atleast_1d(np.asarray(t_s, dtype=float))
    whole = np.floor(jd0)
    frac = jd0 - whole + t_s / 86400.0
    err, r, v = sat.sgp4_array(np.full(t_s.shape, whole), frac)
    if np.any(err != 0):
        raise RuntimeError("SGP4 error code(s): %s" % set(err[err != 0]))
    return r, v


def separation(sat1, sat2, jd0, t_s):
    r1, v1 = state_at(sat1, jd0, t_s)
    r2, v2 = state_at(sat2, jd0, t_s)
    return r2 - r1, v2 - v1


def find_closest_approach(sat1, sat2, jd0, hours):
    """Return (t_ca seconds after jd0, distance km) for the global minimum in the window."""
    t = np.arange(0, hours * 3600.0, STEP_S)
    rel_r, _ = separation(sat1, sat2, jd0, t)
    d = np.linalg.norm(rel_r, axis=1)
    i = int(np.argmin(d))
    lo, hi = t[max(i - 1, 0)], t[min(i + 1, len(t) - 1)]
    f = lambda x: float(np.linalg.norm(separation(sat1, sat2, jd0, x)[0][0]))
    res = minimize_scalar(f, bounds=(lo, hi), method="bounded", options={"xatol": 1e-4})
    return t, d, float(res.x), float(res.fun)


def predict_cpa(sat1, sat2, jd0, t_from):
    """Straight-line prediction made at time t_from: returns (tcpa seconds from t_from, dcpa km)."""
    r, v = separation(sat1, sat2, jd0, t_from)
    r, v = r[0], v[0]
    tcpa = -np.dot(r, v) / np.dot(v, v)
    dcpa = np.linalg.norm(r + v * tcpa)
    return tcpa, dcpa


def analyse(sat1, sat2, name1, name2, hours=24.0, make_plot=True, out="closest_approach.png"):
    jd0 = max(epoch_jd(sat1), epoch_jd(sat2))   # start at the later TLE epoch
    t, d, t_ca, d_ca = find_closest_approach(sat1, sat2, jd0, hours)
    ep = lambda sat: datetime(2000, 1, 1, 12, tzinfo=timezone.utc) + timedelta(days=epoch_jd(sat) - 2451545.0)
    age = lambda sat: (jd0 + t_ca / 86400.0) - epoch_jd(sat)
    print(f"{name1} vs {name2}: window {hours:g} h from the later TLE epoch")
    print(f"GP data epochs (UTC): {name1} {ep(sat1):%Y-%m-%d %H:%M}, {name2} {ep(sat2):%Y-%m-%d %H:%M}")
    print(f"GP data age at closest approach: {name1} {age(sat1):.2f} days, {name2} {age(sat2):.2f} days (accuracy degrades with age)")
    print(f"Numerical closest approach (SGP4): t = {t_ca:.1f} s ({t_ca/60:.1f} min) after start, distance = {d_ca:.2f} km")
    print("\nStraight-line TCPA/DCPA prediction made some time BEFORE the closest approach:")
    print(f"{'lead time':>10} {'predicted TCPA error (s)':>26} {'predicted DCPA (km)':>22} {'DCPA error (km)':>17}")
    rows = []
    for lead in (60, 300, 600, 1200, 1800):
        t_from = t_ca - lead
        if t_from < 0:
            continue
        tcpa, dcpa = predict_cpa(sat1, sat2, jd0, t_from)
        rows.append((lead, tcpa - lead, dcpa, dcpa - d_ca))
        print(f"{lead:>8} s {tcpa - lead:>26.2f} {dcpa:>22.2f} {dcpa - d_ca:>17.2f}")
    if make_plot:
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
        ax[0].plot(t / 3600, d, lw=0.8)
        ax[0].plot(t_ca / 3600, d_ca, "ro", label=f"closest: {d_ca:.1f} km")
        ax[0].set_xlabel("hours after start"); ax[0].set_ylabel("separation (km)")
        ax[0].set_title(f"{name1} vs {name2}: separation (SGP4)"); ax[0].legend()
        tz = np.linspace(max(t_ca - 1200, 0), t_ca + 1200, 800)
        rz, _ = separation(sat1, sat2, jd0, tz)
        ax[1].plot((tz - t_ca) / 60, np.linalg.norm(rz, axis=1), label="SGP4 (true)")
        for lead, col in ((300, "tab:orange"), (1200, "tab:green")):
            if t_ca - lead < 0:
                continue
            r0, v0 = separation(sat1, sat2, jd0, t_ca - lead)
            line = np.linalg.norm(r0[0] + v0[0] * (tz - (t_ca - lead))[:, None], axis=1)
            ax[1].plot((tz - t_ca) / 60, line, "--", color=col, label=f"straight-line, predicted {lead//60} min before")
        ax[1].axvline(0, color="gray", lw=0.5)
        ax[1].set_ylim(0, max(4 * d_ca, d_ca + 50))
        ax[1].set_xlabel("minutes from closest approach"); ax[1].set_ylabel("separation (km)")
        ax[1].set_title("Zoom: true vs straight-line prediction"); ax[1].legend(fontsize=8)
        plt.tight_layout(); plt.savefig(out, dpi=120)
    return t_ca, d_ca, rows


def main():
    from real_satellites import load_tles
    name1 = sys.argv[1] if len(sys.argv) > 1 else "ISS"
    name2 = sys.argv[2] if len(sys.argv) > 2 else "Tiangong (CSS)"
    hours = float(sys.argv[3]) if len(sys.argv) > 3 else 24.0
    tles = load_tles()
    sat1 = Satrec.twoline2rv(*tles[name1]); sat2 = Satrec.twoline2rv(*tles[name2])
    analyse(sat1, sat2, name1, name2, hours)


if __name__ == "__main__":
    main()
