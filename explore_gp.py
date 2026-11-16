"""Look inside general perturbation (GP) data: decode every TLE in tle_cache.txt.

A TLE holds MEAN elements that only make sense together with SGP4 (the general perturbation
theory they were fitted with). Do not treat them as ordinary Keplerian orbital elements.

Run:  python explore_gp.py        (run real_satellites.py first so tle_cache.txt exists)
Needs: numpy, sgp4
"""
from datetime import datetime, timedelta, timezone
import numpy as np
from sgp4.api import Satrec

MU = 398600.4418      # km^3/s^2
R_EARTH = 6378.137    # km


def jd_to_utc(jd):
    return datetime(2000, 1, 1, 12, tzinfo=timezone.utc) + timedelta(days=jd - 2451545.0)


def describe(name, sat, now=None):
    now = now or datetime.now(timezone.utc)
    epoch = jd_to_utc(sat.jdsatepoch + sat.jdsatepochF)
    n = sat.no_kozai / 60.0                       # rad/s (mean motion)
    a = (MU / n**2) ** (1 / 3)                    # semi-major axis from mean motion, km
    print(f"{name}  (catalogue number {sat.satnum})")
    print(f"  epoch (UTC)            : {epoch:%Y-%m-%d %H:%M:%S}   age now: {(now - epoch).total_seconds()/86400:.2f} days")
    print(f"  inclination            : {np.degrees(sat.inclo):.4f} deg")
    print(f"  RAAN                   : {np.degrees(sat.nodeo):.4f} deg")
    print(f"  eccentricity           : {sat.ecco:.7f}")
    print(f"  argument of perigee    : {np.degrees(sat.argpo):.4f} deg")
    print(f"  mean anomaly           : {np.degrees(sat.mo):.4f} deg")
    print(f"  mean motion            : {sat.no_kozai * 1440 / (2*np.pi):.8f} rev/day  (period {2*np.pi/sat.no_kozai:.2f} min)")
    print(f"  BSTAR (drag term)      : {sat.bstar:.3e}")
    print(f"  approx. altitude       : {a - R_EARTH:.1f} km (from mean motion, spherical Earth)")
    print()


if __name__ == "__main__":
    from real_satellites import load_tles
    for name, (l1, l2) in load_tles().items():
        describe(name, Satrec.twoline2rv(l1, l2))
