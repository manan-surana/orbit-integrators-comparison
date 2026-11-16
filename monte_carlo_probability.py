"""Monte Carlo collision probability at the time of closest approach (TCA).

Follows the same idea as IIIT Delhi's OrCo platform (Chowdhury et al., IPSC 2025):
SGP4 gives the nominal (best-estimate) trajectory from TLE data; since a TLE carries no
stated uncertainty, we assume a position-error standard deviation (sigma) for each object
and sample many perturbed relative positions at the nominal TCA. The collision probability
is the fraction of samples landing inside the combined hard-body radius.

This is the simplest of the three methods described in that paper (Monte Carlo). It does
NOT implement Patera's or Alfano's analytic methods; it is meant as a first, honest estimate
while actual position-uncertainty values are pending from the Space Systems Laboratory.

IMPORTANT: the default sigma_km and hard_body_radius_km below are PLACEHOLDERS, not
measured values. Treat every probability from this script as illustrative until real
values are confirmed.

Run:  python monte_carlo_probability.py                      (ISS vs Tiangong, 24 h window)
      python monte_carlo_probability.py "Hubble" "Terra" 12   (names from real_satellites.py)
Needs: numpy, scipy, matplotlib, sgp4
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sgp4.api import Satrec

from closest_approach import epoch_jd, separation, find_closest_approach

# --- PLACEHOLDERS: ask Partha / the lab what to use before quoting any probability ---
SIGMA_KM = 1.0              # assumed 1-sigma position error per object, each axis (independent, Gaussian)
HARD_BODY_RADIUS_KM = 0.05  # assumed combined physical size of the two objects (50 m)
N_TRIALS = 200_000


def monte_carlo_pc(r_rel_ca, sigma_km=SIGMA_KM, hard_body_km=HARD_BODY_RADIUS_KM,
                    n_trials=N_TRIALS, rng=None):
    """Monte Carlo collision probability given the nominal relative position at TCA.

    Each object's position error is modelled as independent N(0, sigma_km^2) per axis,
    so the relative position error has std sigma_km * sqrt(2) per axis.
    """
    rng = rng or np.random.default_rng()
    sigma_rel = sigma_km * np.sqrt(2)
    samples = r_rel_ca + rng.normal(0, sigma_rel, size=(n_trials, 3))
    dist = np.linalg.norm(samples, axis=1)
    hits = dist <= hard_body_km
    pc = hits.mean()
    se = np.sqrt(pc * (1 - pc) / n_trials)   # standard error of the MC estimate itself
    return pc, se, dist


def analyse(sat1, sat2, name1, name2, hours=24.0, sigma_km=SIGMA_KM,
            hard_body_km=HARD_BODY_RADIUS_KM, n_trials=N_TRIALS, out="monte_carlo_probability.png"):
    jd0 = max(epoch_jd(sat1), epoch_jd(sat2))
    t, d, t_ca, d_ca = find_closest_approach(sat1, sat2, jd0, hours)
    r_rel_ca, _ = separation(sat1, sat2, jd0, t_ca)
    r_rel_ca = r_rel_ca[0]

    pc, se, dist = monte_carlo_pc(r_rel_ca, sigma_km, hard_body_km, n_trials)

    print(f"{name1} vs {name2}")
    print(f"Nominal TCA: t = {t_ca:.1f} s after start, DCA = {d_ca:.3f} km (SGP4, no uncertainty)")
    print(f"Assumed sigma (per object, per axis): {sigma_km} km   |  hard-body radius: {hard_body_km} km")
    print(f"Monte Carlo collision probability at TCA: Pc = {pc:.3e}  (N = {n_trials:,} trials, MC std error {se:.1e})")
    print("These numbers depend entirely on the two assumed inputs above. They are placeholders")
    print("until the lab confirms real position-error and hard-body-radius values.")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ax[0].hist(dist, bins=80, color="tab:blue", alpha=0.8)
    ax[0].axvline(hard_body_km, color="red", ls="--", label=f"hard-body radius ({hard_body_km} km)")
    ax[0].set_xlabel("sampled miss distance at TCA (km)"); ax[0].set_ylabel("count")
    ax[0].set_title(f"Monte Carlo miss-distance distribution\n(Pc = {pc:.2e})"); ax[0].legend()

    sigmas = np.array([0.1, 0.25, 0.5, 1.0, 2.0, 5.0])
    pcs = np.array([monte_carlo_pc(r_rel_ca, s, hard_body_km, n_trials // 4)[0] for s in sigmas])
    # Plot real zeros honestly (as "< 1 hit observed") instead of hiding them under a log-scale floor.
    zero = pcs == 0
    if zero.any():
        min_detectable = 1.0 / (n_trials // 4)  # smallest Pc this many trials could even detect (1 hit)
        ax[1].semilogy(sigmas[~zero], pcs[~zero], "o-", color="tab:blue", label="Pc observed (>0 hits)")
        ax[1].scatter(sigmas[zero], [min_detectable] * zero.sum(), marker="v", color="tab:red",
                      label=f"0 hits in {n_trials//4:,} trials (Pc below this mark)")
        ax[1].legend(fontsize=7)
    else:
        ax[1].semilogy(sigmas, pcs, "o-")
    ax[1].set_xlabel("assumed sigma per object (km)"); ax[1].set_ylabel("Pc (log scale)")
    ax[1].set_title("Sensitivity: Pc vs assumed position-error sigma")
    plt.tight_layout(); plt.savefig(out, dpi=120)
    return pc, se, t_ca, d_ca


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
