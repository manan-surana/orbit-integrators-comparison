# Comparing Euler, RK4 and Leapfrog Integrators for Orbital Motion

## What this project does
Simulates a satellite in a circular orbit (units where GM = 1) using three
numerical methods and compares how well each conserves energy.

## Methods
- Euler (first-order)
- Runge-Kutta 4 (fourth-order)
- Leapfrog / velocity Verlet (second-order, symplectic)

## Results
![Energy error](energy_error.png)

- Euler's energy error grows quickly and the orbit spirals outward.
- RK4 is very accurate at first but its energy drifts slowly over long runs.
- Leapfrog's energy error stays bounded.
(Replace these with what YOUR plots actually show.)

## How to run
pip install numpy matplotlib
python orbit_integrators.py

## What I learned / next steps
- Add convergence tests (error vs step size).
- Compare against SciPy's solve_ivp.
- Test elliptical orbits.

## Author
[Your Name], Class 12, [School], Delhi
