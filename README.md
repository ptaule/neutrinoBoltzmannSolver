# Exact Boltzmann solution for massive neutrinos

This repo implements analytic solutions for the Boltzmann eqs. for massive
neutrinos, as given in Shoji&Komatsu '08. Some assumptions are made:

- psi = phi for the grav. potentials
- The grav. potential is assumed constant during matter domination, which makes
  it cancel out in ratios of perturbations in this period
- At some redshift (e.g. z_lambda=5) the time dependence of the grav.potential is taken into account

The output directory includes results for various neutrino masses, using
z_lambda=5. The time-axis (rows) of the data in this directory is reversed so
that it corresponds to decreasing redshift (increasing eta_D = ln D).  Note
that for high k integrands become oscillatory, and it is difficult to achieve
sub-percent precision.
