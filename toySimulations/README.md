# Toy simulation

The simulation in this folder help to understand the effect a hole in the detector acceptance causes to a mixed event distribution. The following files are included:

**acceptanceHoleMixing.C**

Toy simulation for the hole in acceptance. It is assumed that jets and tracks are uniformly distributed within their acceptances and that there is a hole, in which no tracks and no jets are seen.

**acceptanceHolePlotter.C**

Plotter macro for the above simulation.

**distributionMixing.C**

Create mixed event from measured jet and track eta-phi maps.

**distributionMixingPlotter.C**

Plotter macro for the above simulation.

**toySmearing.C**

Toy simulation to study the effect of jet pT smearing to jet spectra.

**toySmearingPlotter.C**

Plotter macro for the above simulation.
