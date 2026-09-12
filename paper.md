---
title: 'A Compact Thermo-Optic Fault Detection Framework for Silicon Microring Resonators'
tags:
  - Python
  - silicon photonics
  - microring resonators
  - sequential change detection
  - fault detection
  - CUSUM
  - digital twin
authors:
  - name: Dimitar Kretski
    orcid: 0000-0001-5108-2243
    affiliation: 1
affiliations:
  - name: Center for Hydro- and Aerodynamics (CHA), Institute of Chemical Engineering, Bulgarian Academy of Sciences, Varna, Bulgaria
    index: 1
date: 12 September 2026
bibliography: paper_microring.bib
---

# Summary

Silicon photonic microring resonators are core building blocks of optical
transceivers, sensing systems, and wavelength-division-multiplexed links.
Their resonance wavelength is highly sensitive to temperature through the
thermo-optic effect, making early detection of thermal faults (heater
degradation, drift, contact issues) directly relevant to system
reliability. `thermo-optic-benchmark` is a Python software package that
provides a reproducible simulation and evaluation framework for comparing
sequential change-point detection methods applied to this problem. The
package includes a physically parameterized first-order compact thermal
model, a colored-noise simulator (Gaussian, 1/f pink noise, and
random-walk drift), and implementations of three primary detectors
(fixed threshold, CUSUM, and a domain-adapted W-formula heuristic), with
an additional windowed generalized-likelihood-ratio (GLR) detector
evaluated separately through a window-size sweep. Detectors are compared
via full Delay-versus-False-Alarm-Rate operating frontiers, obtained by
sweeping each detector's decision threshold, rather than at a single,
potentially cherry-picked, operating point. The source code is available
at `github.com/Kretski/thermo-optic-benchmark`, with an archived release
on Zenodo [@kretski2025zenodo].

# Statement of need

Comparative studies of fault-detection methods for photonic thermal
monitoring often report performance at a single operating point, making
it difficult to distinguish a genuinely better detector from one that is
simply better-tuned for that specific point [@poor2008quickest]. This
package addresses that gap by sweeping detector thresholds to construct
full operating frontiers, and by validating the underlying thermal model
against literature-reported device parameters [@bahadori2018thermal] and,
separately, against a real oscilloscope measurement of a physical
thermal transient. The real-signal validation (fitted time constant
$\tau = 6.56 \pm 0.03\ \mu s$, consistent with device-averaged literature
values) demonstrates that the simulation-calibrated detectors transfer to
a real measured signal without modification, complementing — though not
replacing — the Monte Carlo simulation results that form the bulk of the
benchmark. The software is aimed at researchers and engineers evaluating
sequential detection methods for photonic thermal management who need a
transparent, reproducible starting point rather than a black-box
comparison.

# Real-signal validation

Raw oscilloscope data for a real thermal step response, corresponding to
Fig. 2(e) of @sun2025integrated, was kindly provided by the paper's
corresponding author for validation purposes. The three primary detectors
were applied to this real transient using the same calibration as the
Monte Carlo results, correctly identifying the thermal transient with no
false alarms on the pre-trigger baseline, at delays consistent in order
of magnitude with the simulated results at comparable fault severity.
This is reported as a single real measurement, not a statistically
powered experimental trial series, and is explicitly scoped as such in
the accompanying manuscript.

# Acknowledgements

Raw oscilloscope data kindly provided by Q. Li (Carnegie Mellon
University), corresponding to Fig. 2(e) of @sun2025integrated.

# References
