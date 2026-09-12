# Thermo-Optic Fault Detection Benchmark for Silicon Microring Resonators

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21782606.svg)](https://doi.org/10.5281/zenodo.21782606)

This software package provides a reproducible benchmark framework for
evaluating sequential change detection algorithms applied to thermo-optic
fault monitoring in silicon microring resonators.

The framework includes: (1) a physically calibrated first-order compact
thermal model (k₀ = 72.9 pm/°C for Si SOI; k₀ = 11.7 pm/mW for 4H-SiC,
calibrated to Sun et al. arXiv:2506.15035); (2) a realistic colored-noise
simulation (Gaussian + 1/f + random-walk drift); (3) implementations of
three primary sequential change detection methods — fixed threshold,
CUSUM, and the W-formula — with an additional windowed GLR detector
evaluated separately via a window-size sweep (see the manuscript's
Discussion section for why GLR is not part of the main comparison); and
(4) four synthetic fault scenarios — heater efficiency loss, thermal
time constant increase, thermal drift, and contact degradation.

Detectors are compared via full Delay–False Alarm Rate operating
frontiers obtained by threshold sweep, rather than at single operating
points. Simulation-based results are complemented by a real-signal
validation using a genuine oscilloscope transient (see
`experimental_validation_corrected.md` and `data/STEP-20DB-V0.4-*.csv`),
provided by Q. Li (Carnegie Mellon University).

**The W = Q·D − T formulation is an empirical heuristic**; the variable
definitions are domain-specific adaptations. No theoretical derivation
from first principles is claimed.

## Repository structure

```
thermal_model_v3.py                    Core simulation, detectors, Monte Carlo sweep
digital_twin_sic.py                    SiC digital-twin model
run_real_detectors.py                  Runs the three detectors on the real transient
seven_detectors.py                     Earlier detector-count exploration (superseded; kept for history)
glr_sweep_exact.py                     Windowed GLR window-size sweep
experimental_validation_corrected.md   Write-up of the real-transient validation
microring_manuscript_v4_corrected.md   Corrected manuscript text
data/                                  Real oscilloscope data + digitized figure data
figures/                                Generated result figures
tests/test_detectors.py                Test suite (pytest)
```

## Installation and usage

See `INSTALL.md` for setup instructions and `CONTRIBUTING.md` for how to
propose changes.

## Dataset attribution

Raw oscilloscope data kindly provided by Q. Li (Carnegie Mellon
University), corresponding to Fig. 2(e) of Sun et al., arXiv:2506.15035.

## Keywords

microring resonator, silicon photonics, sequential change detection,
CUSUM, GLR, thermo-optic, fault detection, digital twin, health
monitoring

## Related work

This benchmark's simulation results are software-only; the accompanying
manuscript places these results relative to established practice
(industrial digital diagnostic monitoring, and the sequential-detection
literature) — see the manuscript's Discussion section.

## License

MIT License (see `LICENSE`).

## Citation

If you use this software, please cite the Zenodo deposit:
DOI: [10.5281/zenodo.21782606](https://doi.org/10.5281/zenodo.21782606)
