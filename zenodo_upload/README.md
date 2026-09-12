# Thermo-Optic Fault Detection Benchmark for Silicon Microring Resonators

**A reproducible benchmark framework for evaluating sequential change detection algorithms on a physically calibrated digital twin.**

## Authors

Dimitar Kretski  
Center for Hydro- and Aerodynamics (CHA), Varna, Bulgaria  
ORCID: [0000-0001-5108-2243](https://orcid.org/0000-0001-5108-2243)

## Overview

This repository contains the simulation code, detector implementations, and figures for the paper:

> "A Compact Thermo-Optic Health Monitoring Framework for Silicon Microring Resonators: Benchmark of Sequential Change Detection Methods"

The framework evaluates **seven sequential change detection algorithms** on a physically calibrated first-order thermal model of a silicon microring resonator:

1. Fixed threshold (µ + m·σ)
2. CUSUM (one-sided, Lorden 1971)
3. EWMA (exponentially weighted moving average)
4. SPRT (sequential probability ratio test, Wald 1947)
5. Shiryaev-Roberts
6. Kalman innovation detector
7. W-formula safety layer: W(t) = Q(t)·D(t) − T(t) [Kretski 2025, DOI: 10.5281/zenodo.19553825]

## Physical Model Parameters

Calibrated to Sun et al. arXiv:2506.15035 (4H-SiC microring, NiCr heater):

| Parameter | Value | Source |
|-----------|-------|--------|
| k₀ (heater efficiency) | 11.7 pm/mW | Sun et al. Fig. 2(d) |
| τ_th (thermal time constant) | 7 µs | Sun et al. Fig. 2(e) |
| k_TO (thermo-optic coefficient) | 26.8 pm/K | Sun et al. Eq. 2 |

Additional Si SOI parameters (Bahadori et al. JLT 2018):

| Parameter | Value | Source |
|-----------|-------|--------|
| k₀ | 72.9 pm/°C | Literature |
| τ_th | 7 µs | Literature [4] |

**All results are simulation-based.** Experimental validation on a fabricated device is the identified next step.

## Files

| File | Description |
|------|-------------|
| `thermal_model_v3.py` | Core Si SOI compact model, noise simulation, CUSUM + W-formula detectors, Monte Carlo sweep |
| `detectors_extended.py` | EWMA, SPRT, Shiryaev-Roberts, Recursive GLR, Kalman innovation, BOCPD |
| `digital_twin_sic.py` | SiC digital twin (Sun et al. parameters), 4 fault scenarios, all 7 detectors |
| `digital_twin_results.csv` | Summary results: Pd, delay, FAR for all detectors × all fault scenarios |
| `sun2025_fig2d.csv` | Digitized λ(P) curve from Sun et al. Fig. 2(d) |
| `sun2025_fig2e.csv` | Digitized thermal step response from Sun et al. Fig. 2(e) |
| `fig1_delay_far_frontiers.png` | Publication figure: Delay–FAR operating frontiers (Si SOI, fault ×2.2) |
| `fig2_all_detectors_frontier.png` | All 9 detectors compared on Delay–FAR space |
| `fig3_robustness_sweep.png` | Robustness across 6 noise configurations |

## Quickstart

```bash
pip install numpy scipy matplotlib

# Run Si SOI benchmark (200 trials × 5 fault severities)
python thermal_model_v3.py

# Run SiC digital twin (4 fault scenarios)
python digital_twin_sic.py
```

## Fault Scenarios (digital_twin_sic.py)

1. **efficiency_loss** — heater efficiency degrades (R_th × 2.0)
2. **tau_increase** — thermal time constant increases (R_th × 1.8, slower response)
3. **thermal_drift** — slow ambient temperature drift
4. **contact_degradation** — random resistance spikes (brief power bursts)

## Key Results

At canonical operating points (CUSUM k=2.0, h=6.0; Fixed µ+5σ; W 0.01th percentile):

- CUSUM achieves 7–17% shorter detection delay than fixed threshold at matched FAR
- W-formula occupies a distinct operating region (shorter delay, higher FAR)
- Shiryaev-Roberts is fastest for step-type faults at Pd≥0.95
- W-formula shows relative advantage for impulsive contact degradation
- No single detector dominates all fault scenarios

## Methodological Note

All results are from physically calibrated Monte Carlo simulation. The W-formula threshold is adaptive (pre-fault percentile calibration). Comparison uses operating frontiers (threshold sweeps), not single operating points.

## Related Work

W-formula framework: DOI [10.5281/zenodo.19553825](https://doi.org/10.5281/zenodo.19553825)  
ORAC-QNode (embedded W-formula for quantum hardware): [github.com/Kretski/ORAC-QNode](https://github.com/Kretski/ORAC-QNode)

## License

MIT License — see LICENSE file.

## Citation

If you use this code, please cite:

```bibtex
@software{kretski2026microring,
  author    = {Kretski, Dimitar},
  title     = {Thermo-Optic Fault Detection Benchmark for Silicon Microring Resonators},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.21782606},
  url       = {https://doi.org/10.5281/zenodo.21782606}
}
```
