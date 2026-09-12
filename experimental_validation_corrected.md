# Experimental Validation (corrected)

## Thermal time constant

The thermal time constant was fitted directly from the raw oscilloscope
step-response data (`STEP-20DB-V0.4-8.csv`, provided by Q. Li, CMU),
following the normalization procedure in the accompanying MATLAB script
(bias adjustment of −0.03, fit region restricted to t > 0):

    τ = 6.56 ± 0.03 µs   (R² = 0.966, n = 1744 points)

This is consistent with the device-average value reported in the original
publication (~7 µs) and with values obtained independently by the original
authors under alternate normalization choices (6.5 µs without bias
adjustment; 5.6 µs with an additional 1.05× scaling factor; Q. Li, personal
communication, 11 August 2026). An earlier draft of this section reported
τ = 4.78 µs, R² = 0.9966; that figure was obtained from a hand-digitized
reconstruction of Fig. 2(e) rather than from the raw measurement, and has
been superseded by the value above, which uses the actual oscilloscope
trace.

## Detection delay on the real transient

The four detectors implemented in the accompanying code
(`thermal_model_v3.py`: fixed threshold, CUSUM, windowed GLR, and the
W-formula detector), run with the same calibration used throughout the
Monte Carlo simulation results reported elsewhere in the manuscript
(CUSUM k=2.0σ, h=6.0σ; GLR window=20, threshold=3.5σ), were applied
directly to the real optical step response. All four detectors correctly
identified the thermal transient, with no false alarms on the pre-trigger
baseline and no missed detections:

| Detector                     | Delay (µs) |
|-------------------------------|-----------:|
| Fixed threshold (µ+5σ)        |      0.200 |
| CUSUM (k=2.0σ, h=6.0σ)         |      0.175 |
| Windowed GLR (window=20, 3.5σ) |      0.425 |
| W-formula (MicroSafe-RL)       |      0.200 |

The dissipated-power term Q(t) in the W-formula detector was constructed
from the real electrical drive signal (`STEP-20DB-V0.4-IN8.csv`),
resampled onto the optical channel's time base; the single missing sample
at the start of that file (an empty oscilloscope reading) was replaced
with the channel minimum before use.

An earlier draft of this section described seven detectors (CUSUM, EWMA,
SPRT, Shiryaev-Roberts, Kalman innovation, W-formula, and fixed threshold)
with delays of 0.05–0.18 µs. That description did not match the four
detectors actually implemented in the accompanying code (fixed, CUSUM,
windowed GLR, W-formula — no EWMA, SPRT, Shiryaev-Roberts, or Kalman
variant exists in the codebase) and has been corrected above to describe
only the methods that were actually run.

## Dataset attribution

Raw oscilloscope data kindly provided by Q. Li (Carnegie Mellon
University), corresponding to Fig. 2(e) of Sun et al., arXiv:2506.15035.
