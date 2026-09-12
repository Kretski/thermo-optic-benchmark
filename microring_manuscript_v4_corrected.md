# A Compact Thermo-Optic Fault Detection Framework for Silicon Microring Resonators

Dimitar Kretski
Center for Hydro- and Aerodynamics (CHA), Varna, Bulgaria
ORCID: 0000-0001-5108-2243

## Abstract

We present a compact thermo-optic fault detection framework for silicon
microring resonators comprising a physics-based first-order thermal model,
a realistic colored-noise simulation (Gaussian + 1/f + random-walk drift),
and a Monte Carlo evaluation methodology that compares detectors through
full Delay–False Alarm Rate (FAR) operating frontiers rather than isolated
operating points. Three structurally distinct detection strategies are
evaluated: fixed threshold, CUSUM, and the W-formula safety layer W(t) =
Q(t)·D(t) − T(t). Under N = 200 Monte Carlo trials per configuration, CUSUM
consistently improves upon the fixed-threshold detector. Within the
evaluated operating range, the W-formula detector achieves shorter
detection delays at the cost of higher FAR, occupying a distinct region of
the Delay–FAR space. A robustness sweep across six noise configurations
shows that the relative ranking of the detectors did not change across the
investigated conditions. **The simulation-based framework is complemented
by a first experimental validation: the thermal time constant and the
detection delay of all three methods were confirmed against a real optical
step-response measurement.** No claim of optimality is made; several
established methods — including BOCPD, Shiryaev–Roberts, and recursive
GLR — were not included in the main comparison (see Discussion).

## I. Introduction

*[unchanged from the submitted version]*

Silicon photonic microring resonators are core components of optical
transceivers, sensing systems, and wavelength-division multiplexed (WDM)
links. Their resonance wavelength shifts at approximately 70–80 pm/°C due
to the thermo-optic effect, making them sensitive to thermal disturbances
and susceptible to failure of thermal management components such as
thermo-electric coolers (TEC).

Early detection of thermal faults directly limits the consequences of
degradation: a faster alarm reduces thermal overshoot, shortens recovery
time, and prevents resonance excursion beyond the tracking range of an
active controller. The detection delay is therefore a first-order quality
metric — yet it must be balanced against the false alarm rate (FAR),
because spurious alarms disrupt system operation.

Prior comparative studies typically evaluate detectors at a single
operating point, making it difficult to distinguish between a
fundamentally superior detector and one that is simply better-tuned. We
address this by sweeping each detector's decision threshold to obtain a
full operating frontier in Delay–FAR space, following a methodology
adapted from sequential detection theory [5].

This paper contributes: (1) a parameterised compact thermo-optic model
with explicit provenance annotations; (2) a realistic colored-noise
simulation; (3) operating frontier evaluation of three structurally
distinct detectors; (4) adaptation of the deterministic W-formula safety
layer [2] — originally a domain-independent health metric — to
thermo-optic monitoring, with analysis of its Delay–FAR trade-off and
robustness across noise conditions; **and (5) a first experimental
validation of the model's thermal time constant and of all three
detectors' delay performance against a real measured transient.**

## II. Thermo-Optic Compact Model

*[unchanged]*

The resonance wavelength shift is modelled as Δλ = k(T)·(T − Tamb), with
k(T) = k₀ + (dk/dT)·(T − Tref). Thermal dynamics follow a first-order RC
model:

    Cth · dT/dt = Pdiss − (T − Tamb) / Rth

R_th and C_th are derived from the heater efficiency and thermal time
constant. The heater response supports arbitrary polynomial order. A
calibration function accepts measured (power, wavelength-shift) pairs and
returns R², RMSE, degrees of freedom, and per-parameter standard errors.

**Table I. Model parameters and provenance.**

| Parameter | Value | Unit | Provenance |
|---|---|---|---|
| k₀ | 72.9 | pm/°C | Literature [3] |
| dk/dT | 0.0 | pm/°C² | Estimated (no lit. value found) |
| Heater eff. | 100.0 | pm/mW | Literature range [3] |
| τ_th | 7 × 10⁻⁶ | s | Literature [4] |
| σ_Gaussian | 2.0 | pm | Estimated |
| σ_pink | 1.0 | pm | Estimated |
| σ_drift | 0.5 | pm/√s | Estimated |

## III. Noise Model and Simulation

*[unchanged]*

Each trial simulates a total time of 2.5 × fault_onset at 70 ns resolution
(τ_th / 100). A fault is introduced as a ramp in thermal resistance over
20 µs followed by a sustained step, parameterised by a multiplier on
R_th. The baseline noise model superimposes: (i) white Gaussian noise
(σ = 2 pm), (ii) 1/f pink noise via the Voss–McCartney algorithm
(amplitude 1 pm), and (iii) a random-walk drift (0.5 pm/√s). The drift
term is the principal driver of false alarms for fixed-threshold
detectors.

## IV. Detector Formulations

*[unchanged]*

### A. Fixed Threshold

Alarm when r(t) = T_implied(t) − T_model(t) exceeds µ_pre + m·σ_pre.
Multiplier m swept over {3, 4, 5, 6, 7, 8, 10}.

### B. CUSUM

S_i = max(0, S_{i−1} + (r_i − µ_pre) − k), alarm when S_i > h. Fixed
k_σ = 2.0; h_σ swept over {2, 3, 4, 5, 6, 7, 8, 10, 15}. CUSUM is
asymptotically optimal under i.i.d. Gaussian residuals (Page 1954 [1]).

*[Corrected: the submitted version cited "Lorden 1971" in-text against a
reference-list entry for Page 1954. Page's 1954 paper introduced the
CUSUM procedure; the asymptotic optimality result is due to Lorden
(1971), a separate paper not currently in the reference list. Either cite
Page 1954 for the procedure itself (as corrected above), or add Lorden's
1971 paper as its own numbered reference if the optimality claim
specifically is being attributed.]*

### C. W-Formula Safety Layer

The W-formula W = Q·D − T was originally proposed in [2] as a
domain-independent adaptive health metric. The two detectors belong to
different mathematical classes: CUSUM linearly accumulates signed
deviations of r(t); W(t) is a nonlinear function of drive power, model
fidelity, and cumulative thermal stress. They are not algebraic
reformulations of one another. The domain-specific variable redefinitions
for thermo-optic monitoring are:

    Q(t) = P_diss(t)/mean(P_diss)      [normalised drive power]
    D(t) = exp(−|r(t)|/σ_pre)          [model fidelity, ∈ (0,1]]
    T(t) = ∫|r|dt/(σ_pre·T_window)     [cumulative thermal stress]

Alarm when W(t) < w_thresh, where w_thresh is set at the p-th percentile
of the pre-fault W trace computed on each trial independently — making
the threshold adaptive to the noise level of that trial. Percentile p
swept over {0.01%, 0.05%, 0.1%, 0.25%, 0.5%, 1%, 2%, 5%, 10%}. The minimum
achievable FAR (~3.5 × 10⁻⁴) reflects the noise correlation structure
under this calibration strategy.

## V. Results

*[unchanged — Monte Carlo simulation results, Sections A–C exactly as
submitted]*

### A. Operating Frontiers (fault ×2.2)

... *(Fig. 1, Table II, as submitted)*

### B. Fault Severity Dependence

... *(Table III, as submitted)*

### C. Robustness Across Noise Configurations

... *(Fig. 2, Table IV, as submitted)*

---

## VI. Experimental Validation

The simulation-based framework above is complemented by a first
experimental validation using a real optical step-response measurement
(the transient underlying Fig. 2(e) of Sun et al. [4]), obtained directly
from the authors as raw oscilloscope data.

### A. Thermal time constant

The thermal time constant was fitted from the raw trace
(`STEP-20DB-V0.4-8.csv`), following the normalization procedure of the
original measurement (bias adjustment of −0.03, single-exponential fit
restricted to t > 0):

    τ = 6.56 ± 0.03 µs   (R² = 0.966, n = 1744 points)

This is consistent with the device-average value reported in [4]
(≈7 µs) and with values obtained by the original authors under alternate
normalization choices (6.5 µs without bias adjustment; 5.6 µs with an
additional 1.05× scaling factor; W. Sun / Q. Li, personal communication).
The literature value of τ_th = 7 µs used as the model's provenance
parameter (Table I) is retained; the fitted value above is reported as an
independent experimental cross-check, not used to recalibrate the
simulation parameters.

### B. Detection delay on the real transient

The three detectors evaluated in Section V (fixed threshold, CUSUM, and
the W-formula safety layer), calibrated identically to the Monte Carlo
results above (CUSUM k = 2.0σ, h = 6.0σ; fixed threshold µ+5σ), were
applied directly to the real optical step response. The W-formula
detector's drive-power term Q(t) was constructed from the real electrical
heater-drive channel (`STEP-20DB-V0.4-IN8.csv`), resampled onto the
optical channel's time base. All three detectors correctly identified the
thermal transient, with no false alarms on the pre-trigger baseline:

| Detector | Delay (µs) |
|---|---:|
| Fixed threshold (µ+5σ) | 0.200 |
| CUSUM (k=2.0σ, h=6.0σ) | 0.175 |
| W-formula (MicroSafe-RL) | 0.200 |

These delays are consistent in order of magnitude with the simulated
values at comparable severity (Table III), and preserve the same
qualitative ordering (CUSUM ≤ Fixed at matched configuration) observed in
the Monte Carlo results, on a single real measurement.

### C. Scope of this validation

This is a single real transient (n = 1), not a repeated experimental
trial series, and does not carry the statistical weight of the N = 200
Monte Carlo results reported in Section V. It demonstrates that the
model's calibration and the three detectors' qualitative behavior
transfer from simulation to a real measured signal; it does not establish
detection performance under real, repeated fault conditions, real sensor
noise statistics, or real device-to-device variation. Hardware-in-the-loop
validation across repeated real fault events remains the appropriate next
step (see Limitations).

---

## VII. Discussion

*(renumbered from VI)*

**Structural distinctness of W and CUSUM.** The two detectors belong to
different mathematical classes. CUSUM linearly accumulates signed
deviations of the residual r(t) and is asymptotically optimal under i.i.d.
Gaussian noise. W(t) = Q(t)·D(t) − T(t) is a nonlinear function of drive
power, an exponential model-fidelity term, and cumulative thermal stress.
The detectors are not algebraic reformulations of one another; no
algebraic transformation W = f(CUSUM) or CUSUM = g(W) is known. A formal
proof of non-equivalence in the measure-theoretic sense is beyond the
scope of this paper.

**CUSUM vs Fixed threshold.** The consistent improvement of CUSUM
confirms that sequential evidence accumulation is beneficial even under
colored noise. The improvement is modest (7–17% in delay at matched FAR)
because the 1/f and drift components are non-stationary, partially
degrading CUSUM's advantage relative to purely Gaussian conditions.

**W-formula operating region and adaptive calibration.** The W-formula is
not a direct replacement for CUSUM: it achieves shorter delays but at
higher FAR. Its threshold is calibrated per trial from the pre-fault W
trace, which means the FAR is partly determined by the noise correlation
structure rather than a fixed design parameter. The robustness of the
relative ranking across noise configurations is partly attributable to
this adaptive calibration.

**Scope of comparisons and optimality.** No claim of optimality is made
for any evaluated detector. Several established sequential change-point
methods — including Bayesian Online Changepoint Detection (BOCPD),
Shiryaev–Roberts, EWMA, and Kalman innovation monitors — were not
evaluated in this study.

*[A parameter sweep of the windowed GLR detector across window sizes
{5, 8, 10, 15, 20} at fixed severity (×2.2), using the same simulation
methodology as Sections III–V (N = 200 trials), shows that its Delay–FAR
operating point depends strongly on window choice:]*

| Window | Threshold (σ) | Pd | Delay (µs) | FAR |
|---:|---:|---:|---:|---:|
| 20 (as tested above) | 3.5 | 1.00 | 2.790 | 0 |
| 8 | 3.0 | 1.00 | 2.215 | 0 |
| 8 | 2.5 | 1.00 | 2.040 | 3.5×10⁻⁶ |
| 5 | 3.0 | 1.00 | 2.092 | 0 |
| 5 | 2.5 | 1.00 | 1.918 | 5.3×10⁻⁶ |

*[For reference, CUSUM (k=2.0σ, h=6.0σ) achieves 2.087 µs at FAR≈0 under
identical trials — reproducing Table II's reported 2.09 µs to within
Monte Carlo noise, confirming the sweep methodology is consistent with
the main results.]*

The originally reported window (20) underperforms CUSUM, as noted in the
accompanying code. However, smaller windows (≤8) bring GLR's Delay–FAR
operating point into direct competition with CUSUM — at window=8,
threshold=2.5σ, GLR achieves 2.04 µs at FAR=3.5×10⁻⁶, matching or
exceeding CUSUM's 2.09 µs at FAR≈0. This indicates that GLR's apparent
disadvantage in the main comparison is a consequence of the specific
window size tested, not a structural limitation of the method under this
noise model. A full operating-frontier sweep of windowed GLR (analogous
to Fig. 1 for the other three detectors) is identified as a natural
extension of this work.

Future work should include these methods to obtain a more complete
picture of the achievable Delay–FAR trade-off.

**Limitations.** The Monte Carlo results (Sections III–V) remain
simulation-based, with the thermal model calibrated to literature values
rather than measured device data. Section VI provides a first real-signal
cross-check of the fitted time constant and of detector delay on a single
transient; it does not substitute for repeated hardware-in-the-loop
trials under real fault conditions. FAR = 0 entries carry an upper bound
of 1.8 × 10⁻³ at N = 200. W-formula FAR stability is reported within the
investigated simulated conditions only.

## VIII. Conclusion

*(renumbered from VII)*

We have presented a compact thermo-optic fault detection framework for
silicon microring resonators. By sweeping detector thresholds to obtain
full Delay–FAR operating frontiers, we demonstrate that CUSUM improves
upon fixed-threshold detection consistently within the evaluated
operating range, while the W-formula safety layer extends the observable
frontier toward shorter detection delays. A robustness sweep across six
noise configurations shows that the relative ranking of the three
detectors did not change across the investigated conditions. **A first
experimental validation against a real optical step-response measurement
confirms the fitted thermal time constant (τ = 6.56 ± 0.03 µs, consistent
with the literature value and with the original authors' own
re-analysis) and shows that all three detectors correctly identify the
real thermal transient with delays consistent in order of magnitude with
the simulated results.** The presented framework provides a reproducible
benchmark, now with an initial real-signal cross-check, for future
hardware-in-the-loop studies and for comparison with additional detection
methods not evaluated here.

## References

[1] E. S. Page, "Continuous inspection schemes," Biometrika, vol. 41, pp.
100–115, 1954.

[2] D. Kretski, "ORAC-NT Embedded Framework," Zenodo, DOI:
10.5281/zenodo.19553825, 2025.

[3] A. Bahadori et al., "Thermal Heterogeneous Integration of Silicon
Photonics," J. Lightwave Technol., vol. 36, pp. 3919–3930, 2018.

[4] W. Sun, R. Wang, J. Li, H. Zhang, Z. Jia, and Q. Li, "Integrated
microheater on the 4H-silicon-carbide-on-insulator platform and its
applications," arXiv:2506.15035, 2025.

[5] H. V. Poor and O. Hadjiliadis, *Quickest Detection*. Cambridge
University Press, 2009.

[6] D. Kretski, "Thermo-Optic Fault Detection Benchmark for Silicon
Microring Resonators," Zenodo, DOI: 10.5281/zenodo.21782606, 2026.
[Code]

## Acknowledgements

Raw oscilloscope data kindly provided by Q. Li (Carnegie Mellon
University), corresponding to Fig. 2(e) of Sun et al., arXiv:2506.15035.
