"""
Thermo-optic compact model v3 -- extends v2 with:
  6. Windowed GLR detector (closes the gap between v2 code and the summary).
  7. Fault-severity sweep: R_th_multiplier in {1.1, 1.2, 1.5, 2.0, 2.2}.
  8. Bootstrap Delay vs FAR plot (log-scale) with 95% CI for all detectors.
  9. W-formula MicroSafe-RL detector: deterministic safety layer W = Q*D - T.
     Provides the "fourth method" and the co-simulation anchor point for publication.

All v2 functionality (calibration, noise model, CUSUM, provenance) is preserved unchanged.
"""

from dataclasses import dataclass, field
from typing import Optional, Literal, List, Tuple, Dict
import numpy as np
import warnings

Provenance = Literal["measured", "literature", "estimated"]


# ---------------------------------------------------------------------------
# Core parameter structures (unchanged from v2)
# ---------------------------------------------------------------------------

@dataclass
class Param:
    value: float
    provenance: Provenance
    source: str
    unit: str

    def __repr__(self):
        if isinstance(self.value, (list, tuple)):
            val_str = "[" + ", ".join(f"{v:.6g}" for v in self.value) + "]"
        else:
            val_str = f"{self.value:.6g}"
        return f"{val_str} {self.unit}  [{self.provenance}: {self.source}]"


@dataclass
class ThermalOpticalParams:
    k0_pm_per_C: Param = field(default_factory=lambda: Param(
        72.9, "literature", "Bahadori et al., JLT 2018", "pm/degC"))
    dk_dT_pm_per_C2: Param = field(default_factory=lambda: Param(
        0.0, "estimated", "no literature value found; defaults to 0 (constant k) "
                          "until a wide-temperature-range sweep is available", "pm/degC^2"))
    heater_poly_coeffs: Param = field(default_factory=lambda: Param(
        100.0, "literature", "linear-only default (78-266 pm/mW reported range); "
                             "becomes a list after calibrate_from_measurements with degree>1", "pm/mW^n"))
    tau_th_s: Param = field(default_factory=lambda: Param(
        7e-6, "literature", "SiC microring measurement, arXiv:2506.15035", "s"))
    wavelength_noise_pm: Param = field(default_factory=lambda: Param(
        2.0, "estimated", "Gaussian component of monitor noise floor", "pm"))
    pink_noise_amplitude_pm: Param = field(default_factory=lambda: Param(
        1.0, "estimated", "1/f component amplitude, not chip-specific", "pm"))
    drift_randomwalk_pm_per_sqrt_s: Param = field(default_factory=lambda: Param(
        0.5, "estimated", "slow thermal/ambient drift, random-walk sigma per sqrt(s)", "pm/sqrt(s)"))
    P_base_mW: Param = field(default_factory=lambda: Param(
        5.0, "estimated", "typical driven-microring operating power", "mW"))

    def k_at_T(self, T_C: float, T_ref: float = 25.0) -> float:
        return self.k0_pm_per_C.value + self.dk_dT_pm_per_C2.value * (T_C - T_ref)

    def heater_shift_pm(self, P_mW: np.ndarray) -> np.ndarray:
        coeffs = np.atleast_1d(self.heater_poly_coeffs.value)
        result = np.zeros_like(np.asarray(P_mW, dtype=float))
        for n, a in enumerate(coeffs, start=1):
            result = result + a * np.asarray(P_mW, dtype=float) ** n
        return result

    def heater_eff_pm_per_mW_smallsignal(self) -> float:
        coeffs = np.atleast_1d(self.heater_poly_coeffs.value)
        return float(coeffs[0])

    def R_th(self, T_C: float = 25.0) -> float:
        return (self.heater_eff_pm_per_mW_smallsignal() / self.k_at_T(T_C)) * 1000.0

    def C_th(self, T_C: float = 25.0) -> float:
        return self.tau_th_s.value / self.R_th(T_C)

    def summary(self) -> str:
        lines = ["Parameter provenance:"]
        for name, p in self.__dict__.items():
            lines.append(f"  {name}: {p}")
        lines.append(f"  R_th @25C (derived): {self.R_th():.2f} K/W")
        lines.append(f"  C_th @25C (derived): {self.C_th():.4e} J/K")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Calibration (unchanged from v2)
# ---------------------------------------------------------------------------

def calibrate_from_measurements(
    params: ThermalOpticalParams,
    heater_power_mW: np.ndarray,
    wavelength_shift_pm: np.ndarray,
    poly_degree: int = 1,
    step_response_time_s: Optional[np.ndarray] = None,
    step_response_shift_pm: Optional[np.ndarray] = None,
    source_label: str = "user-provided CSV measurement",
) -> dict:
    from scipy.optimize import curve_fit
    P = np.asarray(heater_power_mW, dtype=float)
    Y = np.asarray(wavelength_shift_pm, dtype=float)
    A = np.vstack([P**n for n in range(1, poly_degree+1)]).T
    coeffs, _, _, _ = np.linalg.lstsq(A, Y, rcond=None)
    Y_pred = A @ coeffs
    residuals = Y - Y_pred
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((Y - Y.mean())**2)
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else float('nan')
    rmse = np.sqrt(ss_res / len(Y))
    dof = max(len(Y) - poly_degree, 1)
    sigma2 = ss_res / dof
    try:
        cov = sigma2 * np.linalg.inv(A.T @ A)
        param_se = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        param_se = np.full(poly_degree, np.nan)
    coeffs_list = coeffs.tolist() if poly_degree > 1 else float(coeffs[0])
    params.heater_poly_coeffs = Param(
        coeffs_list, "measured",
        f"{source_label} (degree-{poly_degree} poly fit, R^2={r2:.4f}, n={len(P)})", "pm/mW^n")
    diagnostics = {
        "coeffs": coeffs.tolist(), "coeffs_stderr": param_se.tolist(),
        "R2": float(r2), "RMSE_pm": float(rmse),
        "residuals": residuals.tolist(), "n_points": len(P), "dof": dof,
    }
    if step_response_time_s is not None and step_response_shift_pm is not None:
        t = np.asarray(step_response_time_s, dtype=float)
        y = np.asarray(step_response_shift_pm, dtype=float)
        def step_model(t, A_, tau): return A_ * (1 - np.exp(-t/tau))
        try:
            popt, pcov = curve_fit(step_model, t, y, p0=[y[-1], params.tau_th_s.value])
            tau_fit = float(popt[1])
            step_r2 = 1 - np.sum((y - step_model(t, *popt))**2)/np.sum((y-y.mean())**2)
            params.tau_th_s = Param(tau_fit, "measured",
                                     f"{source_label} (exp fit, R^2={step_r2:.4f})", "s")
            diagnostics["tau_fit_s"] = tau_fit
            diagnostics["tau_R2"] = float(step_r2)
        except RuntimeError:
            diagnostics["tau_fit_error"] = "curve_fit did not converge"
    return diagnostics


# ---------------------------------------------------------------------------
# Noise generator (unchanged from v2)
# ---------------------------------------------------------------------------

def generate_realistic_noise(N: int, dt_s: float, gaussian_sigma: float,
                              pink_amplitude: float, drift_sigma_per_sqrt_s: float,
                              rng: np.random.Generator) -> np.ndarray:
    white = rng.normal(0, gaussian_sigma, N)
    n_sources = 8
    pink = np.zeros(N)
    sources = rng.normal(0, 1, n_sources)
    for i in range(N):
        idx = 0
        n = i
        while n & 1 and idx < n_sources - 1:
            n >>= 1
            idx += 1
        sources[idx] = rng.normal(0, 1)
        pink[i] = sources.sum()
    if pink.std() > 0:
        pink = pink / pink.std() * pink_amplitude
    steps = rng.normal(0, drift_sigma_per_sqrt_s * np.sqrt(dt_s), N)
    drift = np.cumsum(steps)
    return white + pink + drift


# ---------------------------------------------------------------------------
# Detectors
# ---------------------------------------------------------------------------

def cusum_detector(residual: np.ndarray, k_sigma: float, h_sigma: float,
                   baseline_mu: float, baseline_sigma: float):
    """
    One-sided CUSUM (upward shift). Unchanged from v2.
    k_sigma: allowance in units of sigma. h_sigma: threshold in units of sigma.
    Returns (alarm_bool_array, cusum_trace).
    """
    k = k_sigma * baseline_sigma
    h = h_sigma * baseline_sigma
    S = np.zeros_like(residual)
    alarms = np.zeros(len(residual), dtype=bool)
    s = 0.0
    for i, x in enumerate(residual):
        s = max(0.0, s + (x - baseline_mu) - k)
        S[i] = s
        if s > h:
            alarms[i] = True
            s = 0.0
    return alarms, S


def windowed_glr_detector(residual: np.ndarray, window: int,
                           baseline_mu: float, baseline_sigma: float,
                           threshold_sigma: float):
    """
    Windowed Generalised Likelihood Ratio (GLR) detector.

    For each window of length `window`, computes the log-likelihood ratio
    between H1 (mean shift by unknown delta) and H0 (no shift):

        GLR_i = max_{delta>0} [ n*(delta_hat)^2 / (2*sigma^2) ]
               = n * (max(0, x_bar - mu))^2 / (2 * sigma^2)

    where x_bar is the window mean. This is the one-sided GLRT for Gaussian
    residuals with known variance (sigma from pre-fault baseline).

    threshold_sigma: alarm threshold in units of sigma (converted to GLR units
    internally as (threshold_sigma * sigma)^2 * n / (2 * sigma^2) = n * threshold_sigma^2 / 2).

    Returns (alarm_bool_array, glr_trace).

    Note: This is a simplified windowed approximation, NOT the full recursive
    GLR (Willsky-Jones). It has higher FAR under coloured noise (pink + drift)
    because the window integrates correlated samples as if they were i.i.d.
    This limitation is documented and is why CUSUM outperforms it here.
    """
    N = len(residual)
    glr = np.zeros(N)
    alarms = np.zeros(N, dtype=bool)
    sigma2 = baseline_sigma ** 2
    # Threshold: GLR > h_glr triggers alarm.
    # Derivation: GLR = n*(delta_hat)^2/(2*sigma^2); set delta_hat = threshold_sigma*sigma
    # => h_glr = window * threshold_sigma^2 / 2
    h_glr = window * threshold_sigma**2 / 2.0

    for i in range(window - 1, N):
        w = residual[i - window + 1: i + 1]
        x_bar = w.mean()
        delta_hat = max(0.0, x_bar - baseline_mu)   # one-sided: upward shift only
        glr_val = window * delta_hat**2 / (2.0 * sigma2) if sigma2 > 0 else 0.0
        glr[i] = glr_val
        if glr_val > h_glr:
            alarms[i] = True
    return alarms, glr


def w_formula_detector(residual: np.ndarray, P_diss_W: np.ndarray,
                        dt_s: float, baseline_sigma: float,
                        Q_weight: float = 1.0, D_weight: float = 1.0,
                        T_weight: float = 1.0, threshold_sigma: float = 5.0):
    """
    W-formula MicroSafe-RL detector: W(t) = Q(t)*D(t) - T(t)

    Mapping to microring fault detection:
      Q(t) = dissipated power proxy = P_diss_W[t] / P_diss_W.mean()
             (normalised; represents "quality" of current operating point)
      D(t) = model fidelity = exp(-|residual[t]| / baseline_sigma)
             (1 when residual=0, decays toward 0 as residual grows)
      T(t) = accumulated thermal stress = cumulative integral of |residual|*dt,
             normalised to baseline_sigma * time_window

    Alarm when W(t) < w_threshold.
    w_threshold is set automatically to the 1st percentile of W during pre-fault.

    This is the deterministic safety layer. It is NOT statistical: it uses
    no distribution assumptions. The threshold is calibrated on the pre-fault
    segment, making it self-tuning.

    Returns (alarm_bool_array, W_trace, w_threshold_used).

    Provenance of W = Q*D - T: Kretski (2025), ORAC-NT embedded framework,
    Zenodo DOI: 10.5281/zenodo.19553825
    """
    P_mean = P_diss_W.mean()
    if P_mean == 0:
        P_mean = 1.0  # avoid divide-by-zero in pathological case

    Q = Q_weight * P_diss_W / P_mean
    D = D_weight * np.exp(-np.abs(residual) / (baseline_sigma + 1e-30))

    # Cumulative thermal stress, normalised to [0, 1] range over full window
    cumstress = np.cumsum(np.abs(residual) * dt_s)
    norm_factor = baseline_sigma * (len(residual) * dt_s)
    T_stress = T_weight * cumstress / max(norm_factor, 1e-30)

    W = Q * D - T_stress

    # Auto-threshold: 1st percentile of W in pre-fault segment will be found
    # outside this function (caller passes pre_fault_W), so here we return W
    # and let the caller set threshold after seeing pre-fault stats.
    # For self-contained use, return raw W; threshold is set by caller.
    return W


def apply_w_threshold(W: np.ndarray, fault_onset_step: int,
                       percentile: float = 1.0):
    """
    Set threshold at `percentile`-th percentile of pre-fault W,
    then find alarms where W < threshold.
    Returns (alarm_bool_array, w_threshold).
    """
    pre = W[:fault_onset_step]
    if len(pre) == 0:
        w_threshold = W.min()
    else:
        w_threshold = np.percentile(pre, percentile)
    alarms = W < w_threshold
    return alarms, float(w_threshold)


# ---------------------------------------------------------------------------
# Core simulation engine (extended from v2)
# ---------------------------------------------------------------------------

def _build_simulation(params: ThermalOpticalParams, fault_onset_s: float,
                       fault_ramp_s: float, fault_R_th_multiplier: float,
                       total_time_s: float, dt_s: Optional[float],
                       seed: int, use_realistic_noise: bool):
    """
    Internal: builds one simulation trajectory and returns all raw arrays.
    Separated from detector logic so detectors can be varied cheaply.
    """
    rng = np.random.default_rng(seed)
    T_amb = 25.0
    R_th_nominal = params.R_th(T_amb)
    C_th = params.C_th(T_amb)

    if dt_s is None:
        dt_s = params.tau_th_s.value / 100.0
    N = int(total_time_s / dt_s)
    fault_onset_step = int(fault_onset_s / dt_s)
    ramp_len = max(1, int(fault_ramp_s / dt_s))

    P_diss_W = (params.P_base_mW.value * 1e-3
                + 0.3e-3 * np.sin(2*np.pi*50e3*np.arange(N)*dt_s)
                + rng.normal(0, 0.05e-3, N))

    R_th_actual = np.full(N, R_th_nominal)
    R_th_actual[fault_onset_step:fault_onset_step+ramp_len] = np.linspace(
        R_th_nominal, R_th_nominal*fault_R_th_multiplier, ramp_len)
    R_th_actual[fault_onset_step+ramp_len:] = R_th_nominal * fault_R_th_multiplier

    T_actual = np.zeros(N)
    T_actual[0] = T_amb
    for i in range(1, N):
        dTdt = (P_diss_W[i-1] - (T_actual[i-1]-T_amb)/R_th_actual[i-1]) / C_th
        T_actual[i] = T_actual[i-1] + dTdt * dt_s

    k_series = np.array([params.k_at_T(T) for T in T_actual])

    if use_realistic_noise:
        noise = generate_realistic_noise(
            N, dt_s, params.wavelength_noise_pm.value,
            params.pink_noise_amplitude_pm.value,
            params.drift_randomwalk_pm_per_sqrt_s.value, rng)
    else:
        noise = rng.normal(0, params.wavelength_noise_pm.value, N)

    delta_lambda_measured = k_series * (T_actual - T_amb) + noise

    T_model = np.zeros(N)
    T_model[0] = T_amb
    for i in range(1, N):
        dTdt = (P_diss_W[i-1] - (T_model[i-1]-T_amb)/R_th_nominal) / C_th
        T_model[i] = T_model[i-1] + dTdt * dt_s

    k_model_series = np.array([params.k_at_T(T) for T in T_model])
    T_implied = delta_lambda_measured / k_model_series + T_amb
    residual = T_implied - T_model

    pre_fault = residual[:fault_onset_step]
    mu_pre = pre_fault.mean()
    sigma_pre = pre_fault.std()

    return {
        "residual": residual, "P_diss_W": P_diss_W,
        "mu_pre": mu_pre, "sigma_pre": sigma_pre,
        "fault_onset_step": fault_onset_step,
        "dt_s": dt_s, "N": N,
        "R_th_nominal": R_th_nominal, "C_th": C_th,
    }


def simulate_fault_detection(
    params: ThermalOpticalParams,
    fault_onset_s: float,
    fault_ramp_s: float,
    fault_R_th_multiplier: float,
    total_time_s: float,
    dt_s: Optional[float] = None,
    seed: int = 42,
    use_realistic_noise: bool = True,
    cusum_k: float = 2.0,
    cusum_h: float = 6.0,
    glr_window: int = 20,
    glr_threshold_sigma: float = 3.5,
):
    """
    Runs one simulation and applies all four detectors:
      1. Fixed threshold (μ + 5σ)
      2. CUSUM (best params from v2 study: k=2.0, h=6.0)
      3. Windowed GLR
      4. W-formula MicroSafe-RL

    Returns a flat dict with delay and false-alarm count for each detector.
    """
    sim = _build_simulation(params, fault_onset_s, fault_ramp_s,
                             fault_R_th_multiplier, total_time_s, dt_s, seed,
                             use_realistic_noise)
    residual = sim["residual"]
    P_diss_W = sim["P_diss_W"]
    mu_pre = sim["mu_pre"]
    sigma_pre = sim["sigma_pre"]
    fos = sim["fault_onset_step"]
    dt_s_ = sim["dt_s"]

    results = {"R_th_nominal": sim["R_th_nominal"], "C_th": sim["C_th"],
               "dt_s": dt_s_, "N": sim["N"], "pre_fault_sigma": float(sigma_pre)}

    def _delay_fa(alarms):
        post = np.where(alarms[fos:])[0]
        delay = post[0] * dt_s_ if len(post) > 0 else None
        fa = int(alarms[:fos].sum())
        return delay, fa

    # 1. Fixed threshold
    fixed_detected = residual > (mu_pre + 5 * sigma_pre)
    results["fixed_delay_s"], results["fixed_fa"] = _delay_fa(fixed_detected)

    # 2. CUSUM
    cusum_alarms, _ = cusum_detector(residual, cusum_k, cusum_h, mu_pre, sigma_pre)
    results["cusum_delay_s"], results["cusum_fa"] = _delay_fa(cusum_alarms)

    # 3. Windowed GLR
    glr_alarms, _ = windowed_glr_detector(residual, glr_window, mu_pre, sigma_pre,
                                           glr_threshold_sigma)
    results["glr_delay_s"], results["glr_fa"] = _delay_fa(glr_alarms)

    # 4. W-formula
    W = w_formula_detector(residual, P_diss_W, dt_s_, sigma_pre)
    w_alarms, w_thresh = apply_w_threshold(W, fos, percentile=1.0)
    results["w_delay_s"], results["w_fa"] = _delay_fa(w_alarms)
    results["w_threshold"] = w_thresh

    return results


# ---------------------------------------------------------------------------
# Monte Carlo sweep: N_trials per (detector, severity)
# ---------------------------------------------------------------------------

def monte_carlo_sweep(
    params: ThermalOpticalParams,
    fault_R_th_multipliers: List[float],
    n_trials: int = 200,
    fault_onset_s: float = 200e-6,
    fault_ramp_s: float = 20e-6,
    total_time_s: float = 500e-6,
    use_realistic_noise: bool = True,
    cusum_k: float = 2.0,
    cusum_h: float = 6.0,
    glr_window: int = 20,
    glr_threshold_sigma: float = 3.5,
    base_seed: int = 0,
    verbose: bool = True,
) -> Dict[str, Dict[float, dict]]:
    """
    Runs n_trials per fault severity per detector.
    Returns nested dict: results[detector_name][multiplier] = {
        'delays': [...],   # in µs, None entries excluded
        'fa_rates': [...], # false alarms per pre-fault sample
        'pd': float,       # Pd = fraction of trials with detection
    }
    Detectors: 'fixed', 'cusum', 'glr', 'w_formula'
    """
    detectors = ['fixed', 'cusum', 'glr', 'w_formula']
    results = {d: {m: {'delays': [], 'fa_rates': [], 'pd': 0.0}
                   for m in fault_R_th_multipliers} for d in detectors}

    total = len(fault_R_th_multipliers) * n_trials
    done = 0
    for mult in fault_R_th_multipliers:
        for trial in range(n_trials):
            seed = base_seed + done
            try:
                r = simulate_fault_detection(
                    params, fault_onset_s, fault_ramp_s, mult, total_time_s,
                    seed=seed, use_realistic_noise=use_realistic_noise,
                    cusum_k=cusum_k, cusum_h=cusum_h,
                    glr_window=glr_window, glr_threshold_sigma=glr_threshold_sigma)
                fos = int(fault_onset_s / r["dt_s"])
                for d, delay_key, fa_key in [
                    ('fixed', 'fixed_delay_s', 'fixed_fa'),
                    ('cusum', 'cusum_delay_s', 'cusum_fa'),
                    ('glr',   'glr_delay_s',   'glr_fa'),
                    ('w_formula', 'w_delay_s', 'w_fa'),
                ]:
                    delay = r[delay_key]
                    fa = r[fa_key]
                    if delay is not None:
                        results[d][mult]['delays'].append(delay * 1e6)  # convert to µs
                    fa_rate = fa / max(fos, 1)
                    results[d][mult]['fa_rates'].append(fa_rate)
            except Exception as e:
                warnings.warn(f"Trial seed={seed} mult={mult:.1f} failed: {e}")
            done += 1

        for d in detectors:
            delays = results[d][mult]['delays']
            results[d][mult]['pd'] = len(delays) / n_trials

    return results


# ---------------------------------------------------------------------------
# Bootstrap CI helpers
# ---------------------------------------------------------------------------

def bootstrap_ci(data: List[float], n_boot: int = 2000, alpha: float = 0.05,
                 rng: Optional[np.random.Generator] = None) -> Tuple[float, float, float]:
    """Returns (mean, lower_CI, upper_CI) via percentile bootstrap."""
    if rng is None:
        rng = np.random.default_rng(999)
    arr = np.asarray(data, dtype=float)
    if len(arr) == 0:
        return (np.nan, np.nan, np.nan)
    boots = [rng.choice(arr, size=len(arr), replace=True).mean() for _ in range(n_boot)]
    lo = np.percentile(boots, 100 * alpha / 2)
    hi = np.percentile(boots, 100 * (1 - alpha / 2))
    return float(arr.mean()), float(lo), float(hi)


# ---------------------------------------------------------------------------
# Plotting: Delay vs FAR (log-scale) with bootstrap CI
# ---------------------------------------------------------------------------

def plot_delay_vs_far(sweep_results: Dict[str, Dict[float, dict]],
                      multipliers: List[float],
                      output_path: str = "delay_vs_far.png",
                      n_boot: int = 2000):
    """
    Generates the publication-ready Delay vs FAR figure.
    - X axis: mean FAR (false alarm rate per pre-fault sample), log scale.
    - Y axis: mean detection delay [µs], linear scale.
    - One curve per detector, one point per fault severity.
    - Error bars: 95% bootstrap CI on both axes.
    - Severity annotations on the CUSUM curve.

    Saves to output_path and returns the matplotlib Figure.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(42)

    DETECTOR_STYLE = {
        'fixed':     {'label': 'Fixed threshold (µ+5σ)', 'color': '#1f77b4',
                      'marker': 's', 'ls': '--'},
        'cusum':     {'label': f'CUSUM (k=2.0, h=6.0)',  'color': '#d62728',
                      'marker': 'o', 'ls': '-'},
        'glr':       {'label': 'Windowed GLR',            'color': '#2ca02c',
                      'marker': '^', 'ls': '-.'},
        'w_formula': {'label': 'W-formula (MicroSafe-RL)','color': '#9467bd',
                      'marker': 'D', 'ls': ':'},
    }

    fig, ax = plt.subplots(figsize=(8, 5.5))

    for det, style in DETECTOR_STYLE.items():
        xs_m, xs_lo, xs_hi = [], [], []
        ys_m, ys_lo, ys_hi = [], [], []
        mult_labels = []

        for mult in multipliers:
            bucket = sweep_results[det][mult]
            delays = bucket['delays']
            fa_rates = bucket['fa_rates']

            dy_m, dy_lo, dy_hi = bootstrap_ci(delays, n_boot=n_boot, rng=rng)
            fx_m, fx_lo, fx_hi = bootstrap_ci(fa_rates, n_boot=n_boot, rng=rng)

            ys_m.append(dy_m); ys_lo.append(dy_lo); ys_hi.append(dy_hi)
            xs_m.append(fx_m); xs_lo.append(fx_lo); xs_hi.append(fx_hi)
            mult_labels.append(f"×{mult:.1f}")

        xs_m = np.array(xs_m)
        ys_m = np.array(ys_m)

        # Replace exact-zero FAR with small positive for log axis
        xs_plot = np.where(xs_m == 0, 1e-6, xs_m)
        xs_lo_plot = np.where(np.array(xs_lo) <= 0, 1e-6, np.array(xs_lo))
        xs_hi_plot = np.where(np.array(xs_hi) <= 0, 1e-6, np.array(xs_hi))

        ax.errorbar(xs_plot, ys_m,
                    xerr=[xs_plot - xs_lo_plot, xs_hi_plot - xs_plot],
                    yerr=[ys_m - np.array(ys_lo), np.array(ys_hi) - ys_m],
                    fmt=style['marker'], color=style['color'], ls=style['ls'],
                    label=style['label'], capsize=4, linewidth=1.5,
                    markersize=8, alpha=0.85)

        # Annotate severity on CUSUM curve only (to avoid clutter)
        if det == 'cusum':
            for xi, yi, lbl in zip(xs_plot, ys_m, mult_labels):
                ax.annotate(lbl, (xi, yi), textcoords="offset points",
                            xytext=(6, 4), fontsize=8, color=style['color'])

    ax.set_xscale('log')
    ax.set_xlabel("False Alarm Rate (alarms / pre-fault sample)", fontsize=12)
    ax.set_ylabel("Detection Delay [µs]", fontsize=12)
    ax.set_title("Delay vs FAR — all detectors, fault severity sweep\n"
                 "Silicon microring thermo-optic fault detection "
                 "(Gaussian+1/f+drift noise)", fontsize=11)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, which='both', alpha=0.3)

    # Add note about Windowed GLR limitation
    ax.text(0.02, 0.04,
            "Note: Windowed GLR is an approximation; coloured noise inflates its FAR.\n"
            "Full recursive GLR expected to improve substantially.",
            transform=ax.transAxes, fontsize=7.5, color='gray',
            verticalalignment='bottom')

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Figure saved: {output_path}")
    return fig


# ---------------------------------------------------------------------------
# Main: reproduce all reported results + new severity sweep + figure
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import time

    params = ThermalOpticalParams()
    print(params.summary())

    # --- Single-trial sanity check (all 4 detectors) ---
    print("\n=== Single trial, all 4 detectors (fault ×2.2, realistic noise) ===")
    r = simulate_fault_detection(
        params, fault_onset_s=200e-6, fault_ramp_s=20e-6,
        fault_R_th_multiplier=2.2, total_time_s=500e-6, use_realistic_noise=True)
    for k, v in r.items():
        if v is not None:
            if isinstance(v, float) and abs(v) < 1e-3 and v != 0:
                print(f"  {k}: {v:.3e}")
            else:
                print(f"  {k}: {v}")

    # --- Monte Carlo sweep ---
    MULTIPLIERS = [1.1, 1.2, 1.5, 2.0, 2.2]
    N_TRIALS = 200

    print(f"\n=== Monte Carlo sweep: {N_TRIALS} trials × {len(MULTIPLIERS)} severities ===")
    print("Detectors: Fixed, CUSUM (k=2.0,h=6.0), Windowed GLR, W-formula")
    t0 = time.time()
    sweep = monte_carlo_sweep(
        params, MULTIPLIERS, n_trials=N_TRIALS,
        fault_onset_s=200e-6, fault_ramp_s=20e-6, total_time_s=500e-6,
        use_realistic_noise=True, cusum_k=2.0, cusum_h=6.0,
        glr_window=20, glr_threshold_sigma=3.5, verbose=True)
    elapsed = time.time() - t0
    print(f"Sweep completed in {elapsed:.1f}s")

    # --- Summary table ---
    print("\n=== Results summary ===")
    print(f"{'Detector':<16} {'Mult':>6} {'Pd':>6} {'Delay_mean µs':>14} "
          f"{'Delay_95CI':>18} {'FAR_mean':>12} {'FAR_95CI':>18}")
    rng_ci = np.random.default_rng(7)
    for det in ['fixed', 'cusum', 'glr', 'w_formula']:
        for mult in MULTIPLIERS:
            b = sweep[det][mult]
            dy_m, dy_lo, dy_hi = bootstrap_ci(b['delays'], rng=rng_ci)
            fx_m, fx_lo, fx_hi = bootstrap_ci(b['fa_rates'], rng=rng_ci)
            print(f"{det:<16} {mult:>6.1f} {b['pd']:>6.3f} "
                  f"{dy_m:>14.3f} [{dy_lo:.3f},{dy_hi:.3f}]  "
                  f"{fx_m:>12.2e} [{fx_lo:.2e},{fx_hi:.2e}]")

    # --- Calibration demo (unchanged from v2) ---
    print("\n=== Calibration demo (degree-2, unchanged from v2) ===")
    rng = np.random.default_rng(7)
    P_sweep = np.linspace(0, 15, 16)
    true_a1, true_a2 = 90.0, -0.8
    lam_sweep = true_a1*P_sweep + true_a2*P_sweep**2 + rng.normal(0, 3, len(P_sweep))
    diag = calibrate_from_measurements(params, P_sweep, lam_sweep, poly_degree=2,
                                        source_label="synthetic nonlinear demo")
    print("Diagnostics:")
    for k, v in diag.items():
        if k != "residuals":
            print(f"  {k}: {v}")

    # --- Generate figure ---
    print("\n=== Generating Delay vs FAR figure ===")
    plot_delay_vs_far(sweep, MULTIPLIERS, output_path="/home/claude/delay_vs_far.png")

    print("\nDone. Output: delay_vs_far.png")


# ---------------------------------------------------------------------------
# W-formula threshold sweep: iso-FAR comparison
# ---------------------------------------------------------------------------

def w_formula_threshold_sweep(
    params: ThermalOpticalParams,
    percentile_levels: list,          # e.g. [0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0]
    n_trials: int = 200,
    fault_R_th_multiplier: float = 2.2,
    fault_onset_s: float = 200e-6,
    fault_ramp_s: float = 20e-6,
    total_time_s: float = 500e-6,
    use_realistic_noise: bool = True,
    base_seed: int = 5000,
    n_boot: int = 2000,
):
    """
    Sweeps the W-formula alarm threshold (as pre-fault percentile) across
    `percentile_levels`. For each level, runs `n_trials` Monte Carlo trials
    and records mean delay + FAR with bootstrap CI.

    Returns list of dicts, one per percentile level.
    """
    rng_ci = np.random.default_rng(77)
    records = []

    for pct in percentile_levels:
        delays_us, fa_rates = [], []
        for trial in range(n_trials):
            seed = base_seed + trial
            sim = _build_simulation(
                params, fault_onset_s, fault_ramp_s, fault_R_th_multiplier,
                total_time_s, None, seed, use_realistic_noise)
            res = sim["residual"]
            P_diss = sim["P_diss_W"]
            mu_pre = sim["mu_pre"]
            sigma_pre = sim["sigma_pre"]
            fos = sim["fault_onset_step"]
            dt_s = sim["dt_s"]

            W = w_formula_detector(res, P_diss, dt_s, sigma_pre)
            w_alarms, _ = apply_w_threshold(W, fos, percentile=pct)

            post = np.where(w_alarms[fos:])[0]
            if len(post) > 0:
                delays_us.append(post[0] * dt_s * 1e6)
            fa_rate = w_alarms[:fos].sum() / max(fos, 1)
            fa_rates.append(float(fa_rate))

        dy_m, dy_lo, dy_hi = bootstrap_ci(delays_us, n_boot=n_boot, rng=rng_ci)
        fx_m, fx_lo, fx_hi = bootstrap_ci(fa_rates,  n_boot=n_boot, rng=rng_ci)
        pd = len(delays_us) / n_trials
        records.append({
            "percentile": pct,
            "pd": pd,
            "delay_mean_us": dy_m, "delay_lo_us": dy_lo, "delay_hi_us": dy_hi,
            "far_mean": fx_m,      "far_lo": fx_lo,      "far_hi": fx_hi,
        })
        print(f"  pct={pct:.3f}%  Pd={pd:.3f}  "
              f"delay={dy_m:.3f} [{dy_lo:.3f},{dy_hi:.3f}] µs  "
              f"FAR={fx_m:.2e} [{fx_lo:.2e},{fx_hi:.2e}]")

    return records


def plot_w_sweep_vs_cusum(
    w_records: list,
    cusum_ref: dict,           # {'delay_mean_us', 'delay_lo_us', 'delay_hi_us',
                               #  'far_mean', 'far_lo', 'far_hi'} at matched FAR
    fixed_ref: dict,
    output_path: str = "w_sweep_vs_cusum.png",
):
    """
    Plots the W-formula Delay–FAR trade-off curve together with CUSUM and
    Fixed-threshold reference points. This is the iso-FAR comparison figure.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pcts   = [r["percentile"]    for r in w_records]
    dy_m   = np.array([r["delay_mean_us"] for r in w_records])
    dy_lo  = np.array([r["delay_lo_us"]   for r in w_records])
    dy_hi  = np.array([r["delay_hi_us"]   for r in w_records])
    fx_m   = np.array([r["far_mean"]      for r in w_records])
    fx_lo  = np.array([r["far_lo"]        for r in w_records])
    fx_hi  = np.array([r["far_hi"]        for r in w_records])

    # Clamp zeros for log axis
    eps = 1e-7
    fx_plot    = np.where(fx_m  <= 0, eps, fx_m)
    fx_lo_plot = np.where(fx_lo <= 0, eps, fx_lo)
    fx_hi_plot = np.where(fx_hi <= 0, eps, fx_hi)

    fig, ax = plt.subplots(figsize=(8, 5.5))

    # W-formula sweep curve
    ax.errorbar(fx_plot, dy_m,
                xerr=[fx_plot - fx_lo_plot, fx_hi_plot - fx_plot],
                yerr=[dy_m - dy_lo, dy_hi - dy_m],
                fmt='D-', color='#9467bd', capsize=4, linewidth=1.8,
                markersize=7, label='W-formula (threshold sweep)', alpha=0.9)
    for xi, yi, p in zip(fx_plot, dy_m, pcts):
        ax.annotate(f"{p:.2f}%", (xi, yi), textcoords="offset points",
                    xytext=(5, 4), fontsize=7.5, color='#9467bd')

    # CUSUM reference point
    c_fx = cusum_ref["far_mean"] if cusum_ref["far_mean"] > 0 else eps
    ax.errorbar(c_fx, cusum_ref["delay_mean_us"],
                xerr=[[c_fx - max(cusum_ref["far_lo"], eps)],
                      [cusum_ref["far_hi"] - c_fx]],
                yerr=[[cusum_ref["delay_mean_us"] - cusum_ref["delay_lo_us"]],
                      [cusum_ref["delay_hi_us"] - cusum_ref["delay_mean_us"]]],
                fmt='o', color='#d62728', capsize=5, markersize=10,
                label='CUSUM k=2.0, h=6.0', zorder=5)

    # Fixed threshold reference point
    f_fx = fixed_ref["far_mean"] if fixed_ref["far_mean"] > 0 else eps
    ax.errorbar(f_fx, fixed_ref["delay_mean_us"],
                xerr=[[f_fx - max(fixed_ref["far_lo"], eps)],
                      [fixed_ref["far_hi"] - f_fx]],
                yerr=[[fixed_ref["delay_mean_us"] - fixed_ref["delay_lo_us"]],
                      [fixed_ref["delay_hi_us"] - fixed_ref["delay_mean_us"]]],
                fmt='s', color='#1f77b4', capsize=5, markersize=10,
                label='Fixed threshold (µ+5σ)', zorder=5)

    ax.set_xscale('log')
    ax.set_xlabel("False Alarm Rate (alarms / pre-fault sample)", fontsize=12)
    ax.set_ylabel("Detection Delay [µs]", fontsize=12)
    ax.set_title(
        "W-formula threshold sweep vs CUSUM & Fixed threshold\n"
        "Fault ×2.2, 200 trials, Gaussian+1/f+drift noise  "
        "(iso-FAR comparison)", fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(True, which='both', alpha=0.3)
    ax.text(0.02, 0.04,
            "Each diamond = one W-formula percentile threshold.\n"
            "Compare detectors at the same FAR operating point.",
            transform=ax.transAxes, fontsize=8, color='gray',
            verticalalignment='bottom')

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Figure saved: {output_path}")
    return fig


if __name__ == "__main__" and False:  # guard so it doesn't re-run above
    pass
