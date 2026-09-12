"""
Extended detector suite for comparison against W-formula.
All detectors operate on the same residual r(t) and pre-fault baseline (mu, sigma).

Detectors implemented:
  1. EWMA           — exponentially weighted moving average
  2. SPRT           — sequential probability ratio test (one-sided)
  3. Shiryaev-Roberts — SR statistic (CUSUM variant, different reset)
  4. Recursive GLR  — exact recursive GLRT for Gaussian shift (Willsky-Jones)
  5. Kalman innovation — model-based, treats residual as innovation sequence
  6. BOCPD          — Bayesian online changepoint detection (Adams & MacKay 2007),
                      Gaussian likelihood with known variance

Each function signature:
    detector(residual, baseline_mu, baseline_sigma, **params)
    -> (alarm_bool_array, statistic_trace)

All are one-sided (upward shift detection) to match the fault scenario.
"""

import numpy as np
from typing import Tuple


# ---------------------------------------------------------------------------
# 1. EWMA
# ---------------------------------------------------------------------------

def ewma_detector(residual: np.ndarray, baseline_mu: float, baseline_sigma: float,
                  lam: float = 0.2, threshold_sigma: float = 3.0
                  ) -> Tuple[np.ndarray, np.ndarray]:
    """
    EWMA control chart.
    Z_i = (1-lam)*Z_{i-1} + lam*r_i
    Alarm when Z_i > mu + threshold_sigma * sigma * sqrt(lam/(2-lam))
    One-sided upper (fault = increase in residual).
    lam in (0,1]: smaller = more smoothing = slower but more stable.
    """
    N = len(residual)
    Z = np.zeros(N)
    z = baseline_mu
    control_limit = baseline_mu + threshold_sigma * baseline_sigma * np.sqrt(lam / (2 - lam))
    alarms = np.zeros(N, dtype=bool)
    trace = np.zeros(N)
    for i, r in enumerate(residual):
        z = (1 - lam) * z + lam * r
        Z[i] = z
        trace[i] = z
        if z > control_limit:
            alarms[i] = True
    return alarms, trace


# ---------------------------------------------------------------------------
# 2. SPRT (one-sided)
# ---------------------------------------------------------------------------

def sprt_detector(residual: np.ndarray, baseline_mu: float, baseline_sigma: float,
                  delta: float = 1.0, alpha: float = 0.001, beta: float = 0.001
                  ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Sequential Probability Ratio Test (Wald 1947), one-sided.
    H0: mean = mu_0 = baseline_mu
    H1: mean = mu_1 = baseline_mu + delta * baseline_sigma

    Log-likelihood ratio:
        logLR_i += (r_i - mu_0) * delta/sigma - delta^2/2

    Alarm when logLR > log((1-beta)/alpha)  [reject H0]
    Reset after alarm.

    alpha: false alarm probability (Type I)
    beta:  miss probability (Type II)
    delta: expected shift in units of sigma
    """
    mu1 = baseline_mu + delta * baseline_sigma
    sigma = baseline_sigma
    A = np.log((1 - beta) / alpha)   # upper threshold (alarm)
    B = np.log(beta / (1 - alpha))   # lower threshold (accept H0) — not used in sequential

    N = len(residual)
    logLR = np.zeros(N)
    alarms = np.zeros(N, dtype=bool)
    s = 0.0
    for i, r in enumerate(residual):
        # log f(r|H1) - log f(r|H0) for Gaussian
        s += ((r - baseline_mu) * (mu1 - baseline_mu) - 0.5 * (mu1 - baseline_mu)**2) / sigma**2
        logLR[i] = s
        if s >= A:
            alarms[i] = True
            s = 0.0  # reset
        elif s <= B:
            s = 0.0  # reset (accept H0 momentarily)
    return alarms, logLR


# ---------------------------------------------------------------------------
# 3. Shiryaev-Roberts
# ---------------------------------------------------------------------------

def shiryaev_roberts_detector(residual: np.ndarray, baseline_mu: float,
                               baseline_sigma: float,
                               delta: float = 1.0, threshold: float = 500.0
                               ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Shiryaev-Roberts statistic (Roberts 1966, Pollak 1985).
    R_i = (R_{i-1} + 1) * L_i
    where L_i = f(r_i|H1) / f(r_i|H0) is the likelihood ratio.

    Alarm when R_i > threshold.
    Unlike CUSUM, no max(0,...) — accumulates evidence for ALL possible
    changepoint times, giving it minimax-optimal average run length properties
    under certain conditions.

    threshold: scalar alarm level (not in sigma units — SR grows exponentially).
    delta: expected shift in sigma units.
    """
    mu1 = baseline_mu + delta * baseline_sigma
    sigma = baseline_sigma
    N = len(residual)
    R = np.zeros(N)
    alarms = np.zeros(N, dtype=bool)
    r_stat = 0.0
    for i, x in enumerate(residual):
        # Gaussian LR: exp((x-mu0)*(mu1-mu0)/sigma^2 - (mu1-mu0)^2/(2*sigma^2))
        log_lr = ((x - baseline_mu) * (mu1 - baseline_mu)
                  - 0.5 * (mu1 - baseline_mu)**2) / sigma**2
        lr = np.exp(np.clip(log_lr, -50, 50))
        r_stat = (r_stat + 1.0) * lr
        R[i] = r_stat
        if r_stat > threshold:
            alarms[i] = True
            r_stat = 0.0  # reset
    return alarms, R


# ---------------------------------------------------------------------------
# 4. Recursive GLR (Willsky-Jones 1976)
# ---------------------------------------------------------------------------

def recursive_glr_detector(residual: np.ndarray, baseline_mu: float,
                            baseline_sigma: float,
                            threshold_sigma: float = 5.0,
                            min_segment: int = 5
                            ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Recursive (exact) Generalized Likelihood Ratio test.
    At each time i, tests all possible changepoint times k in [i-max_lag, i-min_segment].
    GLR_i = max_k { n_k * delta_hat_k^2 / (2*sigma^2) }
    where n_k = i - k + 1, delta_hat_k = mean(r[k:i+1]) - mu_0

    Alarm when GLR_i > h_glr.

    This is the exact recursive GLR (no windowing approximation), but has
    O(n^2) complexity in the worst case — bounded by max_lag.

    threshold_sigma: alarm threshold calibrated as (threshold_sigma)^2 * n / 2
    min_segment: minimum segment length before testing (avoids trivial alarms)
    max_lag: how far back to test for changepoint (controls complexity)
    """
    sigma2 = baseline_sigma**2
    N = len(residual)
    glr = np.zeros(N)
    alarms = np.zeros(N, dtype=bool)
    max_lag = min(N, 200)  # cap for speed; adjust if needed

    # Running cumulative sum for efficient mean computation
    cumsum = np.zeros(N + 1)
    for i in range(N):
        cumsum[i + 1] = cumsum[i] + residual[i]

    for i in range(min_segment, N):
        best_glr = 0.0
        # Test all candidate changepoints from min_segment samples ago
        k_start = max(0, i - max_lag)
        for k in range(k_start, i - min_segment + 1):
            n_k = i - k + 1
            seg_mean = (cumsum[i + 1] - cumsum[k]) / n_k
            delta_hat = max(0.0, seg_mean - baseline_mu)  # one-sided
            glr_val = n_k * delta_hat**2 / (2.0 * sigma2) if sigma2 > 0 else 0.0
            if glr_val > best_glr:
                best_glr = glr_val
        glr[i] = best_glr
        h_glr = min_segment * threshold_sigma**2 / 2.0
        if best_glr > h_glr:
            alarms[i] = True
    return alarms, glr


# ---------------------------------------------------------------------------
# 5. Kalman innovation detector
# ---------------------------------------------------------------------------

def kalman_innovation_detector(residual: np.ndarray, baseline_mu: float,
                                baseline_sigma: float,
                                process_noise_ratio: float = 0.01,
                                threshold_sigma: float = 5.0
                                ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Kalman filter on the residual signal, treating it as a noisy measurement
    of a slowly-varying state (thermal drift).

    State model:  x_k = x_{k-1} + w_k,  w_k ~ N(0, Q)
    Measurement:  z_k = x_k + v_k,       v_k ~ N(0, R)

    Innovation = z_k - x_hat_{k|k-1}
    Alarm when |innovation| > threshold_sigma * sqrt(S_k)
    where S_k is the innovation variance.

    process_noise_ratio: Q/R — controls how much the state is allowed to drift.
    Small value = Kalman assumes slow drift → innovations are large during fault.
    """
    R_noise = baseline_sigma**2
    Q_noise = process_noise_ratio * R_noise

    N = len(residual)
    innovations = np.zeros(N)
    alarms = np.zeros(N, dtype=bool)

    x_hat = baseline_mu   # state estimate
    P = R_noise            # initial state covariance

    for i, z in enumerate(residual):
        # Predict
        x_pred = x_hat
        P_pred = P + Q_noise

        # Innovation
        innov = z - x_pred
        S = P_pred + R_noise   # innovation covariance

        # Update
        K = P_pred / S         # Kalman gain
        x_hat = x_pred + K * innov
        P = (1 - K) * P_pred

        innovations[i] = innov
        # One-sided alarm: upward innovation only
        if innov > threshold_sigma * np.sqrt(S):
            alarms[i] = True

    return alarms, innovations


# ---------------------------------------------------------------------------
# 6. BOCPD (Adams & MacKay 2007)
# ---------------------------------------------------------------------------

def bocpd_detector(residual: np.ndarray, baseline_mu: float,
                   baseline_sigma: float,
                   hazard: float = 1e-3,
                   threshold_prob: float = 0.5
                   ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Bayesian Online Changepoint Detection (Adams & MacKay 2007).
    Uses Gaussian likelihood with known variance (sigma from baseline).
    Prior on run length: geometric with hazard rate h (= 1/expected_run_length).

    At each step, computes P(r_t | data) — the distribution over current
    run lengths. An alarm is raised when P(run_length = 0 | data) > threshold_prob,
    i.e., when the posterior probability of a changepoint AT THIS STEP is high.

    hazard: 1/expected_run_length between changepoints (1e-3 = expect ~1000 steps)
    threshold_prob: alarm threshold on P(r_t=0 | x_{1:t})

    Note: This is the simplified version with Gaussian conjugate prior (known sigma).
    The full Adams-MacKay uses Student-t predictive with unknown mean and variance.
    """
    N = len(residual)
    sigma2 = baseline_sigma**2

    # Run-length distribution: R[k] = P(run_length = k at time t)
    # Initialise with run_length = 0 (start fresh)
    max_rl = N + 1
    R = np.zeros(max_rl)
    R[0] = 1.0   # at t=0, run_length = 0 with probability 1

    # Sufficient statistics per run-length hypothesis
    # For Gaussian with known variance: predictive mean = baseline_mu (prior),
    # updated by online mean of segment
    # We track: count, sum for each possible run-length
    counts = np.zeros(max_rl, dtype=int)
    sums   = np.zeros(max_rl)

    changepoint_probs = np.zeros(N)
    alarms = np.zeros(N, dtype=bool)

    for t in range(N):
        x = residual[t]

        # Predictive probabilities for each run-length hypothesis
        # Under run-length k: mean = (baseline_mu + sum_{last k}) / (1 + k) approx baseline_mu
        # Simplified: use fixed baseline_mu as prior mean for all hypotheses
        # Predictive: N(x | mu_k, sigma^2) where mu_k = running mean of segment
        pred = np.zeros(max_rl)
        for k in range(t + 1):
            if R[k] > 1e-300:
                if counts[k] == 0:
                    mu_k = baseline_mu
                else:
                    mu_k = sums[k] / counts[k]
                pred[k] = np.exp(-0.5 * (x - mu_k)**2 / sigma2) / np.sqrt(2 * np.pi * sigma2)

        # Growth: run_length k -> k+1 with probability (1-h)
        R_new = np.zeros(max_rl)
        R_new[1:t+2] = R[:t+1] * pred[:t+1] * (1 - hazard)

        # Changepoint: all run-lengths -> 0 with probability h
        R_new[0] = np.sum(R[:t+1] * pred[:t+1]) * hazard

        # Normalise
        norm = R_new[:t+2].sum()
        if norm > 0:
            R_new[:t+2] /= norm

        # Update sufficient statistics
        new_counts = np.zeros(max_rl, dtype=int)
        new_sums   = np.zeros(max_rl)
        new_counts[1:t+2] = counts[:t+1] + 1
        new_sums[1:t+2]   = sums[:t+1] + x
        new_counts[0] = 0
        new_sums[0]   = 0.0

        R = R_new
        counts = new_counts
        sums   = new_sums

        cp_prob = R[0]   # P(changepoint NOW)
        changepoint_probs[t] = cp_prob
        if cp_prob > threshold_prob:
            alarms[t] = True

    return alarms, changepoint_probs


# ---------------------------------------------------------------------------
# Sweep helper: one simulation, all detectors
# ---------------------------------------------------------------------------

def run_all_detectors(residual, P_diss_W, mu_pre, sigma_pre,
                      fault_onset_step, dt_s,
                      # CUSUM
                      cusum_k=2.0, cusum_h=6.0,
                      # EWMA
                      ewma_lam=0.2, ewma_thr=3.5,
                      # SPRT
                      sprt_delta=1.0, sprt_alpha=0.001, sprt_beta=0.001,
                      # SR
                      sr_delta=1.0, sr_thr=500.0,
                      # Recursive GLR
                      rglr_thr=3.5, rglr_min_seg=5,
                      # Kalman
                      kal_q=0.01, kal_thr=5.0,
                      # BOCPD
                      bocpd_h=1e-3, bocpd_thr=0.5,
                      # W-formula (imported from thermal_model_v3)
                      w_pct=0.01):
    import sys
    sys.path.insert(0, '/home/claude')
    from thermal_model_v3 import cusum_detector, w_formula_detector, apply_w_threshold

    fos = fault_onset_step

    def _delay_fa(alarms):
        post = np.where(alarms[fos:])[0]
        delay = post[0] * dt_s * 1e6 if len(post) > 0 else None
        fa = int(alarms[:fos].sum()) / max(fos, 1)
        return delay, fa

    results = {}

    # Fixed µ+5σ
    alarms = residual > (mu_pre + 5 * sigma_pre)
    results['Fixed µ+5σ'] = _delay_fa(alarms)

    # CUSUM k=2,h=6
    c_alarms, _ = cusum_detector(residual, cusum_k, cusum_h, mu_pre, sigma_pre)
    results['CUSUM k=2,h=6'] = _delay_fa(c_alarms)

    # W-formula 0.01%ile
    W = w_formula_detector(residual, P_diss_W, dt_s, sigma_pre)
    w_alarms, _ = apply_w_threshold(W, fos, percentile=w_pct)
    results['W-formula'] = _delay_fa(w_alarms)

    # EWMA
    e_alarms, _ = ewma_detector(residual, mu_pre, sigma_pre, ewma_lam, ewma_thr)
    results['EWMA'] = _delay_fa(e_alarms)

    # SPRT
    sp_alarms, _ = sprt_detector(residual, mu_pre, sigma_pre,
                                  sprt_delta, sprt_alpha, sprt_beta)
    results['SPRT'] = _delay_fa(sp_alarms)

    # Shiryaev-Roberts
    sr_alarms, _ = shiryaev_roberts_detector(residual, mu_pre, sigma_pre,
                                              sr_delta, sr_thr)
    results['Shiryaev-Roberts'] = _delay_fa(sr_alarms)

    # Recursive GLR
    rg_alarms, _ = recursive_glr_detector(residual, mu_pre, sigma_pre,
                                           rglr_thr, rglr_min_seg)
    results['Recursive GLR'] = _delay_fa(rg_alarms)

    # Kalman innovation
    k_alarms, _ = kalman_innovation_detector(residual, mu_pre, sigma_pre,
                                              kal_q, kal_thr)
    results['Kalman'] = _delay_fa(k_alarms)

    # BOCPD
    b_alarms, _ = bocpd_detector(residual, mu_pre, sigma_pre,
                                  bocpd_h, bocpd_thr)
    results['BOCPD'] = _delay_fa(b_alarms)

    return results
