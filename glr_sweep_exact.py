import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict

# ---- Exact reconstruction from thermal_model_v3.py (all viewed sections) ----

@dataclass
class Params:
    k0: float = 72.9
    dk_dT: float = 0.0
    heater_eff: float = 100.0
    tau_th: float = 7e-6
    sigma_gauss: float = 2.0
    sigma_pink: float = 1.0
    sigma_drift: float = 0.5
    P_base_mW: float = 5.0

    def k_at_T(self, T, T_ref=25.0):
        return self.k0 + self.dk_dT*(T-T_ref)
    def R_th(self, T_C=25.0):
        return (self.heater_eff / self.k_at_T(T_C)) * 1000.0
    def C_th(self, T_C=25.0):
        return self.tau_th / self.R_th(T_C)

def generate_realistic_noise(N, dt_s, gaussian_sigma, pink_amplitude,
                              drift_sigma_per_sqrt_s, rng):
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
    return white + pink + np.cumsum(steps)

def cusum_detector(residual, k_sigma, h_sigma, baseline_mu, baseline_sigma):
    k = k_sigma * baseline_sigma
    h = h_sigma * baseline_sigma
    alarms = np.zeros(len(residual), dtype=bool)
    s = 0.0
    for i, x in enumerate(residual):
        s = max(0.0, s + (x - baseline_mu) - k)
        if s > h:
            alarms[i] = True
            s = 0.0
    return alarms

def windowed_glr_detector(residual, window, baseline_mu, baseline_sigma, threshold_sigma):
    N = len(residual)
    alarms = np.zeros(N, dtype=bool)
    sigma2 = baseline_sigma ** 2
    h_glr = window * threshold_sigma**2 / 2.0
    for i in range(window - 1, N):
        w = residual[i - window + 1: i + 1]
        delta_hat = max(0.0, w.mean() - baseline_mu)
        glr_val = window * delta_hat**2 / (2.0 * sigma2) if sigma2 > 0 else 0.0
        if glr_val > h_glr:
            alarms[i] = True
    return alarms

def build_simulation(params, fault_onset_s, fault_ramp_s, fault_mult,
                     total_time_s, dt_s, seed, use_realistic_noise=True):
    rng = np.random.default_rng(seed)
    T_amb = 25.0
    R_th_nominal = params.R_th(T_amb)
    C_th = params.C_th(T_amb)
    N = int(total_time_s / dt_s)
    fos = int(fault_onset_s / dt_s)
    ramp_len = max(1, int(fault_ramp_s / dt_s))

    P_diss_W = (params.P_base_mW*1e-3
               + 0.3e-3*np.sin(2*np.pi*50e3*np.arange(N)*dt_s)
               + rng.normal(0, 0.05e-3, N))

    R_th_actual = np.full(N, R_th_nominal)
    R_th_actual[fos:fos+ramp_len] = np.linspace(R_th_nominal, R_th_nominal*fault_mult, ramp_len)
    R_th_actual[fos+ramp_len:] = R_th_nominal*fault_mult

    T_actual = np.zeros(N); T_actual[0] = T_amb
    for i in range(1, N):
        dTdt = (P_diss_W[i-1] - (T_actual[i-1]-T_amb)/R_th_actual[i-1]) / C_th
        T_actual[i] = T_actual[i-1] + dTdt*dt_s
    k_series = np.array([params.k_at_T(T) for T in T_actual])

    if use_realistic_noise:
        noise = generate_realistic_noise(N, dt_s, params.sigma_gauss, params.sigma_pink,
                                         params.sigma_drift, rng)
    else:
        noise = rng.normal(0, params.sigma_gauss, N)

    delta_lambda_measured = k_series*(T_actual-T_amb) + noise

    T_model = np.zeros(N); T_model[0] = T_amb
    for i in range(1, N):
        dTdt = (P_diss_W[i-1] - (T_model[i-1]-T_amb)/R_th_nominal) / C_th
        T_model[i] = T_model[i-1] + dTdt*dt_s
    k_model_series = np.array([params.k_at_T(T) for T in T_model])
    T_implied = delta_lambda_measured / k_model_series + T_amb
    residual = T_implied - T_model

    pre = residual[:fos]
    return {"residual": residual, "mu_pre": pre.mean(), "sigma_pre": pre.std(),
            "fault_onset_step": fos, "dt_s": dt_s, "N": N}

def run_one(params, mult, seed, glr_window, glr_thr, cusum_k=2.0, cusum_h=6.0,
           fault_onset_s=200e-6, fault_ramp_s=20e-6, total_time_s=500e-6, dt_s=None):
    if dt_s is None:
        dt_s = params.tau_th/100.0
    sim = build_simulation(params, fault_onset_s, fault_ramp_s, mult, total_time_s, dt_s, seed)
    residual, mu_pre, sigma_pre, fos = sim["residual"], sim["mu_pre"], sim["sigma_pre"], sim["fault_onset_step"]

    def delay_fa(alarms):
        post = np.where(alarms[fos:])[0]
        delay = post[0]*sim["dt_s"]*1e6 if len(post) else None
        fa = int(alarms[:fos].sum())
        return delay, fa

    cusum_alarms = cusum_detector(residual, cusum_k, cusum_h, mu_pre, sigma_pre)
    glr_alarms = windowed_glr_detector(residual, glr_window, mu_pre, sigma_pre, glr_thr)
    return delay_fa(cusum_alarms), delay_fa(glr_alarms), fos

params = Params()
N_TRIALS = 200
mult = 2.2  # matching the paper's primary severity for Table II

print(f"Exact-methodology sweep, N={N_TRIALS} trials, fault severity x{mult}\n")
print("CUSUM baseline (k=2.0, h=6.0):")
delays, fas, fos_ref = [], [], None
for seed in range(N_TRIALS):
    (cd, cfa), _, fos = run_one(params, mult, seed, glr_window=20, glr_thr=3.5)
    if cd is not None: delays.append(cd)
    fas.append(cfa)
    fos_ref = fos
print(f"  Pd={len(delays)/N_TRIALS:.3f}  delay={np.mean(delays):.3f}us  FAR={np.mean(fas)/fos_ref:.2e}\n")

print("GLR window/threshold sweep:")
print(f"{'window':>8}{'thr':>6}{'Pd':>7}{'delay(us)':>11}{'FAR':>12}")
results = []
for window in [5, 8, 10, 15, 20]:
    for thr in [2.0, 2.5, 3.0, 3.5]:
        delays, fas = [], []
        for seed in range(N_TRIALS):
            _, (gd, gfa), fos = run_one(params, mult, seed, glr_window=window, glr_thr=thr)
            if gd is not None: delays.append(gd)
            fas.append(gfa)
        pd = len(delays)/N_TRIALS
        md = np.mean(delays) if delays else float('nan')
        far = np.mean(fas)/fos
        results.append((window, thr, pd, md, far))
        print(f"{window:>8}{thr:>6.1f}{pd:>7.2f}{md:>11.3f}{far:>12.2e}")

print(f"\nFor reference, paper's Table II (fault x2.2): CUSUM h=6 -> 2.09us @ FAR~0; "
      f"CUSUM h=2 -> 1.79us @ FAR=9.6e-5")
