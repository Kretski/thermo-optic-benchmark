import numpy as np

# ---- Load real oscilloscope data ----
data1 = np.genfromtxt("/mnt/user-data/uploads/STEP-20DB-V0_4-8.csv", delimiter=",", skip_header=2)
x1 = data1[:,0]*1e6 + 0.8
y1 = data1[:,1]
y1n = (y1-np.min(y1))/(np.max(y1)-np.min(y1)) - 0.03

data2 = np.genfromtxt("/mnt/user-data/uploads/STEP-20DB-V0_4-IN8.csv", delimiter=",", skip_header=2)
x2 = data2[:,0]*1e6 + 553.14
y2 = data2[:,1]
# electrical input trace -> used as P_diss_W proxy for w_formula_detector

dt_us = np.mean(np.diff(x1))
dt_s = dt_us * 1e-6

# ---- Exact detector implementations from thermal_model_v3.py ----
def cusum_detector(residual, k_sigma, h_sigma, baseline_mu, baseline_sigma):
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

def windowed_glr_detector(residual, window, baseline_mu, baseline_sigma, threshold_sigma):
    N = len(residual)
    glr = np.zeros(N)
    alarms = np.zeros(N, dtype=bool)
    sigma2 = baseline_sigma ** 2
    h_glr = window * threshold_sigma**2 / 2.0
    for i in range(window - 1, N):
        w = residual[i - window + 1: i + 1]
        x_bar = w.mean()
        delta_hat = max(0.0, x_bar - baseline_mu)
        glr_val = window * delta_hat**2 / (2.0 * sigma2) if sigma2 > 0 else 0.0
        glr[i] = glr_val
        if glr_val > h_glr:
            alarms[i] = True
    return alarms, glr

def w_formula_detector(residual, P_diss_W, dt_s, baseline_sigma,
                       Q_weight=1.0, D_weight=1.0, T_weight=1.0):
    P_mean = P_diss_W.mean()
    if P_mean == 0:
        P_mean = 1.0
    Q = Q_weight * P_diss_W / P_mean
    D = D_weight * np.exp(-np.abs(residual) / (baseline_sigma + 1e-30))
    cumstress = np.cumsum(np.abs(residual) * dt_s)
    norm_factor = baseline_sigma * (len(residual) * dt_s)
    T_stress = T_weight * cumstress / max(norm_factor, 1e-30)
    return Q * D - T_stress

def apply_w_threshold(W, fault_onset_step, percentile=1.0):
    pre = W[:fault_onset_step]
    w_threshold = np.percentile(pre, percentile) if len(pre) else W.min()
    alarms = W < w_threshold
    return alarms, float(w_threshold)

# ---- Baseline stats (pre-trigger) ----
baseline_mask = x1 < -0.5
baseline = y1n[baseline_mask]
mu0, sigma0 = baseline.mean(), baseline.std()

true_onset_idx = np.searchsorted(x1, 0.0)

# ---- Fixed threshold (mu+5sigma), matching the paper's style label ----
thr = mu0 + 5*sigma0
fixed_idx = np.argmax(y1n > thr)
fixed_alarm = fixed_idx if y1n[fixed_idx] > thr else None

# ---- CUSUM, exact paper calibration k=2.0, h=6.0 ----
alarms_c, _ = cusum_detector(y1n, k_sigma=2.0, h_sigma=6.0, baseline_mu=mu0, baseline_sigma=sigma0)
cusum_idx = np.argmax(alarms_c) if alarms_c.any() else None

# ---- Windowed GLR, exact paper calibration window=20, threshold_sigma=3.5 ----
alarms_g, _ = windowed_glr_detector(y1n, window=20, baseline_mu=mu0, baseline_sigma=sigma0, threshold_sigma=3.5)
glr_idx = np.argmax(alarms_g) if alarms_g.any() else None

# ---- W-formula: need P_diss_W. Use the real electrical input (heater drive) ----
# Resample y2 onto x1's time grid (nearest), since the two channels aren't synchronized
# in time (per Qing Li), but represent the same physical event structurally.
# Use the electrical step itself (offset to positive) as a proxy dissipated-power trace.
y2_clean = np.nan_to_num(y2, nan=np.nanmin(y2))
P_diss_proxy = (y2_clean - np.nanmin(y2_clean))
P_diss_on_grid = np.interp(x1, x2, P_diss_proxy)  # align onto x1's grid, out-of-range clipped
if P_diss_on_grid.sum() == 0 or np.all(P_diss_on_grid == P_diss_on_grid[0]):
    print("WARNING: electrical channel did not align meaningfully onto x1 grid (unsynchronized triggers)")

W = w_formula_detector(y1n, P_diss_on_grid, dt_s, sigma0)
alarms_w, w_thr = apply_w_threshold(W, true_onset_idx, percentile=1.0)
w_idx = np.argmax(alarms_w[true_onset_idx:]) + true_onset_idx if alarms_w[true_onset_idx:].any() else None

print(f"True onset: index {true_onset_idx}, t=0.000 us\n")
print(f"{'Detector':<25}{'Alarm idx':>12}{'Alarm t (us)':>16}{'Delay (us)':>14}")
print("-"*67)
for name, idx in [("Fixed (mu+5sigma)", fixed_alarm), ("CUSUM (k=2.0,h=6.0)", cusum_idx),
                  ("Windowed GLR (w=20,3.5s)", glr_idx), ("W-formula (P=1%ile)", w_idx)]:
    if idx is None:
        print(f"{name:<25}{'NO ALARM':>12}")
    else:
        print(f"{name:<25}{idx:>12}{x1[idx]:>16.4f}{x1[idx]-x1[true_onset_idx]:>14.4f}")
