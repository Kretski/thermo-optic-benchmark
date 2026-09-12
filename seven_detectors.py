import numpy as np

# ---- Load real data exactly as before ----
data1 = np.genfromtxt("/mnt/user-data/uploads/STEP-20DB-V0_4-8.csv", delimiter=",", skip_header=2)
x1 = data1[:,0] * 1e6 + 0.8   # time, microseconds, shifted to zero
y1 = data1[:,1]
y1n = (y1 - np.min(y1)) / (np.max(y1) - np.min(y1)) - 0.03

dt = np.mean(np.diff(x1))
print(f"Sample interval dt = {dt:.5f} us, N = {len(x1)} points, total span = {x1[-1]-x1[0]:.2f} us")

# ---- Baseline (pre-trigger, x<0) statistics for calibrating detectors ----
baseline_mask = x1 < -0.5   # comfortably before the step, avoid edge effects
baseline = y1n[baseline_mask]
mu0, sigma0 = np.mean(baseline), np.std(baseline)
print(f"Baseline: n={len(baseline)}, mean={mu0:.5f}, std={sigma0:.5f}")

signal = y1n
t = x1
n = len(signal)
true_onset_idx = np.searchsorted(t, 0.0)  # index where t crosses 0

def delay_from_onset(alarm_idx):
    if alarm_idx is None:
        return None
    return (t[alarm_idx] - t[true_onset_idx])

results = {}

# 1. Fixed threshold: 5-sigma above baseline mean
thr = mu0 + 5*sigma0
idx = np.argmax(signal > thr)
results['Fixed threshold (5σ)'] = idx if signal[idx] > thr else None

# 2. CUSUM
k = 0.5*sigma0   # allowance
h = 5*sigma0     # decision threshold
cusum = 0.0
alarm = None
for i in range(n):
    cusum = max(0, cusum + (signal[i]-mu0) - k)
    if cusum > h:
        alarm = i
        break
results['CUSUM'] = alarm

# 3. EWMA control chart
lam = 0.2
ewma = mu0
L = 3.0  # control limit multiplier
alarm = None
for i in range(n):
    ewma = lam*signal[i] + (1-lam)*ewma
    sigma_ewma = sigma0*np.sqrt(lam/(2-lam))
    if ewma > mu0 + L*sigma_ewma:
        alarm = i
        break
results['EWMA'] = alarm

# 4. SPRT (sequential probability ratio test) for mean shift mu0 -> mu0+delta
delta = 5*sigma0  # anticipated shift magnitude to detect
alpha_err, beta_err = 0.01, 0.01
A = np.log((1-beta_err)/alpha_err)
llr = 0.0
alarm = None
for i in range(n):
    # log-likelihood ratio increment for Gaussian mean shift, known sigma0
    llr += (delta/(sigma0**2))*(signal[i]-mu0-delta/2)
    if llr > A:
        alarm = i
        break
results['SPRT'] = alarm

# 5. Shiryaev-Roberts
sr = 0.0
sr_thr = 1/alpha_err
alarm = None
for i in range(n):
    lr = np.exp((delta/(sigma0**2))*(signal[i]-mu0-delta/2))
    sr = (1+sr)*lr
    if sr > sr_thr:
        alarm = i
        break
results['Shiryaev-Roberts'] = alarm

# 6. Kalman innovation
# simple constant-level Kalman filter, flag when |innovation| > 5 sigma of innovation std
q = 1e-6   # process noise
r = sigma0**2  # measurement noise
x_est, p_est = mu0, sigma0**2
alarm = None
for i in range(n):
    p_pred = p_est + q
    innov = signal[i] - x_est
    s_innov = p_pred + r
    if abs(innov) > 5*np.sqrt(s_innov):
        alarm = i
        break
    kgain = p_pred / s_innov
    x_est = x_est + kgain*innov
    p_est = (1-kgain)*p_pred
results['Kalman innovation'] = alarm

# 7. "W-formula": simple normalized deviation score W = 1 - |signal-mu0|/(k*sigma0), alarm when W<0
k_w = 5.0
alarm = None
for i in range(n):
    W = 1 - abs(signal[i]-mu0)/(k_w*sigma0)
    if W < 0:
        alarm = i
        break
results['W-formula (deviation score)'] = alarm

print(f"\nTrue onset at t=0.000 us (index {true_onset_idx})\n")
print(f"{'Detector':<28}{'Alarm index':>12}{'Alarm time (us)':>18}{'Delay (us)':>14}")
print("-"*72)
for name, idx in results.items():
    if idx is None:
        print(f"{name:<28}{'NO ALARM':>12}")
    else:
        delay = delay_from_onset(idx)
        print(f"{name:<28}{idx:>12}{t[idx]:>18.4f}{delay:>14.4f}")
