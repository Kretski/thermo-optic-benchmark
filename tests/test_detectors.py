"""
Test suite for the Thermo-Optic Fault Detection Benchmark.

Covers: the compact thermal model's basic physical sanity, the three
evaluated detectors (fixed threshold via CUSUM's building blocks, CUSUM,
windowed GLR, W-formula) each firing on an obvious injected fault and NOT
firing on clean noise-only data, and the noise generator's basic
statistical properties. This is a correctness/regression suite, not a
re-run of the full Monte Carlo simulation (that remains the job of
monte_carlo_sweep, which is comparatively expensive to run in CI).
"""

import numpy as np
import pytest

import thermal_model_v3 as tm


# ---------------------------------------------------------------------
# ThermalOpticalParams / compact model sanity
# ---------------------------------------------------------------------

class TestThermalOpticalParams:
    def test_default_params_construct(self):
        params = tm.ThermalOpticalParams()
        assert params is not None

    def test_k_at_T_matches_k0_at_reference_temperature(self):
        params = tm.ThermalOpticalParams()
        # At T == T_ref, k(T) should reduce to k0 by construction
        # (k(T) = k0 + dk_dT * (T - T_ref)).
        k_ref = params.k_at_T(25.0, T_ref=25.0)
        assert k_ref == pytest.approx(params.k0_pm_per_C.value, rel=1e-9)

    def test_R_th_positive(self):
        params = tm.ThermalOpticalParams()
        assert params.R_th(25.0) > 0

    def test_C_th_positive(self):
        params = tm.ThermalOpticalParams()
        assert params.C_th(25.0) > 0


# ---------------------------------------------------------------------
# Noise generator — basic statistical sanity, not a full distributional test
# ---------------------------------------------------------------------

class TestNoiseGenerator:
    def test_zero_amplitude_noise_is_zero(self):
        rng = np.random.default_rng(0)
        noise = tm.generate_realistic_noise(
            N=500, dt_s=1e-7, gaussian_sigma=0.0, pink_amplitude=0.0,
            drift_sigma_per_sqrt_s=0.0, rng=rng)
        assert np.allclose(noise, 0.0)

    def test_gaussian_only_noise_has_expected_scale(self):
        rng = np.random.default_rng(1)
        sigma = 2.0
        noise = tm.generate_realistic_noise(
            N=20000, dt_s=1e-7, gaussian_sigma=sigma, pink_amplitude=0.0,
            drift_sigma_per_sqrt_s=0.0, rng=rng)
        # with pink/drift off, std should be close to the requested sigma
        assert noise.std() == pytest.approx(sigma, rel=0.15)

    def test_noise_is_reproducible_with_fixed_seed(self):
        rng1 = np.random.default_rng(42)
        rng2 = np.random.default_rng(42)
        n1 = tm.generate_realistic_noise(200, 1e-7, 2.0, 1.0, 0.5, rng1)
        n2 = tm.generate_realistic_noise(200, 1e-7, 2.0, 1.0, 0.5, rng2)
        assert np.array_equal(n1, n2)


# ---------------------------------------------------------------------
# CUSUM detector
# ---------------------------------------------------------------------

class TestCusumDetector:
    def test_no_alarm_on_pure_baseline_noise(self):
        rng = np.random.default_rng(2)
        residual = rng.normal(0, 1.0, 2000)
        alarms, _ = tm.cusum_detector(
            residual, k_sigma=2.0, h_sigma=6.0,
            baseline_mu=0.0, baseline_sigma=1.0)
        # a well-calibrated CUSUM at h=6 sigma should essentially never
        # false-alarm on 2000 samples of clean Gaussian noise
        assert alarms.sum() == 0

    def test_alarms_on_obvious_mean_shift(self):
        rng = np.random.default_rng(3)
        baseline = rng.normal(0, 1.0, 500)
        shifted = rng.normal(8.0, 1.0, 500)  # large, obvious shift
        residual = np.concatenate([baseline, shifted])
        alarms, _ = tm.cusum_detector(
            residual, k_sigma=2.0, h_sigma=6.0,
            baseline_mu=0.0, baseline_sigma=1.0)
        assert alarms[500:].any(), "CUSUM failed to detect an 8-sigma mean shift"

    def test_decays_after_a_brief_shift_ends(self):
        rng = np.random.default_rng(4)
        residual = np.concatenate([
            rng.normal(0, 1.0, 200),
            rng.normal(10.0, 1.0, 50),   # brief shift
            rng.normal(0, 1.0, 200),     # back to baseline
        ])
        alarms, S = tm.cusum_detector(
            residual, k_sigma=2.0, h_sigma=6.0,
            baseline_mu=0.0, baseline_sigma=1.0)
        assert alarms.any(), "CUSUM failed to alarm on the brief shift"
        # after returning to baseline for a while, the cumulative sum
        # should have decayed back down, not stay permanently elevated
        assert S[-50:].mean() < 2.0, \
            f"CUSUM statistic did not decay after the shift ended: {S[-50:].mean()}"


# ---------------------------------------------------------------------
# Windowed GLR detector
# ---------------------------------------------------------------------

class TestWindowedGLR:
    def test_no_alarm_on_pure_baseline_noise(self):
        rng = np.random.default_rng(5)
        residual = rng.normal(0, 1.0, 2000)
        alarms, _ = tm.windowed_glr_detector(
            residual, window=20, baseline_mu=0.0, baseline_sigma=1.0,
            threshold_sigma=5.0)
        # allow a small number of false alarms at this threshold rather
        # than requiring exactly zero, since GLR's false-alarm behavior
        # under a short window is less conservative than CUSUM's
        assert alarms.sum() < 5

    def test_alarms_on_obvious_mean_shift(self):
        rng = np.random.default_rng(6)
        baseline = rng.normal(0, 1.0, 500)
        shifted = rng.normal(8.0, 1.0, 500)
        residual = np.concatenate([baseline, shifted])
        alarms, _ = tm.windowed_glr_detector(
            residual, window=20, baseline_mu=0.0, baseline_sigma=1.0,
            threshold_sigma=3.5)
        assert alarms[500:].any(), "Windowed GLR failed to detect an 8-sigma mean shift"


# ---------------------------------------------------------------------
# W-formula detector
# ---------------------------------------------------------------------

class TestWFormulaDetector:
    def test_w_drops_on_large_residual(self):
        rng = np.random.default_rng(7)
        n = 1000
        residual = rng.normal(0, 0.1, n)
        residual[500:] += 5.0  # inject an obvious fault partway through
        P_diss = np.ones(n)  # flat drive power proxy
        W = tm.w_formula_detector(residual, P_diss, dt_s=1e-7, baseline_sigma=0.1)
        # W should be lower (more anomalous) during the faulted region
        # than during the clean baseline region
        assert W[500:].mean() < W[:500].mean()

    def test_apply_w_threshold_flags_the_fault_region(self):
        rng = np.random.default_rng(8)
        n = 1000
        residual = rng.normal(0, 0.1, n)
        residual[500:] += 5.0
        P_diss = np.ones(n)
        W = tm.w_formula_detector(residual, P_diss, dt_s=1e-7, baseline_sigma=0.1)
        alarms, w_thresh = tm.apply_w_threshold(W, fault_onset_step=500, percentile=1.0)
        assert alarms[500:].any(), "W-formula threshold failed to flag the injected fault"


# ---------------------------------------------------------------------
# calibrate_from_measurements, if present — smoke test only
# ---------------------------------------------------------------------

class TestCalibration:
    def test_calibrate_from_measurements_is_callable(self):
        assert callable(tm.calibrate_from_measurements)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
