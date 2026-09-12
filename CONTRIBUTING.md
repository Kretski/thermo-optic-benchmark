# Contributing

Thank you for considering a contribution to this benchmark.

## Reporting issues

Please open a GitHub issue describing:
- What you expected to happen
- What actually happened
- Python version and OS
- A minimal script reproducing the problem, if applicable

## Proposing changes

1. Fork the repository and create a branch from `main`.
2. Make your change. If it touches detector logic (`cusum_detector`,
   `windowed_glr_detector`, `w_formula_detector`) or the thermal model
   (`ThermalOpticalParams`), add or update a test in `tests/test_detectors.py`
   covering the change.
3. Run the test suite locally before opening a pull request:
   ```bash
   pip install pytest
   pytest tests/ -v
   ```
   All tests must pass.
4. Open a pull request describing what changed and why. Reference any
   related issue.

## Scope of contributions we're especially interested in

- Additional detectors (BOCPD, Shiryaev-Roberts, EWMA, Kalman innovation
  monitors — see the manuscript's Discussion section for why these were
  not included in the original comparison)
- Full recursive GLR (Willsky-Jones), as a proper counterpart to the
  currently-implemented windowed approximation
- Hardware-in-the-loop validation against additional real transients,
  beyond the single real measurement in the Experimental Validation
  section
- Bug reports on the noise model, calibration routine, or Monte Carlo
  sweep

## Code style

Follow the existing style in `thermal_model_v3.py`: type hints on
function signatures, NumPy-style docstrings, and explicit unit
annotations in variable names and comments where physical quantities
are involved (this codebase has previously had bugs caused by unit
mismatches, so being explicit here is not optional style preference).

## Questions

Open an issue, or see the contact details in the main README.
