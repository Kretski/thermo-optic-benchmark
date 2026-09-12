# Installation and Quickstart

This document covers installing and running the code. For the benchmark
methodology, results, and dataset description, see the main manuscript
and the Zenodo deposit (DOI: 10.5281/zenodo.21782606).

## Requirements

- Python 3.9+
- numpy, scipy, matplotlib (see `requirements.txt`)
- pytest (only needed to run the test suite)

## Install

```bash
git clone https://github.com/kretski/thermo-optic-benchmark.git
cd thermo-optic-benchmark
pip install -r requirements.txt
```

No compiled extensions; a plain `pip install` is sufficient on Linux,
macOS, and Windows.

## Running the test suite

```bash
pip install pytest
pytest tests/ -v
```

Expected: all tests pass in well under a second (the test suite checks
detector correctness on small synthetic signals, not the full Monte
Carlo sweep).

## Running the benchmark

```bash
python thermal_model_v3.py
```

This runs the full Monte Carlo sweep (Sections III-V of the manuscript)
and reproduces the operating-frontier and robustness-sweep tables. This
takes several minutes on a typical laptop (N=200 trials per configuration
across the swept threshold ranges).

## Reproducing the experimental validation section

The real oscilloscope traces used in Section VI of the manuscript
(`STEP-20DB-V0.4-8.csv`, `STEP-20DB-V0.4-IN8.csv`) are included in the
Zenodo deposit's supplementary data, with permission and attribution to
Q. Li (Carnegie Mellon University). See `run_real_detectors.py` for the
script that reproduces the fitted time constant and detection-delay
table from Section VI using these files.

## Common issues

**`ModuleNotFoundError` for numpy/scipy**: run
`pip install -r requirements.txt` in the same Python environment you're
using to run the scripts (a virtualenv mismatch is the most common cause).

**Monte Carlo sweep is slow**: this is expected — the sweep runs
hundreds of trials per configuration across multiple detector threshold
values. Reduce `N` in `monte_carlo_sweep` calls for a faster, lower-
precision check while developing, and restore the full `N` before
citing any results.

## Getting help

Open a GitHub issue (see `CONTRIBUTING.md`) or use the contact details
in the main README.
