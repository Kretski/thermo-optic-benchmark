# Cover Letter — Optics Express Submission

---

**To:** The Editors, Optics Express  
**From:** Dimitar Kretski, Center for Hydro- and Aerodynamics (CHA), Varna, Bulgaria  
**ORCID:** 0000-0001-5108-2243  
**Date:** August 2026

---

Dear Editors,

We submit for consideration the manuscript **"A Compact Thermo-Optic Health Monitoring Framework for Silicon Microring Resonators: Benchmark of Sequential Change Detection Methods"**.

**What the paper does.**
We develop a physically calibrated compact thermo-optic model for silicon microring resonators and use it to benchmark seven sequential change detection algorithms — fixed threshold, CUSUM, EWMA, SPRT, Shiryaev-Roberts, Kalman innovation, and a deterministic W-formula safety layer — by sweeping each detector's decision threshold to obtain full Delay–False Alarm Rate (FAR) operating frontiers. The key methodological contribution is this frontier-based evaluation: rather than comparing detectors at a single operating point, we characterise the entire achievable trade-off space for each method under a realistic colored-noise model (Gaussian + 1/f + random-walk drift) across four synthetic fault scenarios.

The proposed framework is intended as a non-invasive software layer operating alongside existing PID-based thermal controllers, monitoring the residual between a physical model prediction and the measured wavelength shift, without modifying the hardware or the control loop.

**Principal findings.**
(1) CUSUM consistently improves upon the fixed-threshold detector across the evaluated operating range, with 7–17% shorter detection delay at matched FAR. (2) Shiryaev-Roberts achieves the shortest detection delay for step-type thermal faults at Pd ≥ 0.95. (3) The W-formula detector shows a relative advantage for impulsive contact degradation scenarios. (4) No single detector dominates all fault scenarios — the relative ranking varies with fault type. (5) To the best of our knowledge, this is the first systematic benchmark of sequential change detection methods for thermo-optic fault monitoring in microring resonators.

**Scope fit.**
Optics Express regularly publishes work on silicon photonic device characterisation, thermal effects in integrated photonics, and photonic system reliability. This manuscript addresses the intersection of thermo-optic modelling and health monitoring — relevant to both researchers designing resonator-based systems and engineers responsible for their long-term operation. Existing photonic systems rely predominantly on PID feedback and fixed-threshold alarms; this work demonstrates that sequential statistical methods offer systematically characterised alternatives.

**Limitations stated.**
All results are simulation-based, using literature-calibrated parameters (k₀ = 72.9 pm/°C for Si SOI; k₀ = 11.7 pm/mW, τ = 7 µs for 4H-SiC from Sun et al., arXiv:2506.15035). Experimental validation on a fabricated device is explicitly identified as the essential next step in the paper.

**No conflicts of interest.** This work has not been submitted elsewhere and is not under review at another journal.

The simulation code, detector implementations, and figure reproduction scripts are openly available at:
https://doi.org/10.5281/zenodo.21782606 (MIT License)

We believe this work will be of interest to the Optics Express readership and welcome the editors' and reviewers' assessment.

Sincerely,

Dimitar Kretski  
Center for Hydro- and Aerodynamics (CHA)  
Varna, Bulgaria  
ORCID: 0000-0001-5108-2243  
kretski1@gmail.com  
https://doi.org/10.5281/zenodo.21782606

---

*Note: If submitting to IEEE Photonics Journal instead, replace the scope paragraph with:*
*"IEEE Photonics Journal publishes work on photonic devices, materials, and systems including reliability and monitoring applications. This manuscript addresses thermo-optic health monitoring in silicon microring resonators — a topic directly relevant to the journal's scope on integrated photonic device characterisation and performance monitoring. The benchmark framework and open-source code may also be of interest to the embedded photonic systems community."*
