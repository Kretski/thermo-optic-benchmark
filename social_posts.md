# Постове за HN и LinkedIn

---

## Hacker News (Ask HN / Show HN)

**Title:**
Show HN: Benchmark of 7 sequential fault detectors on a silicon photonic microring model

**Body:**
I built a physically calibrated simulation framework comparing sequential change detection algorithms (CUSUM, EWMA, SPRT, Shiryaev-Roberts, Kalman, W-formula, fixed threshold) for thermal fault monitoring in silicon microring resonators.

The context: microring resonators are used in silicon photonic transceivers and are sensitive to thermal drift (~75 pm/°C). Existing systems use PID feedback and fixed thresholds. No prior work has systematically compared sequential change-point methods on this problem.

What the framework does:
- First-order thermal RC model calibrated to published device parameters (k₀ = 11.7 pm/mW, τ = 7 µs)
- Realistic colored noise: Gaussian + 1/f (pink) + random-walk drift
- 4 fault scenarios: heater efficiency loss, thermal time constant increase, thermal drift, contact degradation
- Detectors compared via Delay–FAR operating frontiers (threshold sweep, N=200 Monte Carlo trials each), not single operating points
- All results are simulation-based — no real device measurements yet

Key finding: no single detector dominates all fault scenarios. CUSUM is most consistent at near-zero FAR. Shiryaev-Roberts is fastest for step-type faults. W-formula shows advantage for impulsive degradation.

Code + figures: https://doi.org/10.5281/zenodo.21782606

Honest limitations: simulation only, SiC parameters (no Si SOI experimental validation yet). Looking for anyone with heater transient data who wants to collaborate.

---

## LinkedIn

**Post:**

Публикувах benchmark framework за сравнение на 7 алгоритъма за sequential change detection при термична деградация в silicon photonic microring резонатори.

Контекст: microring резонаторите се използват в оптични трансивери и са чувствителни към термичен дрейф. Съществуващите системи ползват PID feedback и fixed threshold. Колкото ми е известно, никой досега не е сравнявал систематично sequential detection методи (CUSUM, EWMA, SPRT, Shiryaev-Roberts, Kalman) за тази задача.

Основни характеристики на framework-а:
→ Физически калибриран термичен модел (параметри от публикувани измервания)
→ Реалистичен шум: Gaussian + 1/f + random-walk drift
→ 4 типа повреди: деградация на нагревател, термичен дрейф, контактна деградация
→ Сравнение чрез пълни Delay–FAR frontiers (не само единични operating points)
→ 200 Monte Carlo trials на конфигурация

Ключов извод: нито един детектор не доминира при всички типове повреди. CUSUM е най-стабилен при нулев FAR. Shiryaev-Roberts е най-бърз при step-type повреди.

Всичко е open source под MIT лиценз.

Кодът и фигурите: https://doi.org/10.5281/zenodo.21782606

Ако имате достъп до реални heater transient данни от silicon photonic устройства — радвам се да се свържем.

#SiliconPhotonics #FaultDetection #SignalProcessing #SequentialDetection #OpenScience

---

## Бележки

- HN: постирай в ~18:00 UTC в делничен ден (понеделник-сряда) за максимален трафик
- LinkedIn: сутринта е по-добре (8-10 ч. местно)
- Не казвай "first" или "novel" — просто описвай какво прави
- Ако някой попита за сравнение с ML методи — честният отговор: "not evaluated, would be interesting"
