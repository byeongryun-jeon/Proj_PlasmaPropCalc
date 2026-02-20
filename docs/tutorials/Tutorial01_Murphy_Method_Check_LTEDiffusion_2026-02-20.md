# Tutorial 01: Murphy Method Conformance Check (via LTEDiffusion)

Date: 2026-02-20  
Scope: Check whether current Tutorial 01 transport workflow is a true Murphy-style implementation.

## 1. LTEDiffusion basis used

Notebook:
- `LTEDiffusion` (`43760a8e-f80a-4800-a105-7bd391e2de73`)

Key sources loaded in the notebook:
- `Murphy_&_Tam_2014_JPhysD.pdf`
- `Murphy_&_Arundell_1994_PlasmaChemPlasmaProcess.pdf`
- `devoto1966.pdf`
- `devoto1967.pdf`
- `mason1967.pdf`
- `aziz1990.pdf`
- `aubreton1986.pdf`
- `barker1964.pdf`
- `1977AuJPh__30___61M.pdf` (Milloy)
- `frost1964.pdf`

Relevant LTEDiffusion conversation IDs used in this check:
- `91a87bc4-7a36-4bd4-b20f-ba1bffdcd0e9`
- `9e8bb74e-fcf6-4fbd-8902-03d0b50bbe2f`

## 2. Murphy-side reference rules (from LTEDiffusion)

LTEDiffusion summarizes Murphy 2014 as using the Murphy 1994 transport framework.

Transport approximation orders:
- Viscosity `mu`: 1st-order approximation.
- Heavy translational thermal conductivity `k_H`: 2nd-order approximation.
- Electron thermal conductivity `k_e`: 3rd-order approximation.
- Internal thermal conductivity `k_int`: Hirschfelder-Eucken type treatment.
- Reactive thermal conductivity `k_reac`: Butler-Brokaw type general expression.
- Electrical conductivity `sigma`: 3rd-order approximation (Devoto-type electron formulation).

Collision-integral handling:
- Ar-Ar: Aziz HFDTCS2 potential-based treatment.
- Ar-Ar+: elastic + charge-exchange combined for `l=1`.
  - `Q_ex = (A - B ln E)^2`
  - `Q^(1) = 2 Q_ex`
  - `Omega^(1,s) = sqrt((Omega_in^(1,s))^2 + (Omega_el^(1,s))^2)`
  - Constants:
    - `2Sigma`: `A=8.921`, `B=0.3960` (weight `1/3`)
    - `2Pi`: `A=6.189`, `B=0.2934` (weight `2/3`)
- e-Ar:
  - Momentum-transfer cross-section data (Milloy), high-energy completion (Frost).
  - High-order terms are not directly measured as standalone tables; derived through numerical treatment.
  - Assumption noted in LTEDiffusion summary: `Omega_ej^(2,2) = Omega_ej^(1,1)`.
- charged-charged:
  - Screened Coulomb model with Mason/Devoto tables.
  - Debye screening length based on electron density (electron-only screening convention).

## 3. Current implementation snapshot (code + metadata)

Evidence from current repo:
- Collision metadata: `data/processed/transport/argon_collision_integrals_metadata.json`
- Transport metadata: `data/processed/transport_properties/argon_transport_metadata.json`
- Collision builder: `scripts/murphy/build_collision_integrals.py`
- Transport builder: `scripts/murphy/build_transport_properties.py`

Current active profile/settings:
- `murphy_profile.enabled = true`
- `ar_ar_model = hfdtcs2_scatter`
- `ar_arp_elastic_model = aubreton_barker`
- `ar_arp_state_mix = boltzmann_normalized`
- `e_ar_high_order_blend = 1.0`
- e-Ar source includes LXCat merge (`Biagi`) in addition to Milloy/Frost.
- `debye_screening_density_model = electron_only`
- `k_reac_model = legacy_composite`

## 4. Conformance matrix

| Item | Murphy reference (LTEDiffusion) | Current implementation | Verdict |
|---|---|---|---|
| `mu` | 1st-order CE/Devoto-family | Devoto 1st-order heavy matrix (`build_transport_properties.py`) | Match |
| `k_H` | 2nd-order | Devoto/Muckenfuss-Curtiss 2nd-order matrix | Match |
| `k_e` | 3rd-order | Devoto 3rd-order electron Lee matrix | Match |
| `sigma` | 3rd-order | Devoto 3rd-order electron conductivity from Lee matrix inverse | Match |
| `k_int` | Hirschfelder-Eucken style | Hirschfelder-Eucken algebraic form | Match |
| `k_reac` | Butler-Brokaw general form | Default is `legacy_composite` (mass-gradient + cp amplification), not pure Butler-only mode | Partial match |
| Ar-Ar | Aziz HFDTCS2 | HFDTCS2 scattering active | Match |
| Ar-Ar+ charge exchange constants | Murphy constants and composition rule | Constants and composition are implemented (`AR_ARP_CX_STATES`, `sqrt(in^2+el^2)`) | Match |
| Ar-Ar+ state mixing detail | Murphy text/equation convention | Current mode uses `boltzmann_normalized`; code also has alternative `murphy_text` mode | Partial match |
| e-Ar data source set | Milloy + Frost basis | Milloy + Frost + LXCat(Biagi) merged | Partial match |
| e-Ar high-order handling | derived numerically from momentum-transfer treatment; not direct standalone measured tables | uses direct high-order moments with blend control; current `blend=1.0` | Partial match |
| charged-charged screening | screened Coulomb + electron-only Debye convention | electron-only Debye + Mason tables; includes practical caps/floors/scales | Match (with numerical stabilizers) |
| Ar2+/Ar3+/Ar4+ ion-neutral extension | not explicitly confirmed as Murphy canonical in LTEDiffusion summary | current default uses `sqrt(z)` scaling from Ar-Ar+ | Project-specific extension |

## 5. Direct code anchors

Collision builder:
- Murphy profile switch: `scripts/murphy/build_collision_integrals.py:390`
- `e_ar_high_order_blend = 1.0` in profile: `scripts/murphy/build_collision_integrals.py:397`
- Ar-Ar+ CX constants (`A,B,weights`): `scripts/murphy/build_collision_integrals.py:202`
- Ar-Ar+ composition (`sqrt(in^2 + el^2)`): `scripts/murphy/build_collision_integrals.py:1553`
- e-Ar high-order moment blending (`Q14/Q15`): `scripts/murphy/build_collision_integrals.py:1620`
- charged-charged screened Coulomb assembly: `scripts/murphy/build_collision_integrals.py:1650`

Transport builder:
- `k_reac_model` options and default (`legacy_composite`): `scripts/murphy/build_transport_properties.py:288`
- `k_reac` model switch logic: `scripts/murphy/build_transport_properties.py:1166`
- 3rd-order electron block and sigma construction: `scripts/murphy/build_transport_properties.py:1023`

## 6. Final assessment

Current Tutorial 01 is **Murphy-inspired and largely Murphy-consistent in transport formulation order**, but it is **not a strict Murphy reproduction** at present.

Main reasons:
- e-Ar input currently uses LXCat-merged data (not only Milloy/Frost basis).
- Ar-Ar+ state-mix is currently `boltzmann_normalized`, not the explicit Murphy-text mode.
- Reactive conductivity default is `legacy_composite`, not a pure Butler-only closure.
- Additional high-charge ion-neutral handling includes project-specific `sqrt(z)` extension.

So, for your question:
- "Is current method Murphy method?" -> **Partially yes (core framework), but not strict/pure Murphy in all data and closures.**
