# Phase 4: Collision Integrals (Transport Input)

## Scope

`scripts/murphy/build_collision_integrals.py` builds collision-integral tables for:

1. `Ar-Ar` (neutral-neutral)
2. `Ar-Ar+` (ion-neutral)
3. `e-Ar` (electron-neutral)
4. charged-charged pairs among `e-, Ar+, Ar2+, Ar3+, Ar4+`

Outputs are written under `data/processed/transport/`.

## Raw Data Files

Raw files expected in `data/raw/murphy/collision_integrals/argon/`:

- `aziz1990_hfdtcs2_constants.csv`
- `aubreton1986_ar_arplus_elastic_potentials.csv`
- `milloy1977_e_ar_momentum_transfer_table1.csv`
- `frost1964_e_ar_high_energy_digitized.csv`
- `mason1967_tableI_t2omega.csv`
- `mason1967_tableII_abc.csv`

## Implemented Models

### Ar-Ar

- Uses Aziz constants (`epsilon_over_k`, `sigma`) from `aziz1990`.
- Current implementation uses a documented LJ-style surrogate (`Omega11*`, `Omega22*` correlations) to generate `Q11`, `Q22`.

### Ar-Ar+

- Charge exchange:
  - `Qex = (A - B ln E)^2`
  - `Q^(1) = 2 * Qex`
  - weighted with Murphy & Tam (2014) constants for `2Sigma`, `2Pi`.
- Elastic component:
  - capture-model cross section for ion-induced dipole (`~ E^(-1/2)`).
- Combination:
  - `Omega(1,s) = sqrt(Omega_el(1,s)^2 + Omega_in(1,s)^2)` surrogate form.

### e-Ar

- Uses `Qm(E)` data from Milloy table (low energy) + Frost extension (high energy).
- Maxwellian energy averaging gives `Q11`.
- `Q22 = Q11` assumption applied.

### Charged-Charged

- Mason table interpolation (Table I/II):
  - Table I provides `(T*)^2 Omega*` terms.
  - Table II provides `A*`, `B*`, `C*`.
- Debye length uses electron density from phase-2 equilibrium output:
  - `lambda_D = sqrt(eps0 k_B T / (n_e e^2))`
- `T*` computed with electron-only screening.
- Low-ionization guardrails in code:
  - `n_e` floor: `1e14 m^-3` (default)
  - `lambda_D` cap: `1e-3 m` (default)

## Run

```bash
cd /home/brjeon/Proj_PlasmaPropCalc
python3 scripts/murphy/build_collision_integrals.py
```

or with project venv:

```bash
cd /home/brjeon/Proj_PlasmaPropCalc
.venv/bin/python scripts/murphy/build_collision_integrals.py
```

## Notes

- `frost1964` values are figure-digitized approximations.
- `Ar-Ar` currently uses a surrogate instead of full HFDTCS2 scattering integration.
