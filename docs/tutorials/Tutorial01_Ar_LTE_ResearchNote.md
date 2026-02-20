# Tutorial 01 Research Note: Ar LTE Property Pipeline

## 1. 목적

Ar LTE 조건에서 아래 물성치를 계산하는 재현 가능한 파이프라인을 정리한다.

- 평형 조성: `n_i(T, P)`
- 열역학 물성: `rho`, `h` (specific enthalpy), `Cp`
- 수송 물성: `mu`, `kappa`, `sigma`

적용 조건:

- 온도 `T = 300 ~ 30000 K`
- 압력 `P = 0.1, 1.0, 4.0 atm`
- Species: `{Ar, Ar+, Ar2+, Ar3+, Ar4+, e-}`

---

## 2. 계산 구성

### 2.1 Internal Partition Function

NIST ASD 에너지 준위 데이터(`Ar0.csv ~ Ar4.csv`)를 사용한다.

\[
Q_{\mathrm{int},s}(T)=\sum_j g_{s,j}\exp\!\left(-\frac{\epsilon_{s,j}}{k_B T}\right)
\]

- `g` 결측 행은 분배함수 합산에서 제외
- 결측 `g` 행의 `Level(eV)`는 해당 단계 이온화 에너지로 추출
- 전자 분배함수: `Q_int,e = 2`

출력:

- `data/processed/partition_functions/argon_partition_functions.csv`
- `data/processed/partition_functions/argon_ionization_energies.csv`

### 2.2 평형 조성 (Saha + EOS + Quasi-neutrality)

화학 퍼텐셜:

\[
\mu_i = k_B T \ln\!\left(\frac{n_i \lambda_i^3}{Q_{\mathrm{int},i}(T)}\right) + E_{\mathrm{ion},i},
\quad
\lambda_i=\frac{h}{\sqrt{2\pi m_i k_B T}}
\]

연립식:

- Saha 평형식 4개
- 전하중성식 1개
- 상태방정식 1개(필요 시 Debye-Huckel 압력 보정 포함)

출력:

- `data/processed/equilibrium/argon_equilibrium_0p1_1_4atm.csv`

### 2.3 열역학 물성

\[
\rho=\sum_i n_i m_i
\]

\[
H_i(T)\,[\mathrm{eV}] = \frac{5}{2}k_B T + k_B T^2\frac{d\ln Q_{\mathrm{int},i}}{dT} + E_{\mathrm{ion},i}
\]

\[
h=\frac{1}{\rho}\sum_i n_i H_i(T)\,[\mathrm{J/particle}],
\quad
C_p=\left(\frac{\partial h}{\partial T}\right)_P
\]

출력:

- `data/processed/thermo/argon_thermo_0p1_1_4atm.csv`

### 2.4 수송 계수

충돌적분:

- Ar-Ar
- Ar-Ar+
- e-Ar
- charged-charged (Coulomb, screening 적용)

수송계수:

- `mu` (점성)
- `kappa = k_H + k_e + k_int + k_reac`
- `sigma` (전기전도도)

출력:

- `data/processed/transport/argon_collision_integrals_all.csv`
- `data/processed/transport_properties/argon_transport_0p1_1_4atm.csv`

---

## 3. 결과 그래프 (GitHub 표시용)

### 3.1 Thermodynamic

`rho vs T`

![rho vs T](../../data/processed/thermo/plots/rho_vs_T.png)

`h vs T` (specific enthalpy)

![h vs T](../../data/processed/thermo/plots/h_vs_T.png)

`Cp vs T` (0~30000 K)

![Cp vs T 0-30000K](../../data/processed/thermo/plots/cp_vs_T_0_30000K.png)

### 3.2 Transport

`mu vs T`

![mu vs T](../../data/processed/transport_properties/plots/mu_vs_T.png)

`kappa vs T`

![kappa vs T](../../data/processed/transport_properties/plots/kappa_vs_T.png)

`sigma vs T`

![sigma vs T](../../data/processed/transport_properties/plots/sigma_vs_T.png)

`kappa components @ 1 atm`

![kappa components 1atm](../../data/processed/transport_properties/plots/kappa_components_1atm.png)

---

## 4. 재현 실행

프로젝트 루트에서 순서대로 실행:

```bash
python3 scripts/murphy/build_partition_functions.py
python3 scripts/murphy/solve_argon_lte_equilibrium.py
python3 scripts/murphy/build_thermo_properties.py
python3 scripts/murphy/build_collision_integrals.py
python3 scripts/murphy/build_transport_properties.py
```

---

## 5. 참고

- NIST ASD Energy Levels: https://physics.nist.gov/PhysRefData/ASD/levels_form.html
- Murphy 관련 기반 논문 및 표 데이터는 `data/raw/murphy/`와 `docs/methods/` 참조
