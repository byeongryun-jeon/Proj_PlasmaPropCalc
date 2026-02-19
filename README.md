# Proj_PlasmaPropCalc

Murphy 방법과 Mutation++를 비교하여 플라즈마 물성치를 계산/검증하는 튜토리얼 프로젝트입니다.

## 목표

아래 6개 주제를 동일한 형식으로 재현 가능하게 정리합니다.

1. Murphy 방법을 이용한 Ar 물성치 획득
2. Mutation++를 이용한 Ar 물성치 획득
3. Murphy 방법을 이용한 Ar-He 물성치 획득
4. Mutation++를 이용한 Ar-He 물성치 획득
5. Diffusion coefficient 계산
6. Ar-He-H2 삼중가스 계산

## 디렉터리 구조

```text
Proj_PlasmaPropCalc/
├─ docs/
│  ├─ references/            # 논문/출처 정리
│  └─ methods/               # 계산식, 가정(LTE, species set, collision data) 정리
├─ data/
│  ├─ raw/
│  │  ├─ murphy/             # 원문 표/수치/스크랩 등 원시 데이터
│  │  └─ mutationpp/         # Mutation++ 입력/원시 출력
│  └─ processed/             # 후처리된 통합 데이터(csv, parquet 등)
├─ cases/
│  ├─ 01_ar_murphy/
│  ├─ 02_ar_mutationpp/
│  ├─ 03_ar_he_murphy/
│  ├─ 04_ar_he_mutationpp/
│  ├─ 05_diffusion_coefficient/
│  └─ 06_ar_he_h2_ternary/
│     ├─ input/              # 케이스 입력 파일
│     ├─ work/               # 중간 산출물
│     └─ output/             # 최종 결과
├─ scripts/
│  ├─ murphy/                # Murphy 방법 계산/추출 스크립트
│  ├─ mutationpp/            # Mutation++ 실행/파싱 스크립트
│  └─ common/                # 공통 유틸(단위변환, plotting, compare)
├─ configs/
│  ├─ species/               # species 목록 및 조성 설정
│  ├─ transport/             # 충돌적분/수송 관련 설정
│  └─ run/                   # 실행 파라미터(T, P, 조성 sweep)
├─ notebooks/                # 탐색/검증용 Jupyter notebook
├─ results/
│  ├─ tables/                # 비교표 (Murphy vs Mutation++)
│  ├─ figures/               # 그래프
│  └─ logs/                  # 실행 로그
├─ tests/                    # 계산값 검증 테스트
└─ env/                      # 환경 설정(예: conda, pip, module notes)
```

## 권장 작업 흐름

1. `docs/references`에 Murphy 논문 및 사용식 출처를 먼저 정리
2. `configs/species`, `configs/run`에 공통 조건(T, P, 조성) 정의
3. `cases/01`~`cases/06` 순서로 입력 생성 및 계산 실행
4. 결과를 `data/processed`와 `results/tables`에 표준 형식으로 저장
5. `scripts/common`으로 Murphy vs Mutation++ 비교 그래프 생성

## 첫 실행 체크리스트

1. 단위 통일: `T [K]`, `P [Pa]`, `rho [kg/m^3]`, `Cp [J/(kg K)]`, `h or en [J/kg]`, `mu [Pa s]`, `kappa [W/(m K)]`, `sigma [S/m]`
2. 조성 기준 통일: mole fraction vs mass fraction
3. species 집합 일치: Ar, Ar+, Ar2+, ... / He, H2, H, H+, e- 등 케이스별 명시
4. 동일한 온도/압력 그리드에서 Murphy와 Mutation++ 비교

## 다음 단계

원하면 다음으로 아래를 바로 생성할 수 있습니다.

1. 케이스별 입력 템플릿(`input.yaml`) 6개
2. 공통 물성치 출력 포맷(`results/tables/schema.md`)
3. 자동 실행 스크립트 골격(`scripts/run_all.sh`)
