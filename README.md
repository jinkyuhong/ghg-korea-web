# 시도별 온실가스 대시보드 (WRF-Chem/DART)

WRF-Chem/DART 대기 역산으로 추정한 CO₂·CH₄ 농도, 배출량, 생태계 순교환량(NEE)을 17개 시도별로 보여 주는 정적 웹페이지입니다.

## 구성

| 경로 | 내용 |
|---|---|
| `index.html` | 대시보드 (단일 페이지, 서버 불필요) |
| `data/sido_data.json` | 시도별 집계 자료 (분기: 2019·2020년 / 월별 배출량: 2015–2020, 2025년) |
| `data/sido.geojson` | 지도용 시도 경계 (단순화) |
| `data/sigungu/<시도코드>.json`, `.geojson` | 시도별 시군구(250개) 자료와 경계. 시도를 고르면 그 시도 파일만 불러옴 |
| `scripts/build_data.py` | 원자료(NetCDF·npz) → `data/` 생성 스크립트 |
| `korea_sigungu.geojson` | 시군구 경계 원본 (2018년 기준, 시도 병합에 사용) |

원자료 폴더(`co2/ ch4/ concentration/ emissions/ uncertainty/`, `*.tar.gz`)는 용량이 커서 `.gitignore`로 제외했습니다.

## 자료 다시 만들기

```bash
pip install numpy netCDF4 shapely
python scripts/build_data.py
```

### N₂O (EDGAR 인벤토리)

N₂O는 역산 대상이 아니어서 EDGAR 격자 배출량(0.1°, 연 단위)만 사용합니다.

```bash
bash scripts/download_edgar_n2o.sh          # EDGAR 2025 GHG N2O 내려받기 → edgar_n2o/ (약 250 MB)
python scripts/build_data.py n2o pack       # 시도별 집계 후 data/sido_data.json 에 반영
```
`data/sido_data.json`에 N₂O 자료가 있으면 대시보드에 N₂O 탭이 자동으로 나타납니다.

## 내 컴퓨터에서 미리 보기

```bash
python3 -m http.server 8000     # 이 폴더에서 실행
# 브라우저에서 http://localhost:8000
```
(`index.html`을 더블클릭해 열면 브라우저 보안 정책 때문에 자료를 읽지 못합니다.)

## GitHub Pages 게시

```bash
git add index.html data scripts README.md .gitignore korea_sigungu.geojson
git commit -m "시도별 CO2·CH4 대시보드로 재구성"
git push origin main
```
GitHub 저장소 → Settings → Pages → Source를 `Deploy from a branch`, Branch를 `main` / `/ (root)`로 지정하면
`https://jinkyuhong.github.io/ghg-korea-web/` 에서 열립니다.

## 집계 방법 요약

- 시도 경계: 2018년 시군구 경계를 17개 시도로 병합
- 배출량·NEE: 격자(9 km) 값을 격자 안 육지에 균등 분포로 보고 시도의 육지 면적 비율로 배분해 합산
- 농도: 격자∩시도 면적 가중 평균 (모형 최하층), ±는 앙상블 표준편차의 지역 평균
- 분기 값은 분기 평균 배출 속도의 연환산(t/yr), 연 값은 네 분기 평균
- CH₄의 CO₂ 환산: GWP₁₀₀ = 28
