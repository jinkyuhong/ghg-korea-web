#!/bin/bash
# EDGAR 2025 GHG — N2O 격자 배출량(0.1°, 연 단위) 내려받기
#   총량(TOTALS): 2015–2024년 / 부문별: 2019·2020년
# 사용법:  bash scripts/download_edgar_n2o.sh      (저장소 루트에서 실행, 약 250 MB)
set -u
cd "$(dirname "$0")/.."
BASE="https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/EDGAR_2025_GHG/N2O"
OUT="edgar_n2o"; mkdir -p "$OUT"
SECTORS="AGS AWB CHE ENE IDE IND MNM N2O PRO_FFF PRU_SOL RCO REF_TRF SWD_INC SWD_LDF TNR_Aviation_CDS TNR_Aviation_CRS TNR_Aviation_LTO TNR_Aviation_SPS TNR_Other TNR_Ship TRO WWT"
get() {  # $1=부문 $2=연도
  f="EDGAR_2025_GHG_N2O_$2_$1_emi_nc.zip"
  [ -s "$OUT/$f" ] && return 0
  if curl -fsSL --retry 3 -o "$OUT/$f.part" "$BASE/$1/emi_nc/$f"; then mv "$OUT/$f.part" "$OUT/$f"; echo "  받음  $f"
  else rm -f "$OUT/$f.part"; echo "  없음  $f (건너뜀)"; fi
}
for y in 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024; do get TOTALS $y; done
for y in 2019 2020; do for s in $SECTORS; do get $s $y; done; done
echo "완료: $(ls "$OUT"/*.zip 2>/dev/null | wc -l) 개 파일 → $OUT/"
