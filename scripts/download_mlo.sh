#!/bin/bash
# 마우나로아(Mauna Loa) 월평균 CO2 농도 내려받기 — NOAA Global Monitoring Laboratory
# 사용법:  bash scripts/download_mlo.sh      (저장소 루트에서 실행)
set -u
cd "$(dirname "$0")/.."
mkdir -p mlo
curl -fsSL --retry 3 -o mlo/co2_mm_mlo.txt "https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_mlo.txt" \
  && echo "받음: mlo/co2_mm_mlo.txt ($(grep -vc '^#' mlo/co2_mm_mlo.txt) 개월)" || echo "내려받기 실패"
