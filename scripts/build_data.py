#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WRF-Chem/DART 역산 산출물 → 시도(광역)·시군구(기초)별 집계 → 웹 대시보드용 JSON 생성

입력 (저장소 루트 기준)
  korea_sigungu.geojson                      시군구 경계(2018년 기준) → 시도로 병합
  emissions/, concentration/, uncertainty/   Zenodo 공개본 분기 평균 NetCDF (2019, 2020)
  co2/, ch4/                                 월별 역산 산출물 npz (선택; 있으면 월별 시계열 생성)
  mlo/                                       마우나로아 월평균 CO2 (선택; scripts/download_mlo.sh 로 내려받음)
  edgar_n2o/                                 EDGAR N2O 격자 배출량 zip (선택; scripts/download_edgar_n2o.sh 로 내려받음)
출력
  data/sido.geojson      단순화한 시도 경계 (지도 표시용)
  data/sido_data.json    시도별 농도·배출량·NEE
  data/sigungu/<시도코드>.json / .geojson   시도별 시군구 자료와 경계

실행
  pip install numpy netCDF4 shapely
  python scripts/build_data.py            # 전체
  python scripts/build_data.py weights    # 단계별: weights | quarterly | monthly | n2o | pack
"""
import os, sys, json, glob, re, calendar, urllib.request
import numpy as np
import shapely
from shapely.geometry import shape, mapping
from shapely.ops import unary_union

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
CACHE = "scripts/.cache"
os.makedirs(CACHE, exist_ok=True)
os.makedirs("data", exist_ok=True)

SIDO = [  # (경계파일 코드 앞 2자리, 국문, 약칭, 영문)
    ("11", "서울특별시", "서울", "Seoul"), ("21", "부산광역시", "부산", "Busan"),
    ("22", "대구광역시", "대구", "Daegu"), ("23", "인천광역시", "인천", "Incheon"),
    ("24", "광주광역시", "광주", "Gwangju"), ("25", "대전광역시", "대전", "Daejeon"),
    ("26", "울산광역시", "울산", "Ulsan"), ("29", "세종특별자치시", "세종", "Sejong"),
    ("31", "경기도", "경기", "Gyeonggi"), ("32", "강원특별자치도", "강원", "Gangwon"),
    ("33", "충청북도", "충북", "Chungbuk"), ("34", "충청남도", "충남", "Chungnam"),
    ("35", "전북특별자치도", "전북", "Jeonbuk"), ("36", "전라남도", "전남", "Jeonnam"),
    ("37", "경상북도", "경북", "Gyeongbuk"), ("38", "경상남도", "경남", "Gyeongnam"),
    ("39", "제주특별자치도", "제주", "Jeju"),
]
YEARS = [2019, 2020]
NOSPLIT = {"29"}   # 하위 구역이 하나뿐인 세종특별자치시는 나누지 않고 시 전체로만 제공
EXTRA_YEARS = [2025]                 # 분기 공개본은 없고 월별 산출물만 있는 해 → 배출량과 MDV 보정 NEE만 제공
CELL_KM2 = 81.0                      # 9 km x 9 km
M_CO2, M_CH4 = 44.01, 16.04          # g/mol
NE_URL = ("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
          "master/geojson/ne_10m_admin_0_countries.geojson")
LAT0 = 36.5                          # 면적비 계산용 국지 투영 기준 위도


def proj(g):
    """경위도 → 국지 등면적 근사 평면(면적 '비율' 계산용)."""
    c = np.cos(np.deg2rad(LAT0))
    return shapely.transform(g, lambda xy: np.column_stack([xy[:, 0] * c, xy[:, 1]]))


def load_units():
    """경계 파일 → 시군구 목록 [{code, name, sido(시도 순번), geom}] (시도 순서, 코드 순서)."""
    gj = json.load(open("korea_sigungu.geojson", encoding="utf-8"))
    order = {c: i for i, (c, *_r) in enumerate(SIDO)}
    sgg = []
    for f in gj["features"]:
        pr = f["properties"]
        name = re.sub(r"^(.+?시)(.+구)$", r"\1 \2", pr["name"])          # '수원시장안구' → '수원시 장안구'
        sgg.append({"code": pr["code"], "name": name, "name_en": pr.get("name_eng", ""),
                    "sido": order[pr["code"][:2]], "geom": shapely.make_valid(shape(f["geometry"]))})
    sgg.sort(key=lambda u: (u["sido"], u["code"]))
    return sgg


def unit_weights(cells, cell_area, tree, sgg, others_geom):
    """격자 폴리곤 → (A, land).  A[행, 격자] = 격자 중 해당 지역이 차지하는 면적비.
    행 순서: 시도 17개 다음에 시군구(load_units 순서). 시도 값은 소속 시군구의 합이라 서로 정확히 맞는다."""
    Ag = np.zeros((len(sgg), len(cells)))
    for j, u in enumerate(sgg):
        gp = proj(u["geom"])
        idx = tree.query(gp, predicate="intersects")
        Ag[j, idx] = shapely.area(shapely.intersection(cells[idx], gp)) / cell_area[idx]
    As = np.zeros((len(SIDO), len(cells)))
    for j, u in enumerate(sgg):
        As[u["sido"]] += Ag[j]
    kor = As.sum(0)
    idx = tree.query(others_geom, predicate="intersects")
    oth = np.zeros(len(cells))
    oth[idx] = shapely.area(shapely.intersection(cells[idx], others_geom)) / cell_area[idx]
    oth = np.clip(oth, 0, 1 - np.minimum(kor, 1))   # 경계자료 간 중복 방지
    return np.vstack([As, Ag]), kor + oth


def other_land(x0, y0, x1, y1):
    """국외 육지(북한·일본·중국·러시아): 접경/연안 격자의 육지 면적 산정용."""
    ne_path = f"{CACHE}/ne_10m_admin_0_countries.geojson"
    if not os.path.exists(ne_path):
        print("  Natural Earth 국가 경계 내려받는 중…", flush=True)
        urllib.request.urlretrieve(NE_URL, ne_path)
    ne = json.load(open(ne_path, encoding="utf-8"))
    others = [shapely.make_valid(shape(f["geometry"])) for f in ne["features"]
              if f["properties"].get("ADM0_A3") in ("PRK", "JPN", "CHN", "RUS")]
    return proj(shapely.clip_by_rect(unary_union(others), x0, y0, x1, y1))


def simplify_geom(g, tol, min_area):
    s = g.simplify(tol, preserve_topology=True)
    parts = [p for p in (s.geoms if s.geom_type == "MultiPolygon" else [s]) if p.geom_type == "Polygon" and p.area > min_area]
    if not parts:
        parts = [max((s.geoms if s.geom_type == "MultiPolygon" else [s]), key=lambda p: p.area)]
    return shapely.set_precision(shapely.MultiPolygon(parts), 0.0001)


# ───────────────────────── 1. 격자-지역 가중치 ─────────────────────────
def build_weights():
    import netCDF4 as nc
    ds = nc.Dataset("emissions/wrfchem_dart_2020_quarterly_co2_emissions.nc")
    lat = np.array(ds["XLAT"][:], dtype=float)
    lon = np.array(ds["XLONG"][:], dtype=float)
    ny, nx = lat.shape

    def corners(a):
        p = np.pad(a, 1, mode="edge")
        p[0, :] = 2 * p[1, :] - p[2, :]; p[-1, :] = 2 * p[-2, :] - p[-3, :]
        p[:, 0] = 2 * p[:, 1] - p[:, 2]; p[:, -1] = 2 * p[:, -2] - p[:, -3]
        return 0.25 * (p[:-1, :-1] + p[1:, :-1] + p[:-1, 1:] + p[1:, 1:])
    clon, clat = corners(lon), corners(lat)
    ring = np.stack([
        np.stack([clon[:-1, :-1], clat[:-1, :-1]], -1), np.stack([clon[:-1, 1:], clat[:-1, 1:]], -1),
        np.stack([clon[1:, 1:], clat[1:, 1:]], -1), np.stack([clon[1:, :-1], clat[1:, :-1]], -1),
        np.stack([clon[:-1, :-1], clat[:-1, :-1]], -1)], axis=2).reshape(ny * nx, 5, 2)
    cells = proj(shapely.polygons(ring))
    cell_area = shapely.area(cells)

    sgg = load_units()
    tree = shapely.STRtree(cells)
    A, land = unit_weights(cells, cell_area, tree, sgg,
                           other_land(lon.min() - 1, lat.min() - 1, lon.max() + 1, lat.max() + 1))
    kor = A[:len(SIDO)].sum(0)

    # 지도용 단순화 경계: 시도(전국 지도) + 시도별 시군구(시도를 고르면 불러옴)
    feats = []
    os.makedirs("data/sigungu", exist_ok=True)
    for j, (code, name, short, en) in enumerate(SIDO):
        mine = [u for u in sgg if u["sido"] == j]
        s = simplify_geom(unary_union([u["geom"] for u in mine]), 0.004, 2e-4)
        c = max(s.geoms, key=lambda p: p.area).representative_point()
        feats.append({"type": "Feature", "properties": {"code": code, "name": name, "short": short, "name_en": en,
                      "cx": round(c.x, 3), "cy": round(c.y, 3)}, "geometry": mapping(s)})
        if code in NOSPLIT:
            continue
        small = max(u["geom"].area for u in mine) < 0.02      # 특별·광역시 자치구는 더 세밀하게
        tol, amin = (0.0006, 3e-6) if small else (0.002, 3e-5)
        sub = [{"type": "Feature", "properties": {"code": u["code"], "name": u["name"]},
                "geometry": mapping(simplify_geom(u["geom"], tol, amin))} for u in mine]
        json.dump({"type": "FeatureCollection", "features": sub}, open(f"data/sigungu/{code}.geojson", "w", encoding="utf-8"),
                  ensure_ascii=False, separators=(",", ":"))
    json.dump({"type": "FeatureCollection", "features": feats}, open("data/sido.geojson", "w", encoding="utf-8"),
              ensure_ascii=False, separators=(",", ":"))
    np.savez(f"{CACHE}/weights.npz", A=A, land=land, kor=kor, ny=ny, nx=nx)
    json.dump([{k: u[k] for k in ("code", "name", "name_en", "sido")} for u in sgg], open(f"{CACHE}/units.json", "w"), ensure_ascii=False)
    print("  가중치 저장:", A.shape, "(시도", len(SIDO), "+ 시군구", len(sgg), ") 한국 육지 격자합", kor.sum().round(1))


def load_weights():
    w = np.load(f"{CACHE}/weights.npz")
    A, land = w["A"], w["land"]
    # 면적 가중(농도 평균용): 격자∩시도 면적비
    # 플럭스 가중(배출량·NEE 합산용): 격자 플럭스를 격자 내 '육지'에 균등 분포한다고 보고 시도 면적비로 배분.
    #   (연안 격자의 발전소·제철소 등 육상 배출이 바다 면적만큼 누락되는 것을 방지)
    #   육지가 5% 미만인 격자(작은 섬)는 해상 배출의 과대 배분을 막기 위해 분모를 0.05로 하한 처리.
    WF = A / np.maximum(land, 0.05)[None, :]
    return A, WF


# ───────────────────────── 2. 분기 자료(Zenodo 공개본) ─────────────────────────
def build_quarterly():
    import netCDF4 as nc
    A, WF = load_weights()
    Asum = A.sum(1)
    out = {"years": YEARS, "sectors": {}, "data": {}}

    def flat(v):
        return np.nan_to_num(np.ma.filled(np.ma.masked_invalid(v), 0.0).astype(float)).reshape(v.shape[0], -1)

    for yr in YEARS:
        ec = nc.Dataset(f"emissions/wrfchem_dart_{yr}_quarterly_co2_emissions.nc")
        em = nc.Dataset(f"emissions/wrfchem_dart_{yr}_quarterly_ch4_emissions.nc")
        cc = nc.Dataset(f"concentration/wrfchem_dart_{yr}_quarterly_mean_co2.nc")
        cm = nc.Dataset(f"concentration/wrfchem_dart_{yr}_quarterly_mean_ch4.nc")
        uc = nc.Dataset(f"uncertainty/wrfchem_dart_{yr}_quarterly_mean_co2_uncertainty.nc")
        um = nc.Dataset(f"uncertainty/wrfchem_dart_{yr}_quarterly_mean_ch4_uncertainty.nc")
        sc = sorted(k[6:] for k in ec.variables if k.startswith("E_CO2_") and k != "E_CO2_total")
        sm = sorted(k[6:] for k in em.variables if k.startswith("E_CH4_") and k != "E_CH4_total")
        out["sectors"] = {"co2": sc, "ch4": sm}
        flux = lambda ds, k: flat(ds[k][:]) @ WF.T                 # (4, nsido) ton/yr (연환산율)
        mean = lambda ds, k: (flat(ds[k][:, 0]) @ A.T) / Asum      # (4, nsido) 최하층 면적가중 평균
        # 배출량 불확실성: 격자별 앙상블 표준편차 σ를 제곱합의 제곱근으로 시도 합산  sqrt(Σ (w·σ)²)
        sd_rss = lambda ds, k: np.sqrt((flat(ds[k][:]) ** 2) @ (WF.T ** 2))
        y = {
            "co2": {"emis": flux(ec, "E_CO2_total"), "nee": flux(ec, "NEE"), "gee": flux(ec, "GEE"),
                    "res": flux(ec, "RES"), "conc": mean(cc, "co2"), "conc_sd": mean(uc, "co2_uncertainty"),
                    "emis_sd": sd_rss(uc, "E_CO2_uncertainty"),
                    "sec": np.stack([flux(ec, "E_CO2_" + s) for s in sc], -1)},
            "ch4": {"emis": flux(em, "E_CH4_total"), "conc": mean(cm, "ch4"), "conc_sd": mean(um, "ch4_uncertainty"),
                    "emis_sd": sd_rss(um, "E_CH4_uncertainty"),
                    "sec": np.stack([flux(em, "E_CH4_" + s) for s in sm], -1)},
        }
        out["data"][str(yr)] = {g: {k: v.tolist() for k, v in d.items()} for g, d in y.items()}
        n = len(SIDO)
        print(f"  {yr}: 17개 시도 합계 CO2 {y['co2']['emis'].mean(0)[:n].sum()/1e6:.1f} Mt/yr, "
              f"NEE {y['co2']['nee'].mean(0)[:n].sum()/1e6:.1f} Mt/yr, CH4 {y['ch4']['emis'].mean(0)[:n].sum()/1e3:.1f} kt/yr"
              f"  (시군구 합 CO2 {y['co2']['emis'].mean(0)[n:].sum()/1e6:.1f})")
    out["area_km2"] = (Asum * CELL_KM2).tolist()
    json.dump(out, open(f"{CACHE}/quarterly.json", "w"))


# ───────────────────────── 3. 월별 자료(npz) ─────────────────────────
def build_monthly():
    A, WF = load_weights()
    res = {"months": [], "co2_emis": [], "nee": [], "ch4_emis": []}
    files = sorted(glob.glob("co2/co2_*_yon.npz"))
    for f in files:
        ym = re.search(r"_(\d{6})_", f).group(1)
        f4 = f"ch4/ch4_{ym}_yon.npz"
        if not os.path.exists(f4):
            continue
        hours = calendar.monthrange(int(ym[:4]), int(ym[4:]))[1] * 24

        def mmean(z, k):
            sel = np.array([os.path.basename(str(p))[:6] == ym for p in z["flist"]])   # 해당 월 분석시각만
            if sel.sum() == 0:
                sel[:] = True
            return np.nan_to_num(z[k][sel].astype(float)).mean(0).ravel(), int(sel.sum())
        zc, zm = np.load(f), np.load(f4)
        e, n = mmean(zc, "ean"); nee, _ = mmean(zc, "nee"); e4, _ = mmean(zm, "ean")
        k = CELL_KM2 * hours / 1e6        # mol km-2 hr-1 → ton/월 (× 분자량)
        res["months"].append(ym)
        res["co2_emis"].append((WF @ e * k * M_CO2).tolist())
        res["nee"].append((WF @ nee * k * M_CO2).tolist())
        res["ch4_emis"].append((WF @ e4 * k * M_CH4).tolist())
        ns = len(SIDO)
        print("  ", ym, n, "steps", f"CO2 {sum(res['co2_emis'][-1][:ns])/1e6:.1f} Mt  NEE {sum(res['nee'][-1][:ns])/1e6:.1f} Mt  "
              f"CH4 {sum(res['ch4_emis'][-1][:ns])/1e3:.1f} kt", flush=True)
    json.dump(res, open(f"{CACHE}/monthly.json", "w"))


# ───────────────────────── 3a. 분기 공개본이 없는 해의 생태계 흡수량(NEE): MDV 보정 ─────────────────────────
def build_mdv():
    """월별 산출물(npz)의 NEE는 6시간 간격 분석시각(하루 4회)만의 평균이라 일변화가 제대로 반영되지 않는다.
    분기 공개본(매시 평균, 2019·2020년)과 같은 해 npz 분기 평균의 차이를 평균 일변화(MDV) 보정량으로 보고,
    두 해의 보정량을 격자·분기별로 평균해 EXTRA_YEARS(2025년)의 npz 분기 평균 NEE에 더한다."""
    import netCDF4 as nc
    A, WF = load_weights()
    fac = CELL_KM2 * 8760 / 1e6 * M_CO2                      # mol km-2 hr-1 → t/격자/yr

    def rawq(yr):
        out = np.zeros((4, WF.shape[1]))
        for q in range(4):
            H = 0
            for mm in range(3 * q + 1, 3 * q + 4):
                z = np.load(f"co2/co2_{yr}{mm:02d}_yon.npz"); ym = f"{yr}{mm:02d}"
                sel = np.array([os.path.basename(str(p))[:6] == ym for p in z["flist"]])
                if sel.sum() == 0:
                    sel[:] = True
                h = calendar.monthrange(yr, mm)[1] * 24; H += h
                out[q] += np.nan_to_num(z["nee"][sel].astype(float)).mean(0).ravel() * h
            out[q] *= fac / H
        return out

    def zen(yr):
        v = nc.Dataset(f"emissions/wrfchem_dart_{yr}_quarterly_co2_emissions.nc")["NEE"][:]
        return np.nan_to_num(np.ma.filled(np.ma.masked_invalid(v), 0.0).astype(float)).reshape(4, -1)

    d = [zen(y) - rawq(y) for y in YEARS]
    mdv = np.mean(d, 0)                                       # (4, 격자) 2019·2020년 평균 MDV 보정량
    ns = len(SIDO); res = {}
    for i, y in enumerate(YEARS):
        print(f"  {y} MDV 보정량(전국, Mt/yr, 분기별):", np.round((d[i] @ WF.T)[:, :ns].sum(1) / 1e6, 1).tolist())
    print("  평균 MDV 보정량:", np.round((mdv @ WF.T)[:, :ns].sum(1) / 1e6, 1).tolist())
    for y in EXTRA_YEARS:
        if not all(os.path.exists(f"co2/co2_{y}{mm:02d}_yon.npz") for mm in range(1, 13)):
            continue
        raw = rawq(y); cor = (raw + mdv) @ WF.T              # (4, 지역) t/yr
        res[str(y)] = cor.tolist()
        print(f"  {y} NEE 보정 전:", np.round((raw @ WF.T)[:, :ns].sum(1) / 1e6, 1).tolist(),
              "→ 보정 후:", np.round(cor[:, :ns].sum(1) / 1e6, 1).tolist(), f"연평균 {cor[:, :ns].sum(1).mean() / 1e6:.1f} Mt/yr")
    json.dump(res, open(f"{CACHE}/nee_extra.json", "w"))

# ───────────────────────── 3b. N2O (EDGAR 인벤토리) ─────────────────────────
EDGAR_DIR = "edgar_n2o"
EDGAR_CODE = {"PRO_FFF": "PRO", "PRU_SOL": "PRU", "REF_TRF": "REF", "SWD_INC": "INC", "SWD_LDF": "LDF",
              "TNR_Aviation_CDS": "CDS", "TNR_Aviation_CRS": "CRS", "TNR_Aviation_LTO": "LTO",
              "TNR_Aviation_SPS": "SPS", "TNR_Other": "TNR", "TNR_Ship": "SHI"}
BBOX = (124.0, 32.5, 132.5, 39.5)     # 한반도 남부 주변만 잘라 사용


def _edgar_read(path):
    """EDGAR emi_nc zip → (lat, lon, 배출량[t/격자/yr]) (BBOX 범위)."""
    import netCDF4 as nc, zipfile
    with zipfile.ZipFile(path) as z:
        name = [n for n in z.namelist() if n.endswith(".nc")][0]
        ds = nc.Dataset("inmem.nc", memory=z.read(name))
    lat = np.array(ds["lat"][:], float); lon = np.array(ds["lon"][:], float)
    var = [v for v in ds.variables.values() if v.ndim == 2][0]
    units = str(getattr(var, "units", "")).lower()
    if not units.startswith("ton"):
        raise ValueError(f"{path}: 예상과 다른 단위 '{units}' (emi_nc의 Tonnes 필요)")
    lon = np.where(lon > 180, lon - 360, lon)
    jj = np.where((lat >= BBOX[1]) & (lat <= BBOX[3]))[0]; ii = np.where((lon >= BBOX[0]) & (lon <= BBOX[2]))[0]
    a = np.ma.filled(var[jj.min():jj.max() + 1, :][:, ii], 0.0).astype(float)
    return lat[jj], lon[ii], np.nan_to_num(a)


def _edgar_weights(lat, lon):
    """0.1° 격자 → 지역별 플럭스 가중치(격자 내 육지에 균등 분포 가정; 역산 격자와 같은 방식). 행 순서는 unit_weights 와 같다."""
    d = 0.05
    LON, LAT = np.meshgrid(lon, lat)
    ring = np.stack([np.stack([LON - d, LAT - d], -1), np.stack([LON + d, LAT - d], -1), np.stack([LON + d, LAT + d], -1),
                     np.stack([LON - d, LAT + d], -1), np.stack([LON - d, LAT - d], -1)], axis=2).reshape(-1, 5, 2)
    cells = proj(shapely.polygons(ring)); area = shapely.area(cells); tree = shapely.STRtree(cells)
    A, land = unit_weights(cells, area, tree, load_units(), other_land(BBOX[0] - 1, BBOX[1] - 1, BBOX[2] + 1, BBOX[3] + 1))
    return A / np.maximum(land, 0.05)[None, :]


def build_n2o():
    files = sorted(glob.glob(f"{EDGAR_DIR}/EDGAR_*_N2O_*_emi_nc.zip"))
    if not files:
        print(f"  {EDGAR_DIR}/ 에 EDGAR 파일이 없습니다 → scripts/download_edgar_n2o.sh 를 먼저 실행하세요. (N2O 건너뜀)")
        return
    parse = lambda f: re.search(r"_N2O_(\d{4})_(.+)_emi_nc\.zip$", os.path.basename(f)).groups()
    version = re.match(r"(EDGAR_\d{4}_GHG|v[\d.]+_FT\d{4}_GHG)", os.path.basename(files[0])).group(1).replace("_", " ")
    WF = None
    tot, sec = {}, {}
    for f in files:
        yr, sname = parse(f)
        lat, lon, a = _edgar_read(f)
        if WF is None:
            WF = _edgar_weights(lat, lon)
        v = WF @ a.ravel()                                    # 시도별 t N2O/yr
        if sname == "TOTALS":
            tot[int(yr)] = v
        else:
            sec.setdefault(yr, {})[EDGAR_CODE.get(sname, sname)] = v
    years = sorted(tot)
    codes = sorted({c for y in sec.values() for c in y})
    out = {"source": version, "years": years, "sectors": codes,
           "total": [np.round(tot[y], 1).tolist() for y in years],
           "sec": {y: np.round(np.stack([d.get(c, np.zeros(len(tot[years[0]]))) for c in codes], -1), 2).tolist() for y, d in sec.items()}}
    json.dump(out, open(f"{CACHE}/n2o.json", "w"))
    for y in years:
        ns = len(SIDO)
        chk = f"  부문 합/총량 {sum(v[:ns].sum() for v in sec[str(y)].values()) / tot[y][:ns].sum():.3f}" if str(y) in sec else ""
        print(f"  {y}: 17개 시도 합계 N2O {tot[y][:ns].sum() / 1e3:.1f} kt/yr{chk}")


# ───────────────────────── 3c. 마우나로아 CO2 (배경대기 비교용) ─────────────────────────
def read_mlo():
    """mlo/ 폴더의 월평균 CO2 파일 → {연도: [분기 평균 4개]} (ppm). 파일이 없으면 None.
    지원: NOAA GML co2_mm_mlo.txt/.csv (4번째 열 = 월평균), Scripps monthly_in_situ_co2_mlo.csv (5번째 열)."""
    files = sorted(glob.glob("mlo/*mlo*.txt") + glob.glob("mlo/*mlo*.csv"))
    if not files:
        return None
    f = files[0]
    scripps = "in_situ" in os.path.basename(f)
    col = 4 if scripps else 3
    mon = {}
    for line in open(f, encoding="utf-8", errors="ignore"):
        t = [x for x in re.split(r"[,\s]+", line.strip()) if x]
        try:
            y, m, v = int(t[0]), int(t[1]), float(t[col])
        except (ValueError, IndexError):
            continue
        if 1 <= m <= 12 and v > 0:
            mon[(y, m)] = v
    out = {}
    for y in YEARS:
        if all((y, m) in mon for m in range(1, 13)):
            out[str(y)] = [round(float(np.mean([mon[(y, 3 * q + k)] for k in (1, 2, 3)])), 2) for q in range(4)]
    if not out:
        return None
    src = "Scripps CO2 Program" if scripps else "NOAA Global Monitoring Laboratory"
    return {"source": src, "file": os.path.basename(f), "co2": out}


# ───────────────────────── 4. 웹용 JSON 묶기 ─────────────────────────
def pack():
    """웹용 JSON 묶기: 전국(시도 17개) 파일 + 시도별 시군구 파일(data/sigungu/<시도코드>.json)."""
    q = json.load(open(f"{CACHE}/quarterly.json"))
    mpath = f"{CACHE}/monthly.json"
    m = json.load(open(mpath)) if os.path.exists(mpath) else None
    n2o = json.load(open(f"{CACHE}/n2o.json")) if os.path.exists(f"{CACHE}/n2o.json") else None
    units = json.load(open(f"{CACHE}/units.json"))
    ns = len(SIDO)
    r = lambda a, n=0: np.round(np.asarray(a, float), n).tolist()

    def bundle(idx, small):
        """지역 행 번호 idx 만 골라 웹용 구조로. small=True(시군구)면 자릿수를 더 살린다."""
        o = {"quarterly": {}}
        for yr, d in q["data"].items():
            o["quarterly"][yr] = {}
            for g in ("co2", "ch4"):
                dd = {}
                for k, v in d[g].items():
                    if k in ("gee", "res") and small:
                        continue                                  # 시군구 파일에는 싣지 않음(화면에서 쓰지 않음)
                    a = np.asarray(v, float)
                    a = a[:, idx] if a.ndim >= 2 else a
                    dp = (2 if g == "co2" else 1) if k.startswith("conc") else (0 if g == "co2" else (2 if small else 1))
                    dd[k] = r(a, dp)
                o["quarterly"][yr][g] = dd
        if m:
            # 월별 NEE는 Zenodo 분기 공개본과 값 차이가 커서 웹 자료에는 싣지 않음(캐시에만 보관)
            mc, m4 = np.asarray(m["co2_emis"], float)[:, idx], np.asarray(m["ch4_emis"], float)[:, idx]
            o["monthly"] = {"months": m["months"], "co2_emis": r(mc), "ch4_emis": r(m4, 2 if small else 1)}
            # 분기 공개본이 없는 해(EXTRA_YEARS)는 월별 산출물을 분기로 묶어 배출량만 제공 (t/yr, 연환산)
            ex = {}
            for y in EXTRA_YEARS:
                if not all(f"{y}{mm:02d}" in m["months"] for mm in range(1, 13)):
                    continue
                yd = 366 if calendar.isleap(y) else 365
                ex[str(y)] = {}
                for key, arr, dp in (("co2_emis", mc, 0), ("ch4_emis", m4, 2 if small else 1)):
                    rows = []
                    for qq in range(4):
                        mms = range(3 * qq + 1, 3 * qq + 4)
                        qsum = sum(arr[m["months"].index(f"{y}{mm:02d}")] for mm in mms)
                        rows.append(qsum * yd / sum(calendar.monthrange(y, mm)[1] for mm in mms))
                    ex[str(y)][key] = r(rows, dp)
            nx = json.load(open(f"{CACHE}/nee_extra.json")) if os.path.exists(f"{CACHE}/nee_extra.json") else {}
            for y, v in nx.items():                               # MDV 보정한 생태계 흡수량(NEE)
                if y in ex:
                    ex[y]["co2_nee"] = r(np.asarray(v, float)[:, idx], 0)
            if ex:
                o["extra_quarterly"] = ex
        if n2o:
            o["n2o"] = {"source": n2o["source"], "years": n2o["years"], "sectors": n2o["sectors"],
                        "total": r(np.asarray(n2o["total"], float)[:, idx], 2 if small else 1),
                        "sec": {y: r(np.asarray(v, float)[idx], 3 if small else 2) for y, v in n2o["sec"].items()}}
        return o

    cells = np.asarray(q["area_km2"], float) / CELL_KM2          # 지역이 덮는 격자 수(농도 평균 가중치)
    out = {
        "meta": {"model": "WRF-Chem/DART", "grid_km": 9, "years": q["years"], "sectors": q["sectors"],
                 "boundary_base_year": 2018,
                 "units": {"flux_quarterly": "t/yr (분기 평균의 연환산율)", "flux_monthly": "t/월",
                           "co2_conc": "ppm", "ch4_conc": "ppb"}},
        "sido": [{"code": c, "name": n, "short": s, "name_en": e, "w": round(cells[i], 2),
                  "n_sub": 0 if c in NOSPLIT else sum(1 for u in units if u["sido"] == i)} for i, (c, n, s, e) in enumerate(SIDO)],
    }
    out.update(bundle(list(range(ns)), False))
    mlo = read_mlo()
    if mlo:
        out["mlo"] = mlo
        print("  마우나로아:", mlo["source"], mlo["co2"])
    json.dump(out, open("data/sido_data.json", "w", encoding="utf-8"), ensure_ascii=False, separators=(",", ":"))
    print("  data/sido_data.json", os.path.getsize("data/sido_data.json") // 1024, "KB /",
          "data/sido.geojson", os.path.getsize("data/sido.geojson") // 1024, "KB")

    os.makedirs("data/sigungu", exist_ok=True)
    tot = 0
    for i, (c, n, s, e) in enumerate(SIDO):
        if c in NOSPLIT:
            continue
        idx = [ns + j for j, u in enumerate(units) if u["sido"] == i]
        sub = {"parent": c, "sido": [{"code": units[j - ns]["code"], "name": units[j - ns]["name"],
                                      "name_en": units[j - ns]["name_en"], "w": round(cells[j], 3)} for j in idx]}
        sub.update(bundle(idx, True))
        path = f"data/sigungu/{c}.json"
        json.dump(sub, open(path, "w", encoding="utf-8"), ensure_ascii=False, separators=(",", ":"))
        tot += os.path.getsize(path) + os.path.getsize(f"data/sigungu/{c}.geojson")
    print(f"  data/sigungu/ 시도별 시군구 자료·경계 {len(SIDO) - len(NOSPLIT)}쌍, 합계 {tot // 1024} KB")


if __name__ == "__main__":
    steps = sys.argv[1:] or ["weights", "quarterly", "monthly", "mdv", "n2o", "pack"]
    for s in steps:
        print(f"=== {s} ===", flush=True)
        {"weights": build_weights, "quarterly": build_quarterly, "monthly": build_monthly, "mdv": build_mdv, "n2o": build_n2o, "pack": pack}[s]()
