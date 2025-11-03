# import requests, csv, io
# from collections import defaultdict

# import json, os, time
# from datetime import datetime, timezone

# CACHE_DIR  = ".macro_cache"
# CACHE_FILE = "cpi_2015eq100.json"

# #%% Get CPI Data for EUR and USD

# def euro_hicp_annual(start=1990, end=2025, partial_ok=True, min_months=1):
#     """
#     Fetch Euro area HICP index (2015=100) annual averages from ECB.
#     - partial_ok=True: include partial years (e.g., 2025 with available months)
#     - partial_ok=False: require 12 months (full year)
#     """
#     url = "https://data-api.ecb.europa.eu/service/data/ICP/M.U2.N.000000.4.INX"
#     params = {"format": "csvdata", "startPeriod": str(start), "endPeriod": str(end)}
#     r = requests.get(url, params=params, timeout=15)
#     r.raise_for_status()

#     rows = list(csv.reader(io.StringIO(r.text)))
#     header = rows[0]
#     ti, vi = header.index("TIME_PERIOD"), header.index("OBS_VALUE")

#     by_year = defaultdict(list)
#     for row in rows[1:]:
#         if len(row) <= vi: continue
#         t, v = row[ti], row[vi]
#         if not t or not v: continue
#         y = int(t[:4])
#         by_year[y].append(float(v))

#     out = {}
#     for y, vals in by_year.items():
#         needed = 12 if not partial_ok else max(min_months, 1)
#         if len(vals) >= needed:
#             out[y] = round(sum(vals)/len(vals), 2)
#     return dict(sorted(out.items()))

# def us_cpi_annual(start=1990, end=2025, rebase_to_2015=True, partial_ok=True, min_months=1):
#     """
#     U.S. CPI (All items, Index 2015=100) via OECD monthly -> annual averages.

#     Args:
#         start, end (int): year range to return.
#         rebase_to_2015 (bool): rebase to 2015 = 100 (OECD already is 2015=100, so this is minor).
#         partial_ok (bool): include partial years (e.g., 2025) if at least `min_months` months exist.
#         min_months (int): minimum months required to include a partial year.

#     Returns:
#         dict[int, float]: {year: annual_average_index_2015eq100}
#     """
#     url = "https://stats.oecd.org/sdmx-json/data/PRICES_CPI/USA.CPALTT01.IXOB.M/all"
#     params = {"contentType": "csv"}

#     r = requests.get(url, params=params, timeout=15)
#     r.raise_for_status()

#     rows = list(csv.reader(io.StringIO(r.text)))
#     if not rows:
#         return {}

#     header = rows[0]

#     # --- Make column detection bulletproof ---
#     time_keys = ["TIME", "Time", "TIME_PERIOD"]
#     value_keys = ["Value", "OBS_VALUE"]

#     ti = next((header.index(k) for k in time_keys if k in header), None)
#     vi = next((header.index(k) for k in value_keys if k in header), None)

#     if ti is None or vi is None:
#         raise RuntimeError(f"Could not find time/value columns in OECD header: {header}")

#     # --- Group monthly values by year ---
#     by_year = defaultdict(list)
#     for row in rows[1:]:
#         if len(row) <= max(ti, vi):
#             continue
#         t, v = row[ti], row[vi]
#         if not t or not v:
#             continue
#         try:
#             y = int(t[:4])
#             val = float(v)
#         except ValueError:
#             continue
#         if y < start or y > end:
#             continue
#         by_year[y].append(val)

#     # --- Compute annual averages ---
#     out = {}
#     for y, vals in by_year.items():
#         needed = 12 if not partial_ok else max(1, int(min_months))
#         if len(vals) >= needed:
#             out[y] = round(sum(vals) / len(vals), 2)

#     # --- Optional rebasing ---
#     if rebase_to_2015 and 2015 in out and out[2015] not in (0, None):
#         factor = 100.0 / out[2015]
#         out = {y: round(v * factor, 2) for y, v in out.items()}

#     return dict(sorted(out.items()))

# #%% Create/update cache

# def _ensure_dir(path: str):
#     os.makedirs(path, exist_ok=True)

# def _cache_path(dirpath: str = CACHE_DIR, filename: str = CACHE_FILE) -> str:
#     _ensure_dir(dirpath)
#     return os.path.join(dirpath, filename)

# def load_cpi_cache(dirpath: str = CACHE_DIR, filename: str = CACHE_FILE) -> dict:
#     path = _cache_path(dirpath, filename)
#     if not os.path.exists(path):
#         return {}
#     with open(path, "r", encoding="utf-8") as f:
#         return json.load(f)

# def save_cpi_cache(payload: dict, dirpath: str = CACHE_DIR, filename: str = CACHE_FILE) -> str:
#     path = _cache_path(dirpath, filename)
#     with open(path, "w", encoding="utf-8") as f:
#         json.dump(payload, f, ensure_ascii=False, indent=2, sort_keys=True)
#     return path

# def update_cpi_cache(
#     start_year: int = 1990,
#     end_year: int | None = None,
#     include_partial_latest: bool = True,
#     dirpath: str = CACHE_DIR,
#     filename: str = CACHE_FILE,
# ) -> tuple[str, int]:
#     """
#     Create/update JSON cache with:
#       - EUR: ECB HICP (2015=100) annual averages (monthly -> annual)
#       - USD: World Bank CPI rebased to 2015=100
#     Tries current year, falls back to current_year-1 if needed.
#     Returns (cache_path, chosen_end_year).
#     """
#     cy = datetime.now(timezone.utc).year
#     targets = [end_year or cy, (end_year or cy) - 1]

#     last_err = None
#     for cand in targets:
#         try:
#             eur = euro_hicp_annual(start=start_year, end=cand, partial_ok=include_partial_latest)
#             usd = us_cpi_annual(start=start_year, end=cand, rebase_to_2015=True)

#             # merge with any existing cache
#             cache = load_cpi_cache(dirpath, filename)
#             payload = {
#                 "meta": {
#                     "base": "2015=100",
#                     "eur_source": "ECB: ICP.M.U2.N.000000.4.INX (monthly -> annual avg)",
#                     "usd_source": "World Bank: FP.CPI.TOTL (rebased to 2015=100)",
#                     "include_partial_latest": bool(include_partial_latest),
#                     "last_updated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
#                     "range": {"start": start_year, "end": cand},
#                 },
#                 "EUR": {**cache.get("EUR", {}), **{str(y): float(v) for y, v in eur.items()}},
#                 "USD": {**cache.get("USD", {}), **{str(y): float(v) for y, v in usd.items()}},
#             }
#             path = save_cpi_cache(payload, dirpath, filename)
#             return path, cand
#         except Exception as e:
#             last_err = e
#             time.sleep(0.4)  # brief backoff and try previous year
#             continue

#     raise RuntimeError(f"Failed to update CPI cache for {targets}. Last error: {last_err}")

# #%% Check if update is needed

# def _today_ymd() -> str:
#     """Current date as 'DD/MM/YYYY' (UTC)."""
#     now = datetime.now(timezone.utc)
#     return now.strftime("%d/%m/%Y")

# def _current_ym() -> str:
#     """Current year-month as 'YYYY-MM'."""
#     now = datetime.now(timezone.utc)
#     return f"{now.year:04d}-{now.month:02d}"

# def _fetch_latest_month_ecb() -> str:
#     """Fetch latest published month for Euro area HICP (ECB), fast."""
#     url = "https://data-api.ecb.europa.eu/service/data/ICP/M.U2.N.000000.4.INX"
#     params = {"format": "csvdata", "lastNObservations": "1"}
#     r = requests.get(url, params=params, timeout=(4, 8))
#     r.raise_for_status()
#     rows = list(csv.reader(io.StringIO(r.text)))
#     header = rows[0]
#     ti = next((header.index(k) for k in ("TIME_PERIOD", "TIME", "Time") if k in header), None)
#     if ti is None:
#         raise RuntimeError("No TIME column in ECB response.")
#     return rows[-1][ti][:7]

# def _fetch_latest_month_oecd() -> str:
#     """Fetch latest published month for US CPI (OECD), fast."""
#     url = "https://stats.oecd.org/sdmx-json/data/PRICES_CPI/USA.CPALTT01.IXOB.M/all"
#     params = {"contentType": "csv", "lastNObservations": "1"}
#     r = requests.get(url, params=params, timeout=(4, 8))
#     r.raise_for_status()
#     rows = list(csv.reader(io.StringIO(r.text)))
#     header = rows[0]
#     ti = next((header.index(k) for k in ("TIME", "TIME_PERIOD", "Time") if k in header), None)
#     if ti is None:
#         raise RuntimeError("No TIME column in OECD response.")
#     return rows[-1][ti][:7]

# def check_and_update_if_needed(
#     start_year: int = 1990,
#     include_partial_latest: bool = True,
#     dirpath: str = CACHE_DIR,
#     filename: str = CACHE_FILE,
# ) -> dict:
#     """
#     Checks the last update date in the local cache.
#     If it's not from the current month, update CPI data (EUR+USD) automatically.
#     Returns a dict describing status and dates.
#     """
#     today = _today_ymd()
#     today_ym = _current_ym()

#     # Load existing metadata
#     cache = load_cpi_cache(dirpath, filename)
#     meta = cache.get("meta", {})

#     last_update_str = meta.get("last_update_date")
#     eur_latest = meta.get("eur_latest_month")
#     usd_latest = meta.get("usd_latest_month")

#     # Determine whether last update is from the current month
#     needs_update = True
#     if last_update_str:
#         try:
#             last_dt = datetime.strptime(last_update_str, "%d/%m/%Y")
#             if f"{last_dt.year:04d}-{last_dt.month:02d}" == today_ym:
#                 needs_update = False
#         except Exception:
#             pass  # if parsing fails, just force an update

#     if needs_update:
#         # Quick check from APIs to confirm new data may exist
#         eur_new = _fetch_latest_month_ecb()
#         usd_new = _fetch_latest_month_oecd()

#         update_cpi_cache(
#             start_year=start_year,
#             include_partial_latest=include_partial_latest,
#             dirpath=dirpath,
#             filename=filename,
#         )

#         # Reload cache and patch metadata
#         cache = load_cpi_cache(dirpath, filename)
#         meta = cache.get("meta", {})
#         meta.update({
#             "eur_latest_month": eur_new,
#             "usd_latest_month": usd_new,
#             "last_update_date": today,
#             "last_checked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
#         })
#         cache["meta"] = meta
#         save_cpi_cache(cache, dirpath, filename)

#         return {
#             "updated": True,
#             "last_update_date": today,
#             "eur_latest_month": eur_new,
#             "usd_latest_month": usd_new,
#         }

#     # No update needed
#     return {
#         "updated": False,
#         "last_update_date": last_update_str,
#         "eur_latest_month": eur_latest,
#         "usd_latest_month": usd_latest,
#     }

# #%% Get CPI Data

# def get_cpi_value(currency: str, year: int,
#                   dirpath: str = CACHE_DIR, filename: str = CACHE_FILE) -> float:
    
#     check_and_update_if_needed(start_year=1990, include_partial_latest=True)
    
#     """Read a single CPI value (2015=100) from the local cache."""
#     data = load_cpi_cache(dirpath, filename)
#     bucket = data.get(currency.upper(), {})
#     key = str(year)
#     if key not in bucket:
#         raise KeyError(f"{currency} CPI for {year} not in cache. Run update_cpi_cache() online first.")
#     return float(bucket[key])

# def get_cpi_series(dirpath: str = CACHE_DIR, filename: str = CACHE_FILE):
    
#     check_and_update_if_needed(start_year=1990, include_partial_latest=True)
    
#     """Return (eur_dict:int->float, usd_dict:int->float, meta_dict) from the cache."""
#     data = load_cpi_cache(dirpath, filename)
#     eur = {int(k): float(v) for k, v in data.get("EUR", {}).items()}
#     usd = {int(k): float(v) for k, v in data.get("USD", {}).items()}
#     return eur, usd, data.get("meta", {})

# #%% Actualize the price

# def actualize_price(amount: float, year: int, currency: str = "EUR", target_year: int | None = None) -> float:
#     """
#     Inflate `amount` from `year` to `target_year` (default: current year) using CPI (2015=100).
#     - If currency is 'EUR', uses Euro HICP.
#     - If currency is 'USD', uses US CPI.
#     - No FX conversion. Returns a float rounded to 2 decimals.
#     """
#     # Ensure cache is up-to-date (fast; your helpers already do lightweight checks)
#     check_and_update_if_needed(start_year=1990, include_partial_latest=True)

#     # Load CPI series from your cache
#     eur_series, usd_series, _ = get_cpi_series()

#     # Pick target year (default: current year UTC)
#     from datetime import datetime, timezone
#     if target_year is None:
#         target_year = datetime.now(timezone.utc).year

#     currency = currency.upper()
#     if currency not in {"EUR", "USD"}:
#         raise ValueError("currency must be 'EUR' or 'USD'")

#     # Helper: pick latest CPI year <= ref_year
#     def _cpi_at_or_before(series: dict[int, float], ref_year: int) -> tuple[int, float]:
#         yrs = [y for y in series.keys() if y <= ref_year]
#         if not yrs:
#             raise ValueError(f"No CPI available up to {ref_year}")
#         y_used = max(yrs)
#         return y_used, series[y_used]

#     # Choose the right CPI series
#     series = eur_series if currency == "EUR" else usd_series

#     # Compute inflation factor and return actualized amount
#     _, cpi_base = _cpi_at_or_before(series, year)
#     _, cpi_ref  = _cpi_at_or_before(series, target_year)
#     return round(float(amount * (cpi_ref / cpi_base)), 2)

# #%%

# # ==== example usage ====
# if __name__ == "__main__":
#     # path, used_end = update_cpi_cache(start_year=1990, include_partial_latest=True)
#     # print(f"Cache written to: {path} (up to {used_end})")

#     year_EU = 2023
#     year_US = 2023

#     # Offline lookups:
#     print(f"EUR {year_EU}:", get_cpi_value("EUR", year_EU))
#     print(f"USD {year_US}:", get_cpi_value("USD", year_US))

#     # eur, usd, meta = get_cpi_series()
#     # print("Meta:", meta)

#     # €2,500 in 2010 -> today's € (inflation only)
#     print(actualize_price(2500, 2010, "EUR"))  
        
#     # $3,000 in 2012 -> today's $ (inflation only)
#     print(actualize_price(3000, 2012, "USD"))

# -*- coding: utf-8 -*-
# CPI cache + updater + inflation helper
# Cache is ALWAYS stored in: labothappy/toolbox/economics/.macro_cache/cpi_2015eq100.json

import requests, csv, io, json, os, time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

# ======================================================
# Centralized cache path (always under labothappy/toolbox/economics)
# ======================================================
PROJECT_ROOT = Path(__file__).resolve().parent          # -> labothappy/toolbox/finance
ECON_PATH    = PROJECT_ROOT.parent / "economics"        # -> labothappy/toolbox/economics
CACHE_DIR    = ECON_PATH / ".macro_cache"               # -> labothappy/toolbox/economics/.macro_cache
CACHE_FILE   = "cpi_2015eq100.json"                     # single authoritative file

def _ensure_dir(p: Path): p.mkdir(parents=True, exist_ok=True)

def _cache_path(filename: str = CACHE_FILE) -> Path:
    _ensure_dir(CACHE_DIR)
    return CACHE_DIR / filename

def load_cpi_cache(filename: str = CACHE_FILE) -> dict:
    path = _cache_path(filename)
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def save_cpi_cache(payload: dict, filename: str = CACHE_FILE) -> str:
    path = _cache_path(filename)
    tmp  = path.with_suffix(".tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2, sort_keys=True)
    os.replace(tmp, path)   # atomic write
    return str(path)

# ======================================================
# 1) Fetch CPI data (EUR from ECB; USD from OECD)
# ======================================================

def euro_hicp_annual(start=1990, end=2025, partial_ok=True, min_months=1):
    """
    Euro area HICP (All items, 2015=100) monthly from ECB -> annual averages.
    Optionally includes partial last year.
    """
    url = "https://data-api.ecb.europa.eu/service/data/ICP/M.U2.N.000000.4.INX"
    params = {"format": "csvdata", "startPeriod": str(start), "endPeriod": str(end)}
    r = requests.get(url, params=params, timeout=15)
    r.raise_for_status()

    rows = list(csv.reader(io.StringIO(r.text)))
    header = rows[0]
    ti, vi = header.index("TIME_PERIOD"), header.index("OBS_VALUE")

    by_year = defaultdict(list)
    for row in rows[1:]:
        if len(row) <= vi: continue
        t, v = row[ti], row[vi]
        if not t or not v: continue
        y = int(t[:4])
        by_year[y].append(float(v))

    out = {}
    for y, vals in by_year.items():
        needed = 12 if not partial_ok else max(1, int(min_months))
        if len(vals) >= needed:
            out[y] = round(sum(vals)/len(vals), 2)
    return dict(sorted(out.items()))

def us_cpi_annual(start=1990, end=2025, rebase_to_2015=True, partial_ok=True, min_months=1):
    """
    U.S. CPI (All items, Index 2015=100) via OECD monthly -> annual averages.
    Returns {year: annual_avg_2015eq100}. Keeps a robust header parse.
    """
    url = "https://stats.oecd.org/sdmx-json/data/PRICES_CPI/USA.CPALTT01.IXOB.M/all"
    params = {"contentType": "csv"}

    r = requests.get(url, params=params, timeout=15)
    r.raise_for_status()

    rows = list(csv.reader(io.StringIO(r.text)))
    if not rows: return {}

    header = rows[0]
    time_keys  = ["TIME", "Time", "TIME_PERIOD"]
    value_keys = ["Value", "OBS_VALUE"]
    ti = next((header.index(k) for k in time_keys  if k in header), None)
    vi = next((header.index(k) for k in value_keys if k in header), None)
    if ti is None or vi is None:
        raise RuntimeError(f"OECD header missing time/value cols: {header}")

    by_year = defaultdict(list)
    for row in rows[1:]:
        if len(row) <= max(ti, vi): continue
        t, v = row[ti], row[vi]
        if not t or not v: continue
        try:
            y = int(t[:4]); val = float(v)
        except ValueError:
            continue
        if y < start or y > end: continue
        by_year[y].append(val)

    out = {}
    for y, vals in by_year.items():
        needed = 12 if not partial_ok else max(1, int(min_months))
        if len(vals) >= needed:
            out[y] = round(sum(vals)/len(vals), 2)

    # Optional rebase (OECD is already 2015=100, but this enforces exact 2015==100.00)
    if rebase_to_2015 and 2015 in out and out[2015] not in (0, None):
        factor = 100.0 / out[2015]
        out = {y: round(v * factor, 2) for y, v in out.items()}

    return dict(sorted(out.items()))

# ======================================================
# 2) Fast “latest month” probes (so we only update when needed)
# ======================================================

def _today_ymd() -> str:
    return datetime.now(timezone.utc).strftime("%d/%m/%Y")

def _current_ym() -> str:
    now = datetime.now(timezone.utc)
    return f"{now.year:04d}-{now.month:02d}"

def _fetch_latest_month_ecb() -> str:
    url = "https://data-api.ecb.europa.eu/service/data/ICP/M.U2.N.000000.4.INX"
    params = {"format": "csvdata", "lastNObservations": "1"}
    r = requests.get(url, params=params, timeout=(4,8))
    r.raise_for_status()
    rows = list(csv.reader(io.StringIO(r.text)))
    header = rows[0]
    ti = next((header.index(k) for k in ("TIME_PERIOD","TIME","Time") if k in header), None)
    if ti is None: raise RuntimeError("ECB response missing TIME column.")
    return rows[-1][ti][:7]  # 'YYYY-MM'

def _fetch_latest_month_oecd() -> str:
    url = "https://stats.oecd.org/sdmx-json/data/PRICES_CPI/USA.CPALTT01.IXOB.M/all"
    params = {"contentType": "csv", "lastNObservations": "1"}
    r = requests.get(url, params=params, timeout=(4,8))
    r.raise_for_status()
    rows = list(csv.reader(io.StringIO(r.text)))
    header = rows[0]
    ti = next((header.index(k) for k in ("TIME","TIME_PERIOD","Time") if k in header), None)
    if ti is None: raise RuntimeError("OECD response missing TIME column.")
    return rows[-1][ti][:7]  # 'YYYY-MM'

# ======================================================
# 3) Update/refresh cache (single canonical file)
# ======================================================

def update_cpi_cache(start_year: int = 1990, end_year: int | None = None,
                     include_partial_latest: bool = True) -> tuple[str, int]:
    """
    Create/update JSON cache with:
      - EUR: ECB HICP (2015=100) annual avg (monthly -> annual)
      - USD: OECD CPI (2015=100) annual avg (monthly -> annual)
    Tries current year; falls back to previous year if the fetch fails.
    Returns (cache_path, chosen_end_year).
    """
    cy = datetime.now(timezone.utc).year
    targets = [end_year or cy, (end_year or cy) - 1]

    last_err = None
    for cand in targets:
        try:
            eur = euro_hicp_annual(start=start_year, end=cand, partial_ok=include_partial_latest)
            usd = us_cpi_annual(start=start_year, end=cand, rebase_to_2015=True,
                                partial_ok=include_partial_latest)

            cache = load_cpi_cache()
            payload = {
                "meta": {
                    "base": "2015=100",
                    "eur_source": "ECB: ICP.M.U2.N.000000.4.INX (monthly -> annual avg)",
                    "usd_source": "OECD: USA.CPALTT01.IXOB.M (monthly -> annual avg, 2015=100)",
                    "include_partial_latest": bool(include_partial_latest),
                    "last_updated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                    "range": {"start": start_year, "end": cand},
                },
                "EUR": {**cache.get("EUR", {}), **{str(y): float(v) for y, v in eur.items()}},
                "USD": {**cache.get("USD", {}), **{str(y): float(v) for y, v in usd.items()}},
            }
            path = save_cpi_cache(payload)
            return path, cand
        except Exception as e:
            last_err = e
            time.sleep(0.4)
            continue

    raise RuntimeError(f"Failed to update CPI cache for {targets}. Last error: {last_err}")

def check_and_update_if_needed(start_year: int = 1990,
                               include_partial_latest: bool = True) -> dict:
    """
    Reads cache meta and last update date.
    If not updated this month, probes sources (fast) and refreshes the cache.
    """
    today = _today_ymd()
    today_ym = _current_ym()

    cache = load_cpi_cache()
    meta  = cache.get("meta", {})

    last_update_str = meta.get("last_update_date")
    eur_latest = meta.get("eur_latest_month")
    usd_latest = meta.get("usd_latest_month")

    needs_update = True
    if last_update_str:
        try:
            last_dt = datetime.strptime(last_update_str, "%d/%m/%Y")
            if f"{last_dt.year:04d}-{last_dt.month:02d}" == today_ym:
                needs_update = False
        except Exception:
            pass

    if needs_update:
        eur_new = _fetch_latest_month_ecb()
        usd_new = _fetch_latest_month_oecd()
        update_cpi_cache(start_year=start_year, include_partial_latest=include_partial_latest)
        # reload & stamp meta
        cache = load_cpi_cache()
        meta = cache.get("meta", {})
        meta.update({
            "eur_latest_month": eur_new,
            "usd_latest_month": usd_new,
            "last_update_date": today,
            "last_checked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        })
        cache["meta"] = meta
        save_cpi_cache(cache)
        return {"updated": True, "last_update_date": today,
                "eur_latest_month": eur_new, "usd_latest_month": usd_new}

    return {"updated": False, "last_update_date": last_update_str,
            "eur_latest_month": eur_latest, "usd_latest_month": usd_latest}

# ======================================================
# 4) Public getters / inflation helper
# ======================================================

def get_cpi_value(currency: str, year: int) -> float:
    """Read a single CPI value (2015=100) from the canonical cache (auto-refresh monthly)."""
    check_and_update_if_needed(start_year=1990, include_partial_latest=True)
    data = load_cpi_cache()
    bucket = data.get(currency.upper(), {})
    key = str(year)
    if key not in bucket:
        raise KeyError(f"{currency} CPI for {year} not in cache. Run update_cpi_cache() online first.")
    return float(bucket[key])

def get_cpi_series():
    """Return (eur_dict:int->float, usd_dict:int->float, meta_dict) from the canonical cache."""
    check_and_update_if_needed(start_year=1990, include_partial_latest=True)
    data = load_cpi_cache()
    eur = {int(k): float(v) for k, v in data.get("EUR", {}).items()}
    usd = {int(k): float(v) for k, v in data.get("USD", {}).items()}
    return eur, usd, data.get("meta", {})

def actualize_price(amount: float, year: int, currency: str = "EUR", target_year: int | None = None) -> float:
    """
    Inflate `amount` from `year` to `target_year` (default: current year UTC) using CPI (2015=100).
    - If currency is 'EUR', uses Euro HICP.
    - If currency is 'USD', uses US CPI.
    - No FX conversion. Returns rounded float.
    """
    check_and_update_if_needed(start_year=1990, include_partial_latest=True)
    eur_series, usd_series, _ = get_cpi_series()

    if target_year is None:
        target_year = datetime.now(timezone.utc).year

    currency = currency.upper()
    if currency not in {"EUR", "USD"}:
        raise ValueError("currency must be 'EUR' or 'USD'")

    def _cpi_at_or_before(series: dict[int, float], ref_year: int) -> tuple[int, float]:
        yrs = [y for y in series if y <= ref_year]
        if not yrs: raise ValueError(f"No CPI available up to {ref_year}")
        y_used = max(yrs)
        return y_used, series[y_used]

    series = eur_series if currency == "EUR" else usd_series
    _, cpi_base = _cpi_at_or_before(series, year)
    _, cpi_ref  = _cpi_at_or_before(series, target_year)
    
    return round(float(amount * (cpi_ref / cpi_base)), 2)

# ======================================================
# Example usage (optional)
# ======================================================
if __name__ == "__main__":
    # Force an update (optional)
    # path, used_end = update_cpi_cache(start_year=1990, include_partial_latest=True)
    # print(f"Cache: {path} up to {used_end}")

    y_eu = 2023
    y_us = 2023
    print(f"EUR {y_eu}:", get_cpi_value("EUR", y_eu))
    print(f"USD {y_us}:", get_cpi_value("USD", y_us))
    print("€2,500 (2010) → today €:", actualize_price(2500, 2010, "EUR"))
    print("$3,000 (2012) → today $:", actualize_price(3000, 2012, "USD"))

