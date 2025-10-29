# ========= EUR/USD yearly average cache + updater (ECB, no keys) =========
import requests, csv, io, json, os, time
from datetime import datetime, timezone

FX_CACHE_DIR  = ".macro_cache"
FX_CACHE_FILE = "eurusd.json"   # {"meta": {...}, "EURUSD": {"1999": 1.07, ...}}

def _fx_cache_path() -> str:
    os.makedirs(FX_CACHE_DIR, exist_ok=True)
    return os.path.join(FX_CACHE_DIR, FX_CACHE_FILE)

def _fx_load() -> dict:
    p = _fx_cache_path()
    if not os.path.exists(p):
        return {}
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)

def _fx_save(obj: dict) -> str:
    p = _fx_cache_path()
    with open(p, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2, sort_keys=True)
    return p

def _today_str() -> str:
    return datetime.now(timezone.utc).strftime("%d/%m/%Y")

def _current_ym() -> str:
    now = datetime.now(timezone.utc)
    return f"{now.year:04d}-{now.month:02d}"

# ---- FAST freshness checks: fetch just the latest month ----
def _ecb_fx_latest_month(timeout=(4, 8)) -> str:
    """
    Latest month published for EURUSD (ECB EXR, daily spot).
    Series: EXR.D.USD.EUR.SP00.A (USD per 1 EUR), daily.
    Returns 'YYYY-MM'.
    """
    url = "https://data-api.ecb.europa.eu/service/data/EXR/D.USD.EUR.SP00.A"
    params = {"format": "csvdata", "lastNObservations": "1", "detail": "dataonly"}
    r = requests.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    rows = list(csv.reader(io.StringIO(r.text)))
    header = rows[0]
    ti = next((header.index(k) for k in ("TIME_PERIOD","TIME","Time") if k in header), None)
    if ti is None:
        raise RuntimeError(f"ECB FX header missing TIME column: {header}")
    latest = rows[-1][ti].strip()
    return latest[:7]

# ---- Yearly (or YTD) average from ECB daily spot ----
def _ecb_fx_year_avg(year: int, timeout=(6, 15)) -> float:
    """
    Average EURUSD (USD per 1 EUR) for a given year from ECB.
    If current year, this returns the YTD average (available days).
    """
    url = "https://data-api.ecb.europa.eu/service/data/EXR/D.USD.EUR.SP00.A"
    params = {
        "startPeriod": f"{year}-01-01",
        "endPeriod":   f"{year}-12-31",
        "format": "csvdata",
        "detail": "dataonly",
    }
    r = requests.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    rows = list(csv.reader(io.StringIO(r.text)))
    header = rows[0]
    if "OBS_VALUE" not in header:
        raise RuntimeError("ECB CSV missing OBS_VALUE")
    vi = header.index("OBS_VALUE")

    total = 0.0; n = 0
    for row in rows[1:]:
        if len(row) <= vi: continue
        val = row[vi].strip()
        if not val or val.upper() == "NA": continue
        total += float(val); n += 1
    if n == 0:
        raise ValueError(f"No EUR/USD observations for {year}")
    return total / n

# ---- Update cache up to current (or previous) year, store meta ----
def update_eurusd_cache(start_year: int = 1999, end_year: int | None = None) -> tuple[str,int]:
    """
    Builds/updates ./macro_cache/eurusd.json with yearly (or YTD) averages.
    Tries current year; on failure, falls back to current_year-1.
    """
    cy = datetime.now(timezone.utc).year
    targets = [end_year or cy, (end_year or cy) - 1]

    cache = _fx_load()
    series = cache.get("EURUSD", {})
    meta   = cache.get("meta", {})

    last_err = None
    for cand in targets:
        try:
            # compute (or refresh) each year from start_year..cand
            changed = False
            for y in range(start_year, cand + 1):
                try:
                    avg = _ecb_fx_year_avg(y)
                    s = round(float(avg), 6)  # 6 dp for FX cleanliness
                    if series.get(str(y)) != s:
                        series[str(y)] = s
                        changed = True
                except Exception as e:
                    # allow gaps for the earliest years if ECB doesn't have them
                    if y < 1999:
                        continue
                    else:
                        raise

            latest_month = _ecb_fx_latest_month()
            meta.update({
                "eurusd_source": "ECB EXR D.USD.EUR.SP00.A (daily spot), yearly avg/YTD",
                "last_update_date": _today_str(),
                "last_checked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                "eurusd_latest_month": latest_month,
                "range": {"start": start_year, "end": cand},
            })
            cache["EURUSD"] = series
            cache["meta"]   = meta
            path = _fx_save(cache)
            return path, cand
        except Exception as e:
            last_err = e
            time.sleep(0.4)
            continue
    raise RuntimeError(f"Failed to update EURUSD cache for {targets}. Last error: {last_err}")

# ---- Fast check: if last_update_date not current month, refresh ----
def ensure_eurusd_cache_current_month(start_year: int = 1999) -> dict:
    """
    If eurusd.json wasn't updated in the current month, refresh it.
    Returns status dict.
    """
    cache = _fx_load()
    meta = cache.get("meta", {})
    last_update = meta.get("last_update_date")  # "DD/MM/YYYY"
    today_ym = _current_ym()

    needs_update = True
    if last_update:
        try:
            dt = datetime.strptime(last_update, "%d/%m/%Y")
            if f"{dt.year:04d}-{dt.month:02d}" == today_ym:
                needs_update = False
        except Exception:
            needs_update = True

    if needs_update:
        path, used_end = update_eurusd_cache(start_year=start_year)
        cache = _fx_load()
        meta = cache.get("meta", {})
        return {"updated": True, "path": path, "used_end_year": used_end, "meta": meta}
    else:
        return {"updated": False, "meta": meta}

# ========= Replacement euro_dollar() that uses the cache first =========
def _eurusd_year_avg_cached_or_fetch(year: int) -> float:
    """
    Return cached yearly (or YTD) EURUSD average. If missing, fetch + persist.
    """
    cache = _fx_load()
    val = cache.get("EURUSD", {}).get(str(year))
    if val is not None:
        return float(val)
    # Not cached: update just enough (will compute YTD for current year)
    update_eurusd_cache(start_year=min(1999, year), end_year=max(year, datetime.now(timezone.utc).year))
    cache = _fx_load()
    if str(year) not in cache.get("EURUSD", {}):
        raise KeyError(f"EURUSD average for {year} not available after update.")
    return float(cache["EURUSD"][str(year)])

def euro_dollar(amount, year, input_currency="USD"):
    """
    Convert between USD and EUR using the cached average EUR/USD rate for a given year.
    Falls back to ECB fetch when cache is missing and persists results in .macro_cache/eurusd.json.

    Parameters:
        amount (float): The amount in the input currency.
        year (int): The year to use for the average EUR/USD rate.
        input_currency (str): "USD" (convert USD → EUR) or "EUR" (convert EUR → USD).

    Returns:
        float: The converted amount.
    """
    # ensure cache is fresh (only once per month)
    try:
        ensure_eurusd_cache_current_month()
    except Exception:
        # Non-fatal; conversion can still occur using existing cache or on-demand fetch
        pass

    avg_rate = _eurusd_year_avg_cached_or_fetch(year)  # USD per 1 EUR
    if input_currency.upper() == "USD":
        return float(amount / avg_rate)  # USD → EUR
    elif input_currency.upper() == "EUR":
        return float(amount * avg_rate)  # EUR → USD
    else:
        raise ValueError("input_currency must be 'USD' or 'EUR'")


if __name__ == "__main__":
    # Convert 1000 USD to EUR in 2020
    print(euro_dollar(1000, 2020, "USD"))
    
    # Convert 1000 EUR to USD in 2020
    print(euro_dollar(1000, 2020, "EUR"))
    