#!/usr/bin/env python3
"""
seestar_plan.py: observability planner for a smart telescope, with a custom (tree) horizon.

Two commands
------------
1) season: assign every target to its single BEST month window and write, per window,
   (a window = the calendar month + the first --overlap=15 days of the next month, so the
   second half of each file overlaps the start of the next month's file)
     YYYY-MM.csv: every night of that month for each of its targets
                  (best time, usable window, max altitude, Moon), and
     YYYY-MM.png: for those targets vs day of month: altitude at --ref-time (default 21:00),
                  best time of night, and Moon separation; Moon phase shaded.
   plus season_summary.csv (one row per target) and not_observable.txt.

     python seestar_plan.py season 2027-01 --months 12 --horizon horizon.csv

2) staralt: Staralt-style plot (altitude vs local time) for one night and one or more
   targets, with twilight bands, Moon, airmass axis and your horizon, plus a table.

     python seestar_plan.py staralt 2027-01-12 Rosette
     python seestar_plan.py staralt 2027-01-12 Rosette M42 Cone --horizon horizon.csv
     python seestar_plan.py staralt 2027-01-12                      # top 8 targets that night
     python seestar_plan.py staralt 2027-01-12 --coord "Thor,07:18:30,-13:13:48"
     python seestar_plan.py staralt 2027-01-12 "NGC 2174" --resolve  # look up name online (Sesame)

What "best month" means
-----------------------
For each night we take the observing SESSION = astronomical darkness (sun < --twilight)
up to --session-end (default 03:00 local; use "dawn" for all-night observing).
A target's nightly score = sum over session time when it is above your horizon of
sin(altitude), i.e. usable hours weighted by 1/airmass. Best month = highest mean score.
This favours months when the target is high during your evening session, so M81/M82 lands
in spring rather than in December just because December nights are long.

Common options
  --targets targets.csv   name,ra,dec[,type,members] (J2000). Name matching is case-insensitive
                          substring, so "rosette" finds Rosette_NGC2237.
  --horizon horizon.csv   az,alt profile of local obstructions (optional)
  --min-alt 15            altitude floor everywhere (deg); combined with horizon as max()
  --twilight -18          sun altitude defining darkness (-12 = nautical)
  --lat/--lon/--elev/--utc-offset   default Tempe, AZ (MST = UTC-7, no DST)

Method notes
  Target positions: J2000 precessed to date; alt/az from local sidereal time with standard
  spherical-astronomy formulas (no refraction; good to ~0.1 deg). Sun and Moon positions and
  Moon illumination from astropy (built-in ERFA ephemeris). Twilight from the sun's altitude.
"""
import argparse
import calendar
import os
import sys
import warnings
from datetime import date, timedelta

import numpy as np

warnings.filterwarnings("ignore")
from astropy import units as u  # noqa: E402
from astropy.coordinates import FK5, AltAz, EarthLocation, SkyCoord, get_body, get_sun  # noqa: E402
from astropy.time import Time  # noqa: E402
from astropy.utils import iers  # noqa: E402

iers.conf.auto_download = False  # offline is fine; UT1 precision is irrelevant for planning

DEF_LAT, DEF_LON, DEF_ELEV, DEF_UTC = 33.3981, -111.9069, 360.0, -7.0  # McClintock High School, Tempe AZ
GRID_START_LOCAL = 17.0   # grid runs 17:00 local -> 07:00 next morning
GRID_HOURS = 14.0
SIDEREAL = 1.00273790935


# =============================================================================== inputs
class Catalog:
    def __init__(self, names, coords, types):
        self.names, self.coords, self.types = names, coords, types

    def subset(self, idx):
        return Catalog([self.names[i] for i in idx], self.coords[idx], [self.types[i] for i in idx])


def read_targets(path):
    names, ras, decs, types = [], [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.lower().startswith("name,"):
                continue
            p = [x.strip() for x in line.split(",")]
            names.append(p[0]); ras.append(p[1]); decs.append(p[2])
            types.append(p[3] if len(p) > 3 else "")
    return Catalog(names, _coords(ras, decs), types)


def _coords(ras, decs):
    sexa = [(":" in r or " " in r) for r in ras]
    out = [SkyCoord(r, d, unit=(u.hourangle, u.deg) if s else (u.deg, u.deg)) for r, d, s in zip(ras, decs, sexa)]
    return SkyCoord(out)


def select(cat, queries, extra_coords=(), resolve=False):
    """Pick targets by case-insensitive substring; add --coord entries; optionally resolve unknown names."""
    names, coords, types = [], [], []
    for q in queries:
        hits = [i for i, n in enumerate(cat.names) if q.lower().replace(" ", "") in n.lower()]
        if hits:
            for i in hits:
                if cat.names[i] not in names:
                    names.append(cat.names[i]); coords.append(cat.coords[i]); types.append(cat.types[i])
        elif resolve:
            try:
                c = SkyCoord.from_name(q)
            except Exception as e:
                sys.exit(f"could not resolve {q!r}: {e}")
            names.append(q.replace(" ", "")); coords.append(c); types.append("resolved")
            print(f"resolved {q!r} -> {c.to_string('hmsdms', precision=0)}")
        else:
            sys.exit(f"{q!r} not in {len(cat.names)}-target catalog. Use --resolve to look it up online, "
                     f"or --coord 'name,ra,dec'.")
    for spec in extra_coords:
        n, r, d = [x.strip() for x in spec.split(",")[:3]]
        names.append(n); coords.append(_coords([r], [d])[0]); types.append("user")
    return Catalog(names, SkyCoord(coords), types)


def read_horizon(path, min_alt):
    if not path:
        return lambda az: np.full(np.shape(az), float(min_alt))
    pts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.lower().startswith("az"):
                continue
            a, h = (float(x) for x in line.split(",")[:2])
            pts.append((a % 360.0, h))
    pts.sort()
    az = np.array([p[0] for p in pts]); alt = np.array([p[1] for p in pts])
    azp = np.concatenate([az - 360, az, az + 360]); altp = np.tile(alt, 3)
    return lambda a: np.maximum(np.interp(np.mod(a, 360.0), azp, altp), min_alt)


# =============================================================================== astronomy
def night_grid(dates, utc_off, step_min):
    n = int(round(GRID_HOURS * 60 / step_min)) + 1
    offs = np.arange(n) * step_min / 60.0  # hours after 17:00 local
    base = Time([f"{d}T17:00:00" for d in dates], scale="utc") - utc_off * u.hour
    return Time(base.jd[:, None] + offs[None, :] / 24.0, format="jd", scale="utc"), offs


def target_altaz(coords, times, lat, lon):
    """alt, az with shape (n_targets, *times.shape)."""
    mid = Time(np.mean(times.jd), format="jd")
    c = coords.transform_to(FK5(equinox=mid))
    lst = times.sidereal_time("mean", longitude=lon * u.deg).rad
    shp = (-1,) + (1,) * lst.ndim
    ra = c.ra.rad.reshape(shp); d = c.dec.rad.reshape(shp)
    phi = np.deg2rad(lat)
    H = lst[None] - ra
    alt = np.rad2deg(np.arcsin(np.clip(np.sin(phi) * np.sin(d) + np.cos(phi) * np.cos(d) * np.cos(H), -1, 1)))
    az = np.rad2deg(np.arctan2(-np.cos(d) * np.sin(H), np.sin(d) * np.cos(phi) - np.cos(d) * np.cos(H) * np.sin(phi)))
    # local clock time of meridian transit for each night (hours after 17:00), from LST at grid start
    lst0 = lst.reshape(lst.shape[0], -1)[:, 0] if lst.ndim > 1 else np.atleast_1d(lst[0])
    transit = np.mod(c.ra.rad[:, None] - lst0[None, :], 2 * np.pi) * 12 / np.pi / SIDEREAL
    return alt, np.mod(az, 360.0), transit


def sun_moon(times, loc):
    flat = times.ravel()
    frame = AltAz(obstime=flat, location=loc)
    sun = get_sun(flat)
    moon = get_body("moon", flat, loc)
    maa = moon.transform_to(frame)
    illum = (1 - np.cos(sun.separation(moon).rad)) / 2
    r = lambda x: np.asarray(x).reshape(times.shape)  # noqa: E731
    return r(sun.transform_to(frame).alt.deg), r(maa.alt.deg), r(maa.az.deg), r(illum)


def angsep(alt1, az1, alt2, az2):
    a1, z1, a2, z2 = map(np.deg2rad, (alt1, az1, alt2, az2))
    return np.rad2deg(np.arccos(np.clip(np.sin(a1) * np.sin(a2) + np.cos(a1) * np.cos(a2) * np.cos(z1 - z2), -1, 1)))


def clock(h):
    if h is None or not np.isfinite(h):
        return ""
    t = (GRID_START_LOCAL + h) % 24
    m = int(round(t * 60)) % (24 * 60)
    return f"{m // 60:02d}:{m % 60:02d}"


def compass(az):
    pts = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]
    return pts[int((az % 360) / 22.5 + 0.5) % 16]


def parse_session_end(s):
    if s.lower() == "dawn":
        return GRID_HOURS
    hh, mm = (int(x) for x in s.split(":"))
    return ((hh + mm / 60 - GRID_START_LOCAL) % 24)


def nightly(cat, dates, a, loc, limit_fn, session_end):
    """Compute per-(target, night) quantities on the time grid."""
    times, offs = night_grid(dates, a.utc_offset, a.step_min)
    alt, az, transit = target_altaz(cat.coords, times, a.lat, a.lon)
    sun_alt, moon_alt, moon_az, illum = sun_moon(times, loc)
    dark = sun_alt < a.twilight                                   # (N,S)
    session = dark & (offs[None, :] <= session_end + 1e-9)        # (N,S)
    lim = limit_fn(az)                                            # (T,N,S)
    above = alt > lim
    return dict(times=times, offs=offs, alt=alt, az=az, lim=lim, transit=transit, dark=dark, session=session,
                above=above, sun_alt=sun_alt, moon_alt=moon_alt, moon_az=moon_az, illum=illum)


def night_stats(R, t, n, step_h):
    """Summary for target t on night n."""
    offs = R["offs"]; alt = R["alt"][t, n]; az = R["az"][t, n]
    ses_ok = R["above"][t, n] & R["session"][n]
    all_ok = R["above"][t, n] & R["dark"][n]
    idx = np.where(ses_ok)[0]
    pool = idx if idx.size else np.where(R["dark"][n])[0]
    best = pool[np.argmax(alt[pool])] if pool.size else int(np.argmax(alt))
    dk = np.where(R["dark"][n])[0]
    return dict(
        transit=R["transit"][t, n],
        usable_session_h=idx.size * step_h,
        usable_night_h=all_ok.sum() * step_h,
        start=offs[idx[0]] if idx.size else np.nan,
        end=offs[idx[-1]] if idx.size else np.nan,
        best=offs[best], max_alt=alt[best], az=az[best], lim=R["lim"][t, n, best],
        airmass=1 / np.sin(np.deg2rad(alt[best])) if alt[best] > 1 else np.inf,
        score=(np.sin(np.deg2rad(np.clip(alt, 0, 90))) * ses_ok).sum() * step_h,
        moon_illum=R["illum"][n, best], moon_alt=R["moon_alt"][n, best],
        moon_sep=angsep(alt[best], az[best], R["moon_alt"][n, best], R["moon_az"][n, best]),
        dark_start=offs[dk[0]] if dk.size else np.nan, dark_end=offs[dk[-1]] if dk.size else np.nan,
    )


# =============================================================================== season
def run_season(a, cat, loc, limit_fn):
    y, m = (int(x) for x in a.start.split("-"))
    months = []
    for k in range(a.months):
        yy, mm = y + (m - 1 + k) // 12, (m - 1 + k) % 12 + 1
        months.append((yy, mm))
    # each "month" window = the calendar month + the first --overlap days of the next month
    windows = []
    for yy, mm in months:
        d0 = date(yy, mm, 1)
        d1 = date(yy, mm, calendar.monthrange(yy, mm)[1]) + timedelta(days=a.overlap)
        windows.append((f"{yy}-{mm:02d}", d0, d1))
    first, last = windows[0][1], windows[-1][2]
    dates = [(first + timedelta(days=i)).isoformat() for i in range((last - first).days + 1)]
    dobj = np.array([first + timedelta(days=i) for i in range(len(dates))])
    win_idx = {k: np.where((dobj >= d0) & (dobj <= d1))[0] for k, d0, d1 in windows}
    session_end = parse_session_end(a.session_end)
    step_h = a.step_min / 60
    print(f"computing {len(cat.names)} targets x {len(dates)} nights ...", file=sys.stderr)
    R = nightly(cat, dates, a, loc, limit_fn, session_end)

    T, N = len(cat.names), len(dates)
    S = np.zeros((T, N)); stats = [[None] * N for _ in range(T)]
    for t in range(T):
        for n in range(N):
            s = night_stats(R, t, n, step_h); stats[t][n] = s; S[t, n] = s["score"]

    mkeys = [w[0] for w in windows]
    mean_score = np.array([[S[t, win_idx[k]].mean() for k in mkeys] for t in range(T)])
    best_month = {t: mkeys[int(np.argmax(mean_score[t]))] for t in range(T) if mean_score[t].max() > 0}
    never = [cat.names[t] for t in range(T) if t not in best_month]

    os.makedirs(a.out, exist_ok=True)
    cols = ["date", "target", "type", "transit", "best_time", "max_alt", "airmass", "direction",
            "horizon_there", "usable_start", "usable_end", "usable_h_session", "usable_h_night",
            "moon_illum_pct", "moon_up", "moon_sep_deg", "dark_start", "dark_end"]
    summary = []
    for k in mkeys:
        nights = win_idx[k]
        ts = sorted([t for t, bm in best_month.items() if bm == k],
                    key=lambda t: np.median([stats[t][n]["best"] for n in nights]))
        with open(os.path.join(a.out, f"{k}.csv"), "w") as f:
            f.write(",".join(cols) + "\n")
            for t in ts:
                for n in nights:
                    s = stats[t][n]
                    f.write(",".join(str(x) for x in [
                        dates[n], cat.names[t], cat.types[t], clock(s["transit"]), clock(s["best"]),
                        f"{s['max_alt']:.1f}", f"{s['airmass']:.2f}", compass(s["az"]), f"{s['lim']:.0f}",
                        clock(s["start"]), clock(s["end"]), f"{s['usable_session_h']:.2f}",
                        f"{s['usable_night_h']:.2f}", f"{100 * s['moon_illum']:.0f}", int(s["moon_alt"] > 0),
                        f"{s['moon_sep']:.0f}", clock(s["dark_start"]), clock(s["dark_end"])]) + "\n")
                mid = nights[len(nights) // 2]
                summary.append((k, cat.names[t], cat.types[t], clock(stats[t][mid]["best"]),
                                max(stats[t][n]["max_alt"] for n in nights),
                                np.mean([stats[t][n]["usable_session_h"] for n in nights])))
        plot_month(a, k, ts, nights, dates, cat, stats, R)

    with open(os.path.join(a.out, "season_summary.csv"), "w") as f:
        f.write("best_month,target,type,best_time_midmonth,max_alt,mean_usable_h_session\n")
        for r in summary:
            f.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[4]:.1f},{r[5]:.2f}\n")
    with open(os.path.join(a.out, "not_observable.txt"), "w") as f:
        f.write("Never above your horizon during the session in this period:\n" + "\n".join(never) + "\n")

    print(f"\nBest window per target (each month + {a.overlap} days of the next; "
          f"session = dark until {a.session_end}, sun < {a.twilight} deg)\n")
    for k in mkeys:
        rows = [r for r in summary if r[0] == k]
        print(f"{calendar.month_name[int(k[5:])]:>9s} {k[:4]}: " +
              (", ".join(f"{r[1]} ({r[3]})" for r in rows) if rows else "(none)"))
    if never:
        print(f"\nnever usable: {', '.join(never)}")
    print(f"\nwrote {a.out}/YYYY-MM.csv, YYYY-MM.png, season_summary.csv, not_observable.txt")


def plot_month(a, key, ts, nights, dates, cat, stats, R):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    yy, mm = int(key[:4]), int(key[5:])
    days = np.arange(len(nights))  # x = night index within the window
    labels = [dates[n][5:] for n in nights]  # MM-DD
    ref = parse_session_end(a.ref_time)
    iref = int(np.argmin(np.abs(R["offs"] - ref)))
    fig, axs = plt.subplots(3, 1, figsize=(12, 10), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.6, 1.4]))
    ax, axb, axm = axs
    # Moon illumination as background shading (darker = brighter Moon)
    mi = np.array([R["illum"][n][np.argmin(np.abs(R["offs"] - 7))] for n in nights])  # at local midnight
    for d, ill in zip(days, mi):
        for axx in axs:
            axx.axvspan(d - 0.5, d + 0.5, color="0.55", alpha=0.35 * ill, lw=0)
    cmap = plt.get_cmap("tab10" if len(ts) <= 10 else "tab20")
    for k, t in enumerate(ts):
        col = cmap(k % cmap.N)
        y = np.array([R["alt"][t, n, iref] for n in nights])
        lim = np.array([R["lim"][t, n, iref] for n in nights])
        mx = max(stats[t][n]["max_alt"] for n in nights)
        ax.plot(days, np.where(y > 0, y, np.nan), "-", color=col, lw=1, alpha=0.4)
        ax.plot(days, np.where(y > lim, y, np.nan), "-o", ms=3, color=col, lw=2,
                label=f"{cat.names[t]}  (peaks at {mx:.0f}\u00b0)")
        bt = np.array([stats[t][n]["best"] for n in nights]) + GRID_START_LOCAL
        axb.plot(days, bt, "-", color=col, lw=2)
        axm.plot(days, [stats[t][n]["moon_sep"] for n in nights], "-", color=col, lw=1.6)
    marks = [(i, "new Moon") for i in range(len(mi)) if mi[i] == mi[max(0, i - 1):i + 2].min() and mi[i] < 0.1]
    marks += [(i, "full Moon") for i in range(len(mi)) if mi[i] == mi[max(0, i - 1):i + 2].max() and mi[i] > 0.9]
    for i, lab in marks:
        ax.annotate(lab, (days[i], 1.0), xycoords=("data", "axes fraction"), xytext=(0, 3),
                    textcoords="offset points", ha="center", va="bottom", fontsize=8)
    if a.horizon is None:
        ax.axhline(a.min_alt, color="k", ls="--", lw=0.8, alpha=0.5)
    ax.set_ylabel(f"altitude at {a.ref_time} local (deg)")
    ax.set_ylim(0, 90); ax.grid(alpha=0.3)
    ax.set_title(f"{calendar.month_name[mm]} {yy} (+{a.overlap} days): targets at their best in this window\n"
                 f"bold = above your horizon at {a.ref_time}; grey background = Moon brightness", pad=16, fontsize=11)
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
    sess_end = GRID_START_LOCAL + parse_session_end(a.session_end)
    yt = np.arange(18, int(np.ceil(sess_end)) + 1, 1)
    axb.set_yticks(yt); axb.set_yticklabels([f"{h % 24:02d}:00" for h in yt])
    axb.set_ylim(18, sess_end + 0.5); axb.invert_yaxis(); axb.grid(alpha=0.3)
    axb.set_ylabel(f"best time\n(highest before {a.session_end})")
    axm.set_ylabel("Moon separation\nat best time (deg)"); axm.set_ylim(0, 180); axm.grid(alpha=0.3)
    axm.set_xlabel("evening date (MM-DD)")
    nxt = [i for i, n in enumerate(nights) if int(dates[n][5:7]) != mm]
    if nxt:
        for axx in axs:
            axx.axvline(nxt[0] - 0.5, color="k", lw=1, ls="-", alpha=0.6)
        ax.text(nxt[0] - 0.3, 2, "next month →", fontsize=8, va="bottom")
    tick = days[::2]
    axm.set_xlim(-0.5, days[-1] + 0.5); axm.set_xticks(tick)
    axm.set_xticklabels([labels[i] for i in tick], rotation=90, fontsize=8)
    if not ts:
        ax.text(0.5, 0.5, "no targets assigned to this month", transform=ax.transAxes, ha="center")
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, f"{key}.png"), dpi=120); plt.close(fig)


# =============================================================================== staralt
def run_staralt(a, cat, loc, limit_fn):
    session_end = parse_session_end(a.session_end)
    step_h = a.step_min / 60
    chosen = select(cat, a.targets, a.coord, a.resolve) if (a.targets or a.coord) else cat
    R = nightly(chosen, [a.date], a, loc, limit_fn, session_end)
    stats = [night_stats(R, t, 0, step_h) for t in range(len(chosen.names))]
    order = sorted(range(len(stats)), key=lambda t: (-stats[t]["usable_night_h"], -stats[t]["max_alt"]))
    if not (a.targets or a.coord):  # auto mode: top N usable targets
        order = [t for t in order if stats[t]["usable_night_h"] > 0][: a.top]

    s0 = stats[order[0]] if order else night_stats(R, 0, 0, step_h)
    moon_mid = R["illum"][0][np.argmin(np.abs(R["offs"] - 7))]
    print(f"\n{a.date} (evening)  dark {clock(s0['dark_start'])}-{clock(s0['dark_end'])}  "
          f"Moon {100 * moon_mid:.0f}% lit\n")
    hdr = (f"{'target':24s} {'transit':>7s} {'best':>5s} {'max alt':>7s} {'airm':>5s} {'dir':>3s} "
           f"{'trees':>5s} {'usable (all night)':>18s} {'moon sep':>8s}")
    print(hdr); print("-" * len(hdr))
    for t in order:
        s = stats[t]
        ok = R["above"][t, 0] & R["dark"][0]; idx = np.where(ok)[0]
        win = f"{clock(R['offs'][idx[0]])}-{clock(R['offs'][idx[-1]])}" if idx.size else "-"
        print(f"{chosen.names[t]:24s} {clock(s['transit']):>7s} {clock(s['best']):>5s} {s['max_alt']:6.1f}° "
              f"{s['airmass']:5.2f} {compass(s['az']):>3s} {s['lim']:4.0f}° {win:>12s} {s['usable_night_h']:4.1f}h "
              f"{s['moon_sep']:7.0f}°")
    plot_staralt(a, chosen, order, R, stats, limit_fn)


def plot_staralt(a, cat, order, R, stats, limit_fn):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    offs = R["offs"]; sun = R["sun_alt"][0]
    fig, ax = plt.subplots(figsize=(12, 6.5))
    # twilight bands
    bands = [(0, -6, "0.80"), (-6, -12, "0.86"), (-12, -18, "0.93")]
    ax.fill_between(offs, 0, 90, where=sun > 0, color="0.70", lw=0, step="mid", label="day")
    for hi, lo, c in bands:
        ax.fill_between(offs, 0, 90, where=(sun <= hi) & (sun > lo), color=c, lw=0, step="mid")
    for lvl, lab in ((-6, "civil"), (-12, "nautical"), (-18, "astronomical")):
        for i in np.where(np.diff(np.sign(sun - lvl)) != 0)[0]:
            ax.axvline(offs[i], color="0.6", lw=0.6, ls="--")
            ax.text(offs[i], 88, lab, rotation=90, va="top", ha="right", fontsize=7, color="0.4")
    # Moon
    ill = R["illum"][0][np.argmin(np.abs(offs - 7))]
    ax.plot(offs, np.where(R["moon_alt"][0] > 0, R["moon_alt"][0], np.nan), color="0.35", ls=":", lw=2,
            label=f"Moon ({100 * ill:.0f}% lit)")
    cmap = plt.get_cmap("tab10")
    for k, t in enumerate(order):
        col = cmap(k % 10); s = stats[t]
        y = R["alt"][t, 0]; lim = R["lim"][t, 0]
        good = (y > lim) & R["dark"][0]
        ax.plot(offs, np.where(y > 0, y, np.nan), color=col, lw=1, alpha=0.45)
        ax.plot(offs, np.where(good, y, np.nan), color=col, lw=2.6,
                label=f"{cat.names[t]}  transit {clock(s['transit'])}, max {s['max_alt']:.0f}°, "
                      f"Moon sep {s['moon_sep']:.0f}°")
        if a.horizon and len(order) <= 6:
            ax.plot(offs, lim, color=col, lw=1, ls=(0, (1, 2)), alpha=0.8)
    if not a.horizon:
        ax.axhline(a.min_alt, color="k", lw=0.8, ls="--", alpha=0.5)
    # axes
    xt = np.arange(0, GRID_HOURS + 0.01, 1)
    ax.set_xticks(xt); ax.set_xticklabels([clock(x) for x in xt])
    ax.set_xlim(0, GRID_HOURS); ax.set_ylim(0, 90)
    ax.set_xlabel(f"local time (UTC{a.utc_offset:+g}), starting evening of {a.date}")
    ax.set_ylabel("altitude (deg)")
    top = ax.secondary_xaxis("top", functions=(lambda x: x, lambda x: x))
    top.set_xticks(xt); top.set_xticklabels([f"{int((GRID_START_LOCAL + x - a.utc_offset) % 24):02d}" for x in xt])
    top.set_xlabel("UTC (h)")
    am = np.array([1.0, 1.1, 1.2, 1.5, 2.0, 3.0, 5.0])
    right = ax.secondary_yaxis("right", functions=(lambda x: x, lambda x: x))
    right.set_yticks(np.rad2deg(np.arcsin(1 / am))); right.set_yticklabels([f"{v:g}" for v in am])
    right.set_ylabel("airmass")
    hz = "dotted = your horizon in each target's direction" if a.horizon else f"dashed = {a.min_alt:g}° floor"
    ax.set_title(f"Staralt-style plot, {a.date}  ({a.lat:.3f}, {a.lon:.3f})  bold = dark AND above horizon; {hz}",
                 fontsize=10)
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.07, 1.0), frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    os.makedirs(a.out, exist_ok=True)
    tag = "_".join(cat.names[t] for t in order[:3]) + ("_etc" if len(order) > 3 else "")
    p = os.path.join(a.out, f"staralt_{a.date}_{tag or 'none'}.png")
    fig.savefig(p, dpi=130); plt.close(fig)
    print(f"\nwrote {p}")


# =============================================================================== cli
def main():
    here = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="mode", required=True)
    ps = sub.add_parser("season", help="best month per target -> monthly CSVs + plots")
    ps.add_argument("start", help="first month, YYYY-MM")
    ps.add_argument("--months", type=int, default=12)
    ps.add_argument("--overlap", type=int, default=15, help="extend each month window this many days into the next")
    ps.add_argument("--out", default="season_out")
    ps.add_argument("--ref-time", default="21:00", help="clock time for the altitude-vs-day plot (HH:MM)")
    pa = sub.add_parser("staralt", help="altitude vs time for one night")
    pa.add_argument("date", help="YYYY-MM-DD (evening date)")
    pa.add_argument("targets", nargs="*", help="target names (substring match); none = top --top targets")
    pa.add_argument("--coord", action="append", default=[], help="'name,ra,dec' ad-hoc target (repeatable)")
    pa.add_argument("--resolve", action="store_true", help="look up unknown names online via Sesame")
    pa.add_argument("--top", type=int, default=8)
    pa.add_argument("--out", default="staralt_out")
    for sp in (ps, pa):
        sp.add_argument("--targets-file", default=os.path.join(here, "targets.csv"))
        sp.add_argument("--horizon", default=None, help="az,alt CSV of local obstructions")
        sp.add_argument("--min-alt", type=float, default=15.0)
        sp.add_argument("--twilight", type=float, default=-18.0)
        sp.add_argument("--session-end", default="03:00", help="HH:MM local or 'dawn'")
        sp.add_argument("--lat", type=float, default=DEF_LAT)
        sp.add_argument("--lon", type=float, default=DEF_LON)
        sp.add_argument("--elev", type=float, default=DEF_ELEV)
        sp.add_argument("--utc-offset", type=float, default=DEF_UTC)
        sp.add_argument("--step-min", type=float, default=5.0)
    a = p.parse_args()
    cat = read_targets(a.targets_file)
    limit_fn = read_horizon(a.horizon, a.min_alt)
    loc = EarthLocation(lat=a.lat * u.deg, lon=a.lon * u.deg, height=a.elev * u.m)
    (run_season if a.mode == "season" else run_staralt)(a, cat, loc, limit_fn)


if __name__ == "__main__":
    main()
