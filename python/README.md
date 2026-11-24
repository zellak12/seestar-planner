# Seestar planner (defaults: McClintock High School, Tempe AZ: 33.3981 N, -111.9069 E)

Requires: python3, numpy, astropy, matplotlib  (`pip install numpy astropy matplotlib`)

## Monthly workflow
1. Open `season_out/YYYY-MM.png` for the coming month. Each file covers the month **plus the first 15 days of the
   next month** (so files overlap), and shows only the targets whose **best window** is that one: altitude at 21:00 vs. day, best time of night, and Moon separation, with Moon brightness shaded.
   `season_out/YYYY-MM.csv` has the same per-night numbers (transit, best time, usable window, airmass, Moon).
2. Pick a target and a date, then make a Staralt-style plot for that night:

       python seestar_plan.py staralt 2027-01-12 Rosette
       python seestar_plan.py staralt 2027-01-12 Rosette M42 Cone --horizon horizon.csv
       python seestar_plan.py staralt 2027-01-12                 # top 8 targets that night
       python seestar_plan.py staralt 2027-01-12 "NGC 2174" --resolve        # online name lookup
       python seestar_plan.py staralt 2027-01-12 --coord "Thor,07:18:30,-13:13:48"

## Regenerate the season
    python seestar_plan.py season 2026-12 --months 12 [--overlap 15] [--horizon horizon.csv] [--session-end 03:00|dawn] [--ref-time 21:00]

"Best month" = the month where the target spends the most time above your horizon, weighted by 1/airmass,
during darkness up to --session-end (default 03:00). Use `--session-end dawn` if you image all night.

## Files
- `targets.csv`: 63 targets; J2000 positions from OpenNGC. Add your own lines (`name,ra,dec`).
- `horizon.csv`: **placeholder** tree profile (az, alt). Replace with your measurements and pass `--horizon`.
- `build_catalog.py`: regenerates targets.csv from OpenNGC (`pip install pyongc`).
- Defaults: 33.3981 N, -111.9069 E (McClintock HS), 360 m, UTC-7. Change with --lat/--lon/--elev/--utc-offset.
