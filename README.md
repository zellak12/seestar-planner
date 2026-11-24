# Seestar Night Planner

An interactive observing planner for a smart telescope (built around the ZWO Seestar S50 Pro), running entirely in the browser. You don't need to install or run anything.

- **Month view**: each window is a calendar month plus the first 15 days of the next. It shows the targets whose best window this is, with three panels plotted against date: altitude at a chosen time, best time of night, and distance from the Moon. Moon brightness is shaded in the background. You can download the CSV.
- **Night view**: a Staralt-style altitude vs. local time plot for any date and any set of targets. It shows twilight bands, the Moon, an airmass axis, UTC, and your tree horizon, along with a table of transit, best time, airmass, usable window and Moon separation.
- **Year view**: one object's altitude through every night of the year as a height map (darkness only, with your horizon and session end marked), plus max altitude and usable hours per night for any targets you add.
- **Share links**: the URL stores the view, date, targets and any settings you've changed.

Default site: Tempe Town Lake, Tempe AZ (33.4322 N, 111.9319 W, UTC−7). You can change the site, darkness level, session end and tree horizon in the settings panel at the bottom of the page.

## How "best" is defined
For each night the planner counts the time a target spends above your horizon between full darkness (sun below −18°, or −12° if you pick nautical) and the session end (03:00 by default), weighting each moment by 1/airmass. A target is listed in the month window where that average is highest.

## Files
- `index.html`: the page
- `planner-core.js`: the calculations, ported from `python/seestar_plan.py`
- `targets.json`: the 63 targets. To add more, edit this file in GitHub's web editor (`name`, `ra` and `dec` in J2000 degrees, `type`).
- `python/`: the original command-line version, which produces monthly CSVs/PNGs and Staralt plots

## Credits
- Sun and Moon positions: [Astronomy Engine](https://github.com/cosinekitty/astronomy) (MIT), loaded from jsDelivr
- Charts: [Plotly.js](https://plotly.com/javascript/) (MIT), loaded from cdn.plot.ly
- Target positions: [OpenNGC](https://github.com/mattiaverga/OpenNGC) (CC-BY-SA 4.0); composite targets use the mean position of their members
- Target alt/az: J2000 coordinates precessed to the date, standard spherical astronomy, no refraction (good to about 0.1°)
