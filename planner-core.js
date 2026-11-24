/* planner-core.js: observability math for the Seestar planner.
   Port of seestar_plan.py. Sun/Moon from Astronomy Engine (MIT, github.com/cosinekitty/astronomy).
   Works in the browser (global `Astronomy`) and in Node (require). */
(function (root) {
  "use strict";
  const A = root.Astronomy || (typeof require !== "undefined" ? require("astronomy-engine") : null);
  const D2R = Math.PI / 180, R2D = 180 / Math.PI;
  const GRID_START = 17, GRID_HOURS = 14, SIDEREAL = 1.00273790935;

  const DEFAULTS = {
    siteName: "McClintock High School, Tempe AZ",
    lat: 33.3981, lon: -111.9069, elev: 360, utc: -7,
    minAlt: 15, twilight: -18, sessionEnd: "03:00", overlap: 15, refTime: "21:00",
    horizon: [] // [[az, alt], ...]
  };

  // ------------------------------------------------------------------ helpers
  function parseClock(s) { // "HH:MM" or "dawn" -> hours after 17:00
    if (String(s).toLowerCase() === "dawn") return GRID_HOURS;
    const [h, m] = String(s).split(":").map(Number);
    return (((h + (m || 0) / 60 - GRID_START) % 24) + 24) % 24;
  }
  function clock(h) {
    if (h == null || !isFinite(h)) return "";
    const t = (((GRID_START + h) % 24) + 24) % 24;
    const m = Math.round(t * 60) % 1440;
    return String(Math.floor(m / 60)).padStart(2, "0") + ":" + String(m % 60).padStart(2, "0");
  }
  const COMPASS = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"];
  const compass = az => COMPASS[Math.floor(((az % 360) + 360) % 360 / 22.5 + 0.5) % 16];

  function parseRA(s) { // "hh:mm:ss" -> deg, or decimal deg
    s = String(s).trim();
    if (/[:\s]/.test(s)) { const p = s.split(/[:\s]+/).map(Number); return 15 * (p[0] + (p[1] || 0) / 60 + (p[2] || 0) / 3600); }
    return Number(s);
  }
  function parseDec(s) {
    s = String(s).trim();
    if (/[:\s]/.test(s)) {
      const neg = s.startsWith("-"); const p = s.replace(/^[+-]/, "").split(/[:\s]+/).map(Number);
      const v = p[0] + (p[1] || 0) / 60 + (p[2] || 0) / 3600; return neg ? -v : v;
    }
    return Number(s);
  }

  function horizonFn(settings) {
    const pts = (settings.horizon || []).map(([a, h]) => [((a % 360) + 360) % 360, h]).sort((x, y) => x[0] - y[0]);
    const floor = settings.minAlt;
    if (!pts.length) return () => floor;
    const az = [], alt = [];
    for (const k of [-360, 0, 360]) for (const [a, h] of pts) { az.push(a + k); alt.push(h); }
    return function (a) {
      a = ((a % 360) + 360) % 360;
      let i = 1; while (i < az.length && az[i] < a) i++;
      const f = (a - az[i - 1]) / (az[i] - az[i - 1] || 1);
      return Math.max(floor, alt[i - 1] + f * (alt[i] - alt[i - 1]));
    };
  }

  function dateAdd(dateStr, days) {
    const [y, m, d] = dateStr.split("-").map(Number);
    const t = new Date(Date.UTC(y, m - 1, d + days));
    return t.toISOString().slice(0, 10);
  }

  // ------------------------------------------------------------------ core
  /* Compute everything for a list of nights (evening dates) and targets on a time grid. */
  function computeNights(targets, dates, settings, stepMin) {
    const obs = new A.Observer(settings.lat, settings.lon, settings.elev || 0);
    const S = Math.round(GRID_HOURS * 60 / stepMin) + 1;
    const offs = Float64Array.from({ length: S }, (_, i) => i * stepMin / 60);
    const N = dates.length, T = targets.length;
    const phi = settings.lat * D2R, sphi = Math.sin(phi), cphi = Math.cos(phi);
    const lim = horizonFn(settings);
    const sessEnd = parseClock(settings.sessionEnd);

    // precess target coordinates to the middle of the period (Astronomy Engine rotation)
    const [y0, m0, d0] = dates[Math.floor(N / 2)].split("-").map(Number);
    const tMid = A.MakeTime(new Date(Date.UTC(y0, m0 - 1, d0)));
    const rot = A.Rotation_EQJ_EQD(tMid);
    const tRA = new Float64Array(T), tDec = new Float64Array(T);
    targets.forEach((tg, i) => {
      const r = tg.ra * D2R, d = tg.dec * D2R;
      const v = new A.Vector(Math.cos(d) * Math.cos(r), Math.cos(d) * Math.sin(r), Math.sin(d), tMid);
      const e = A.EquatorFromVector(A.RotateVector(rot, v));
      tRA[i] = e.ra * 15 * D2R; tDec[i] = e.dec * D2R;
    });

    const night = [];
    for (let n = 0; n < N; n++) {
      const [y, m, d] = dates[n].split("-").map(Number);
      const base = Date.UTC(y, m - 1, d, GRID_START) - settings.utc * 3600e3;
      const lst = new Float64Array(S), sunAlt = new Float64Array(S), moonAlt = new Float64Array(S), moonAz = new Float64Array(S);
      const dark = new Uint8Array(S), session = new Uint8Array(S);
      for (let s = 0; s < S; s++) {
        const t = A.MakeTime(new Date(base + offs[s] * 3600e3));
        lst[s] = ((A.SiderealTime(t) + settings.lon / 15) * 15) * D2R;
        const sun = A.Equator("Sun", t, obs, true, true);
        sunAlt[s] = altaz(sun.ra * 15 * D2R, sun.dec * D2R, lst[s], sphi, cphi)[0];
        const moon = A.Equator("Moon", t, obs, true, true);
        const aa = altaz(moon.ra * 15 * D2R, moon.dec * D2R, lst[s], sphi, cphi);
        moonAlt[s] = aa[0]; moonAz[s] = aa[1];
        dark[s] = sunAlt[s] < settings.twilight ? 1 : 0;
        session[s] = dark[s] && offs[s] <= sessEnd + 1e-9 ? 1 : 0;
      }
      const tMidnight = A.MakeTime(new Date(base + 7 * 3600e3));
      const illum = A.Illumination("Moon", tMidnight).phase_fraction;
      const alt = [], az = [], lm = [];
      for (let k = 0; k < T; k++) {
        const a1 = new Float32Array(S), z1 = new Float32Array(S), l1 = new Float32Array(S);
        for (let s = 0; s < S; s++) {
          const r = altaz(tRA[k], tDec[k], lst[s], sphi, cphi);
          a1[s] = r[0]; z1[s] = r[1]; l1[s] = lim(r[1]);
        }
        alt.push(a1); az.push(z1); lm.push(l1);
      }
      // transit (hours after 17:00) from LST at grid start
      const transit = Array.from(tRA, ra => ((((ra - lst[0]) % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI)) * 12 / Math.PI / SIDEREAL);
      night.push({ date: dates[n], base, lst, sunAlt, moonAlt, moonAz, dark, session, illum, alt, az, lim: lm, transit });
    }
    return { offs, stepMin, targets, nights: night };
  }

  function altaz(ra, dec, lst, sphi, cphi) {
    const H = lst - ra, sd = Math.sin(dec), cd = Math.cos(dec), cH = Math.cos(H);
    const alt = Math.asin(Math.max(-1, Math.min(1, sphi * sd + cphi * cd * cH))) * R2D;
    let az = Math.atan2(-cd * Math.sin(H), sd * cphi - cd * cH * sphi) * R2D;
    return [alt, (az + 360) % 360];
  }

  function angsep(a1, z1, a2, z2) {
    a1 *= D2R; z1 *= D2R; a2 *= D2R; z2 *= D2R;
    return Math.acos(Math.max(-1, Math.min(1, Math.sin(a1) * Math.sin(a2) + Math.cos(a1) * Math.cos(a2) * Math.cos(z1 - z2)))) * R2D;
  }

  /* Summary of target k on night n (same definitions as the Python script). */
  function nightStats(C, k, n) {
    const N = C.nights[n], offs = C.offs, alt = N.alt[k], az = N.az[k], lim = N.lim[k];
    const step = C.stepMin / 60;
    const ses = [], all = []; let score = 0, darkIdx = [];
    for (let s = 0; s < offs.length; s++) {
      if (N.dark[s]) darkIdx.push(s);
      const up = alt[s] > lim[s];
      if (up && N.session[s]) { ses.push(s); score += Math.sin(Math.max(0, alt[s]) * D2R) * step; }
      if (up && N.dark[s]) all.push(s);
    }
    const pool = ses.length ? ses : darkIdx;
    let best = pool.length ? pool[0] : 0;
    for (const s of (pool.length ? pool : offs.map((_, i) => i))) if (alt[s] > alt[best]) best = s;
    const ma = alt[best];
    return {
      transit: N.transit[k],
      usableSession: ses.length * step, usableNight: all.length * step,
      start: ses.length ? offs[ses[0]] : NaN, end: ses.length ? offs[ses[ses.length - 1]] : NaN,
      nightStart: all.length ? offs[all[0]] : NaN, nightEnd: all.length ? offs[all[all.length - 1]] : NaN,
      best: offs[best], maxAlt: ma, az: az[best], lim: lim[best],
      airmass: ma > 1 ? 1 / Math.sin(ma * D2R) : Infinity,
      score, illum: N.illum, moonAlt: N.moonAlt[best],
      moonSep: angsep(ma, az[best], N.moonAlt[best], N.moonAz[best]),
      darkStart: darkIdx.length ? offs[darkIdx[0]] : NaN, darkEnd: darkIdx.length ? offs[darkIdx[darkIdx.length - 1]] : NaN,
    };
  }

  /* Windows: calendar month + `overlap` days of the next month, for `months` months from startYM. */
  function windows(startYM, months, overlap) {
    let [y, m] = startYM.split("-").map(Number);
    const out = [];
    for (let i = 0; i < months; i++) {
      const yy = y + Math.floor((m - 1 + i) / 12), mm = ((m - 1 + i) % 12) + 1;
      const first = `${yy}-${String(mm).padStart(2, "0")}-01`;
      const daysInMonth = new Date(Date.UTC(yy, mm, 0)).getUTCDate();
      const last = dateAdd(first, daysInMonth - 1 + overlap);
      out.push({ key: `${yy}-${String(mm).padStart(2, "0")}`, year: yy, month: mm, first, last });
    }
    return out;
  }
  function dateRange(first, last) {
    const out = []; let d = first;
    while (d <= last) { out.push(d); d = dateAdd(d, 1); }
    return out;
  }

  /* Assign every target to its best window (coarse 10-min grid). */
  function assignSeason(targets, settings, startYM, months, onProgress) {
    const W = windows(startYM, months || 12, settings.overlap);
    const dates = dateRange(W[0].first, W[W.length - 1].last);
    const C = computeNights(targets, dates, settings, 10);
    const idx = Object.fromEntries(dates.map((d, i) => [d, i]));
    const score = targets.map(() => new Float64Array(dates.length));
    for (let n = 0; n < dates.length; n++) for (let k = 0; k < targets.length; k++) score[k][n] = nightStats(C, k, n).score;
    const best = {}, never = [];
    targets.forEach((tg, k) => {
      let bw = null, bs = 0;
      for (const w of W) {
        let s = 0, c = 0;
        for (let n = idx[w.first]; n <= idx[w.last]; n++) { s += score[k][n]; c++; }
        if (s / c > bs) { bs = s / c; bw = w.key; }
      }
      if (bw) best[tg.name] = bw; else never.push(tg.name);
    });
    return { windows: W, best, never };
  }

  const api = { DEFAULTS, GRID_START, GRID_HOURS, parseClock, clock, compass, parseRA, parseDec, horizonFn,
    dateAdd, dateRange, windows, computeNights, nightStats, assignSeason, angsep };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.Planner = api;
})(typeof window !== "undefined" ? window : globalThis);
