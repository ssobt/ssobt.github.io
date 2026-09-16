/* ============================================================================
   cohort.js — interactive patient-cohort figure
   ----------------------------------------------------------------------------
   Reproduces the analysis in Figure 1C–D of Woo & Sobti et al. (2026): split a
   breast cancer cohort into bottom and top expression quartiles of one gene,
   then ask whether the top quartile is more spread out in transcriptome space.

   Two views over the same patients, with positions lerped between them:

     PC space  — PC1/PC2 scatter, per-group centroid, and a ring at each
                 group's median 2-D distance from its own centroid.
     Distance  — the actual metric: distance to group centroid across the top
                 50 PCs, as a jittered strip with median and IQR.

   The ring is drawn from 2-D distances and the readout reports the 50-PC
   metric; each is labelled for what it is, because a ring drawn at a 50-D
   radius on a 2-D plot would be a lie.

   Only ~1,200 points, so hover is a plain linear scan — no spatial index.
   ========================================================================= */

const TAU = Math.PI * 2;
const MAX_DPR = 2;

const COL = {
  q25:  '#56B4E9',
  q75:  '#E69F00',
  mid:  '#3A4352',
  hot:  '#FDE725',
  axis: '#2C3646',
  text: '#8B98A9',
};

const Q25 = 0, MID = 1, Q75 = 2;

const $ = (sel, root = document) => root.querySelector(sel);
const el = (tag, cls, text) => {
  const n = document.createElement(tag);
  if (cls) n.className = cls;
  if (text != null) n.textContent = text;
  return n;
};
const fmt = (n) => n.toLocaleString('en-US');
const easeInOutCubic = (t) =>
  t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
const reducedMotion = () =>
  matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches;

/* --- stats helpers ------------------------------------------------------- */

function quantile(sorted, q) {
  if (!sorted.length) return 0;
  const pos = (sorted.length - 1) * q;
  const lo = Math.floor(pos), hi = Math.ceil(pos);
  return lo === hi ? sorted[lo] : sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo);
}

/** Round tick values spanning [lo, hi], roughly `count` of them. */
function niceTicks(lo, hi, count) {
  const raw = (hi - lo) / Math.max(1, count);
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag)
    .find((v) => v >= raw) || 10 * mag;
  const out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi; v += step) {
    out.push(Math.round(v * 100) / 100);
  }
  return out;
}

/** Median, Q1 and Q3 of the values at `idx`. */
function summary(values, idx) {
  const s = idx.map((i) => values[i]).sort((a, b) => a - b);
  return { n: s.length, med: quantile(s, 0.5), q1: quantile(s, 0.25), q3: quantile(s, 0.75) };
}

/* --- main ---------------------------------------------------------------- */

export function mountCohort(data, root) {
  const canvas = $('#cohort-canvas', root);
  const tip = $('#cohort-tip', root);
  const geneSeg = $('#cohort-gene', root);
  const viewSeg = $('#cohort-view', root);
  const readout = $('#cohort-readout', root);
  const capEl = $('#cohort-caption', root);
  const tableWrap = $('#cohort-table');

  const ctx = canvas.getContext('2d');
  const n = data.meta.n;
  const pc = data.coords.pc;

  let geneKey = data.geneOrder[0];
  let view = 'pc';                       // 'pc' | 'distance'
  let blend = 0;                         // 0 = PC layout, 1 = distance layout
  let anim = null;
  let hoverIdx = -1;
  let cssW = 0, cssH = 0, dpr = 1, raf = 0, running = true;

  /* -- layout computation -------------------------------------------------- */

  // Normalised PC positions, aspect preserved. Scaled on a robust percentile
  // range rather than min/max: a couple of extreme tumours otherwise squeeze
  // the whole cohort into the middle few percent of the plot, which both
  // hides the structure and shrinks the dispersion rings to nothing.
  const pcPos = (() => {
    const col = (d) => {
      const v = new Array(n);
      for (let i = 0; i < n; i++) v[i] = pc[i * 2 + d];
      return v.sort((a, b) => a - b);
    };
    const xs = col(0), ys = col(1);
    const span = Math.max(
      quantile(xs, 0.99) - quantile(xs, 0.01),
      quantile(ys, 0.99) - quantile(ys, 0.01),
    ) * 1.25 || 1;
    const c = [quantile(xs, 0.5), quantile(ys, 0.5)];
    const out = new Float32Array(n * 2);
    for (let i = 0; i < n; i++) {
      out[i * 2]     = (pc[i * 2] - c[0]) / span + 0.5;
      out[i * 2 + 1] = 0.5 - (pc[i * 2 + 1] - c[1]) / span;
    }
    return out;
  })();

  // Deterministic jitter for the strip view — stable across redraws. A
  // golden-ratio sequence bands into visible stripes here, so hash instead.
  const jitter = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    let h = Math.imul(i ^ 0x9e3779b9, 0x85ebca6b);
    h ^= h >>> 13; h = Math.imul(h, 0xc2b2ae35); h ^= h >>> 16;
    jitter[i] = ((h >>> 0) / 4294967296) - 0.5;
  }

  /** Per-gene derived quantities, computed once and cached. */
  const cache = new Map();
  function gene(key) {
    if (cache.has(key)) return cache.get(key);
    const g = data.genes[key];
    const idx = { [Q25]: [], [MID]: [], [Q75]: [] };
    for (let i = 0; i < n; i++) idx[g.quartile[i]].push(i);

    // 2-D centroid and median 2-D radius per group — what the ring shows.
    const ring = {};
    for (const q of [Q25, Q75]) {
      const ids = idx[q];
      let cx = 0, cy = 0;
      for (const i of ids) { cx += pcPos[i * 2]; cy += pcPos[i * 2 + 1]; }
      cx /= ids.length; cy /= ids.length;
      const d2 = ids
        .map((i) => Math.hypot(pcPos[i * 2] - cx, pcPos[i * 2 + 1] - cy))
        .sort((a, b) => a - b);
      ring[q] = { cx, cy, med: quantile(d2, 0.5) };
    }

    // The published metric: distance across the top 50 PCs.
    const stats = {
      [Q25]: summary(g.distance, idx[Q25]),
      [Q75]: summary(g.distance, idx[Q75]),
    };

    // Strip-view layout, scaled to the pooled distance range.
    const all = [...idx[Q25], ...idx[Q75]].map((i) => g.distance[i]).sort((a, b) => a - b);
    const dLo = Math.min(0, quantile(all, 0)), dHi = quantile(all, 1) * 1.04;

    const stripPos = new Float32Array(n * 2);
    for (let i = 0; i < n; i++) {
      const q = g.quartile[i];
      const col = q === Q25 ? 0.30 : q === Q75 ? 0.70 : 0.5;
      stripPos[i * 2] = col + jitter[i] * 0.20;
      stripPos[i * 2 + 1] = 1 - (g.distance[i] - dLo) / (dHi - dLo || 1);
    }

    const out = { g, idx, ring, stats, stripPos, dLo, dHi };
    cache.set(key, out);
    return out;
  }

  /* -- geometry ------------------------------------------------------------ */

  const PAD = { l: 56, r: 18, t: 18, b: 40 };

  function measure() {
    const r = canvas.getBoundingClientRect();
    if (!r.width || !r.height) return false;
    cssW = r.width; cssH = r.height;
    dpr = Math.min(devicePixelRatio || 1, MAX_DPR);
    const bw = Math.round(cssW * dpr), bh = Math.round(cssH * dpr);
    if (canvas.width !== bw || canvas.height !== bh) { canvas.width = bw; canvas.height = bh; }
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return true;
  }

  const plotW = () => cssW - PAD.l - PAD.r;
  const plotH = () => cssH - PAD.t - PAD.b;

  /** Screen position of point i at the current blend. */
  function pos(i, G) {
    const t = blend;
    const w = plotW(), h = plotH();
    // PC layout keeps a square aspect; strip layout fills the box.
    const side = Math.min(w, h);
    const ox = PAD.l + (w - side) / 2, oy = PAD.t + (h - side) / 2;
    const px = ox + pcPos[i * 2] * side;
    const py = oy + pcPos[i * 2 + 1] * side;
    const sx = PAD.l + G.stripPos[i * 2] * w;
    const sy = PAD.t + G.stripPos[i * 2 + 1] * h;
    return [px + (sx - px) * t, py + (sy - py) * t];
  }

  /* -- drawing ------------------------------------------------------------- */

  function drawAxes(G) {
    ctx.save();
    ctx.font = '500 10px ui-monospace, SF Mono, Menlo, monospace';
    ctx.fillStyle = COL.text;
    ctx.strokeStyle = COL.axis;
    ctx.lineWidth = 1;

    if (blend < 0.5) {
      ctx.globalAlpha = 1 - blend * 2;
      ctx.textAlign = 'center';
      ctx.fillText('PC1', PAD.l + plotW() / 2, cssH - 12);
      ctx.save();
      ctx.translate(14, PAD.t + plotH() / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.fillText('PC2', 0, 0);
      ctx.restore();
    } else {
      ctx.globalAlpha = (blend - 0.5) * 2;
      // y axis: distance from centroid
      ctx.beginPath();
      ctx.moveTo(PAD.l - 6, PAD.t);
      ctx.lineTo(PAD.l - 6, PAD.t + plotH());
      ctx.stroke();
      ctx.textAlign = 'right';
      for (const val of niceTicks(G.dLo, G.dHi, 5)) {
        const frac = 1 - (val - G.dLo) / (G.dHi - G.dLo || 1);
        if (frac < -0.01 || frac > 1.01) continue;
        const y = PAD.t + plotH() * frac;
        ctx.fillText(String(val), PAD.l - 11, y + 3.5);
        ctx.beginPath();
        ctx.moveTo(PAD.l - 9, y); ctx.lineTo(PAD.l - 6, y); ctx.stroke();
      }
      ctx.save();
      ctx.translate(15, PAD.t + plotH() / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.textAlign = 'center';
      ctx.fillText('distance from centroid', 0, 0);
      ctx.restore();

      ctx.textAlign = 'center';
      ctx.fillText('bottom quartile', PAD.l + plotW() * 0.30, cssH - 12);
      ctx.fillText('top quartile', PAD.l + plotW() * 0.70, cssH - 12);
    }
    ctx.restore();
  }

  function drawRings(G) {
    if (blend > 0.02) return;
    const w = plotW(), h = plotH(), side = Math.min(w, h);
    const ox = PAD.l + (w - side) / 2, oy = PAD.t + (h - side) / 2;
    ctx.save();
    ctx.globalAlpha = (1 - blend) * 0.9;
    ctx.setLineDash([4, 4]);
    for (const q of [Q25, Q75]) {
      const r = G.ring[q];
      const cx = ox + r.cx * side, cy = oy + r.cy * side;
      ctx.strokeStyle = q === Q25 ? COL.q25 : COL.q75;
      ctx.lineWidth = 1.4;
      ctx.beginPath();
      ctx.arc(cx, cy, r.med * side, 0, TAU);
      ctx.stroke();
      // centroid marker
      ctx.setLineDash([]);
      ctx.beginPath();
      ctx.moveTo(cx - 5, cy); ctx.lineTo(cx + 5, cy);
      ctx.moveTo(cx, cy - 5); ctx.lineTo(cx, cy + 5);
      ctx.lineWidth = 1.6;
      ctx.stroke();
      ctx.setLineDash([4, 4]);
    }
    ctx.restore();
  }

  function drawBoxes(G) {
    if (blend < 0.98) return;
    const w = plotW(), h = plotH();
    const yOf = (v) => PAD.t + h * (1 - (v - G.dLo) / (G.dHi - G.dLo || 1));
    ctx.save();
    ctx.globalAlpha = blend;
    for (const q of [Q25, Q75]) {
      const s = G.stats[q];
      const cx = PAD.l + w * (q === Q25 ? 0.30 : 0.70);
      const bw = Math.min(104, w * 0.19);
      const top = yOf(s.q3), bot = yOf(s.q1);
      ctx.fillStyle = 'rgba(10,12,16,0.55)';
      ctx.fillRect(cx - bw / 2, top, bw, bot - top);
      ctx.strokeStyle = q === Q25 ? COL.q25 : COL.q75;
      ctx.lineWidth = 1.2;
      ctx.strokeRect(cx - bw / 2, top, bw, bot - top);
      // Median in the foreground colour so it reads against the points.
      ctx.strokeStyle = '#E6EDF3';
      ctx.lineWidth = 2.6;
      ctx.beginPath();
      ctx.moveTo(cx - bw / 2 - 4, yOf(s.med));
      ctx.lineTo(cx + bw / 2 + 4, yOf(s.med));
      ctx.stroke();
    }
    ctx.restore();
  }

  function draw() {
    if (!measure()) return;
    const G = gene(geneKey);
    ctx.clearRect(0, 0, cssW, cssH);

    drawAxes(G);
    ctx.save();
    ctx.beginPath();
    ctx.rect(PAD.l - 2, PAD.t - 2, plotW() + 4, plotH() + 4);
    ctx.clip();
    drawRings(G);

    const r = 2.3;
    // Middle-quartile patients are context, not part of the comparison.
    // Context in PC space; not part of the comparison, so gone in the strip.
    ctx.globalAlpha = 0.46 * (1 - blend);
    ctx.fillStyle = COL.mid;
    ctx.beginPath();
    for (const i of G.idx[MID]) {
      const [x, y] = pos(i, G);
      ctx.moveTo(x + r, y); ctx.arc(x, y, r, 0, TAU);
    }
    ctx.fill();

    ctx.globalAlpha = 0.9;
    for (const q of [Q25, Q75]) {
      ctx.fillStyle = q === Q25 ? COL.q25 : COL.q75;
      ctx.beginPath();
      for (const i of G.idx[q]) {
        const [x, y] = pos(i, G);
        ctx.moveTo(x + r + 0.4, y); ctx.arc(x, y, r + 0.4, 0, TAU);
      }
      ctx.fill();
    }

    drawBoxes(G);
    ctx.restore();

    if (hoverIdx >= 0) {
      const [x, y] = pos(hoverIdx, G);
      ctx.globalAlpha = 1;
      ctx.strokeStyle = COL.hot;
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.arc(x, y, r + 4, 0, TAU);
      ctx.stroke();
    }
    ctx.globalAlpha = 1;
  }

  function schedule() {
    if (!raf && running) raf = requestAnimationFrame(() => { raf = 0; draw(); });
  }

  /* -- view transition ----------------------------------------------------- */

  function setView(next) {
    if (next === view) return;
    view = next;
    const from = blend, to = next === 'pc' ? 0 : 1;
    tip.dataset.on = 'false';
    hoverIdx = -1;

    if (reducedMotion()) { blend = to; draw(); return; }
    const start = performance.now(), dur = 780;
    if (anim) cancelAnimationFrame(anim);
    const step = (now) => {
      const p = Math.min(1, (now - start) / dur);
      blend = from + (to - from) * easeInOutCubic(p);
      draw();
      if (p < 1) anim = requestAnimationFrame(step); else anim = null;
    };
    anim = requestAnimationFrame(step);
  }

  /* -- hover --------------------------------------------------------------- */

  canvas.addEventListener('mousemove', (ev) => {
    const G = gene(geneKey);
    const rect = canvas.getBoundingClientRect();
    const mx = ev.clientX - rect.left, my = ev.clientY - rect.top;
    let best = -1, bestD = 9 * 9;
    for (let i = 0; i < n; i++) {
      const [x, y] = pos(i, G);
      const d = (x - mx) * (x - mx) + (y - my) * (y - my);
      if (d < bestD) { bestD = d; best = i; }
    }
    if (best === hoverIdx) return;
    hoverIdx = best;
    schedule();
    if (best < 0) { tip.dataset.on = 'false'; return; }

    const q = G.g.quartile[best];
    tip.replaceChildren();
    const row = (k, v) => {
      const d = el('div');
      d.append(el('dt', null, `${k}:`), el('dd', null, v));
      tip.append(d);
    };
    row('patient', `#${best + 1}`);
    row(`${geneKey}`, `${G.g.expression[best]} FPKM`);
    row('group', q === Q25 ? 'bottom quartile' : q === Q75 ? 'top quartile' : 'middle 50% (not compared)');
    row('distance', G.g.distance[best].toFixed(1));
    const [x, y] = pos(best, G);
    tip.style.left = `${x}px`;
    tip.style.top = `${y}px`;
    tip.dataset.on = 'true';
  });

  canvas.addEventListener('mouseleave', () => {
    hoverIdx = -1; tip.dataset.on = 'false'; schedule();
  });

  /* -- readout ------------------------------------------------------------- */

  function renderReadout() {
    const G = gene(geneKey);
    const s25 = G.stats[Q25], s75 = G.stats[Q75];
    const ratio = s75.med / s25.med;
    const pub = G.g.published;
    const isControl = G.g.role === 'control';

    readout.replaceChildren();

    const head = el('p', 'readout__head');
    head.append(el('strong', null, geneKey));
    head.append(el('span', 'chip' + (isControl ? '' : ' chip--accent'),
                   isControl ? 'negative control' : 'candidate regulator'));
    readout.append(head);

    const dl = el('dl', 'readout__grid');
    const pair = (k, v, cls) => {
      dl.append(el('dt', null, k));
      dl.append(el('dd', cls, v));
    };
    pair('bottom quartile', `n = ${fmt(s25.n)} · median ${s25.med.toFixed(1)}`);
    pair('top quartile', `n = ${fmt(s75.n)} · median ${s75.med.toFixed(1)}`);
    pair('spread ratio', `${ratio.toFixed(2)}×`, ratio > 1.25 ? 'is-hit' : 'is-null');
    readout.append(dl);

    const verdict = el('p', 'readout__verdict');
    verdict.textContent = isControl
      ? 'The two quartiles are indistinguishable — which is the point of a control.'
      : 'The top quartile is measurably more dispersed at the same mean expression.';
    readout.append(verdict);

    const cite = el('div', 'readout__cite');
    cite.append(el('p', 'readout__citehead', 'Published result for this gene'));
    const cdl = el('dl', 'readout__grid');
    const cpair = (k, v) => { cdl.append(el('dt', null, k)); cdl.append(el('dd', null, v)); };
    cpair('centroid distance', `P = ${pub.centroidP}`);
    cpair('METABRIC survival', `HR = ${pub.hr} · P = ${pub.survP}`);
    cite.append(cdl);
    const src = el('p', 'readout__src');
    const a = el('a', null, 'Woo & Sobti et al., bioRxiv 2026');
    a.href = `https://doi.org/${data.meta.published.doi}`;
    a.rel = 'noopener';
    src.append(document.createTextNode('Fig. 1C–D, 6C · '), a);
    cite.append(src);
    readout.append(cite);

    renderTable(G);
  }

  function renderTable(G) {
    if (!tableWrap) return;
    tableWrap.replaceChildren();
    const t = el('table', 'tbl');
    const cap = el('caption', 'visually-hidden',
      `Spread statistics for ${geneKey} quartile groups.`);
    const thead = el('thead'), hr = el('tr');
    ['Group', 'Patients', 'Median distance', 'IQR'].forEach((h, i) => {
      const th = el('th', i ? 'num' : null, h); th.scope = 'col'; hr.append(th);
    });
    thead.append(hr);
    const tb = el('tbody');
    for (const [q, name] of [[Q25, 'Bottom quartile'], [Q75, 'Top quartile']]) {
      const s = G.stats[q], tr = el('tr');
      tr.append(el('td', null, name));
      tr.append(el('td', 'num', fmt(s.n)));
      tr.append(el('td', 'num', s.med.toFixed(1)));
      tr.append(el('td', 'num', `${s.q1.toFixed(1)}–${s.q3.toFixed(1)}`));
      tb.append(tr);
    }
    t.append(cap, thead, tb);
    const w = el('div', 'tbl-wrap'); w.append(t);
    tableWrap.append(w);
  }

  /* -- controls ------------------------------------------------------------ */

  function buildSeg(container, items, current, onPick) {
    container.replaceChildren();
    for (const it of items) {
      const b = el('button', null, it.label);
      b.type = 'button';
      b.dataset.key = it.key;
      b.setAttribute('aria-pressed', String(it.key === current));
      b.addEventListener('click', () => onPick(it));
      container.append(b);
    }
  }

  const geneItems = data.geneOrder.map((k) => ({
    key: k,
    label: data.genes[k].role === 'control' ? `${k} (ctrl)` : k,
  }));

  function selectGene(it) {
    geneKey = it.key;
    hoverIdx = -1; tip.dataset.on = 'false';
    buildSeg(geneSeg, geneItems, geneKey, selectGene);
    renderReadout();
    draw();
  }

  const viewItems = [
    { key: 'pc', label: 'PC space' },
    { key: 'distance', label: 'Distance' },
  ];
  function selectView(it) {
    buildSeg(viewSeg, viewItems, it.key, selectView);
    setView(it.key);
  }

  buildSeg(geneSeg, geneItems, geneKey, selectGene);
  buildSeg(viewSeg, viewItems, view, selectView);
  renderReadout();

  if (capEl && data.meta.caption) {
    capEl.replaceChildren(el('b', null, 'Figure 1. '),
                          document.createTextNode(data.meta.caption));
  }

  new ResizeObserver(() => { if (measure()) draw(); }).observe(canvas);
  new IntersectionObserver((e) => {
    running = e[0].isIntersecting;
    if (running) schedule();
  }, { threshold: 0 }).observe(canvas);

  measure();
  draw();
}

/** Fetch + validate cohort.json. */
export async function loadCohort(url) {
  const res = await fetch(url, { cache: 'no-cache' });
  if (!res.ok) throw new Error(`${res.status} ${res.statusText} fetching ${url}`);
  const d = await res.json();
  if (!d.coords || !Array.isArray(d.coords.pc)) throw new Error('missing coords.pc');
  const n = d.coords.pc.length >> 1;
  if (!d.meta) d.meta = {};
  d.meta.n = n;
  if (!d.genes || !Object.keys(d.genes).length) throw new Error('no genes in payload');
  if (!Array.isArray(d.geneOrder) || !d.geneOrder.length) d.geneOrder = Object.keys(d.genes);
  for (const k of d.geneOrder) {
    const g = d.genes[k];
    if (!g) throw new Error(`geneOrder lists "${k}" but genes.${k} is missing`);
    for (const f of ['expression', 'quartile', 'distance']) {
      if (!Array.isArray(g[f]) || g[f].length !== n) {
        throw new Error(`genes.${k}.${f} must be an array of length ${n}`);
      }
    }
  }
  return d;
}
