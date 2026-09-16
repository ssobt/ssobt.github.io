/* ============================================================================
   scatter.js — canvas renderer for single-cell embeddings
   ----------------------------------------------------------------------------
   No dependencies. Handles 5k–50k points at 60fps by doing three things:

     1. Canvas 2D, not SVG. SVG stalls past ~2k nodes once hover is involved.
     2. Draw calls batched by color — one path per color, not one per point.
     3. Hover resolved through a uniform spatial grid (CSR layout), so a
        mousemove costs O(points in a few cells) instead of O(n).

   Coordinates are normalized per embedding into a unit box with a *single*
   scale factor for both axes, so the embedding's aspect ratio is preserved —
   a stretched UMAP is a misleading UMAP.

   The hero additionally runs a 3D intro: a rotating PC1/PC2/PC3 cloud that
   collapses into the 2D embedding. The 3D is deliberately transient. A
   persistent 3D scatter occludes its own points, which makes cluster density
   and relative group size unreadable — fine as a two-second gesture depicting
   dimensionality reduction, not acceptable as a figure you'd ask someone to
   read. Everything at rest, and the whole of §02, is 2D.
   ========================================================================= */

const TAU = Math.PI * 2;
const GRID = 64;          // spatial grid resolution (GRID × GRID cells)
const CONT_BINS = 64;     // quantization bins for continuous color batching
const MAX_DPR = 2;        // cap: 3× backing stores cost a lot for no gain
const DEPTH_SLICES = 14;  // back-to-front slices during the 3D intro

/* --- viridis ------------------------------------------------------------- */

const VIRIDIS = [
  [68, 1, 84],    [72, 40, 120],  [62, 74, 137],  [49, 104, 142],
  [38, 130, 142], [31, 158, 137], [53, 183, 121], [109, 205, 89],
  [180, 222, 44], [253, 231, 37],
];

/** Sample the viridis ramp at t ∈ [0,1]. Returns "rgb(r,g,b)". */
export function viridis(t) {
  t = t < 0 ? 0 : t > 1 ? 1 : t;
  const x = t * (VIRIDIS.length - 1);
  const i = Math.min(Math.floor(x), VIRIDIS.length - 2);
  const f = x - i;
  const a = VIRIDIS[i], b = VIRIDIS[i + 1];
  return `rgb(${Math.round(a[0] + (b[0] - a[0]) * f)},` +
         `${Math.round(a[1] + (b[1] - a[1]) * f)},` +
         `${Math.round(a[2] + (b[2] - a[2]) * f)})`;
}

/** CSS gradient string for a colorbar that matches the point colors. */
export function viridisGradient() {
  const stops = VIRIDIS.map((c, i) => {
    const pct = Math.round((i / (VIRIDIS.length - 1)) * 100);
    return `rgb(${c[0]},${c[1]},${c[2]}) ${pct}%`;
  });
  return `linear-gradient(90deg, ${stops.join(', ')})`;
}

/* Okabe–Ito, reordered so the first few are maximally distinct on dark.
   Colorblind-safe; cycles if a label set has more levels than colors. */
export const CATEGORICAL = [
  '#56B4E9', '#E69F00', '#009E73', '#CC79A7',
  '#F0E442', '#0072B2', '#D55E00', '#8FA7B8',
];

export const categoricalColor = (i) => CATEGORICAL[i % CATEGORICAL.length];

/* --- helpers ------------------------------------------------------------- */

const easeInOutCubic = (t) =>
  t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;

const prefersReducedMotion = () =>
  window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;

/**
 * Normalize a flat [x0,y0,…] array into a unit box, preserving aspect.
 */
function normalize2(flat) {
  const n = flat.length >> 1;
  let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
  for (let i = 0; i < n; i++) {
    const x = flat[i * 2], y = flat[i * 2 + 1];
    if (x < minX) minX = x;
    if (x > maxX) maxX = x;
    if (y < minY) minY = y;
    if (y > maxY) maxY = y;
  }
  // One scale for both axes — never distort the embedding.
  const span = Math.max(maxX - minX, maxY - minY) || 1;
  const cx = (minX + maxX) / 2, cy = (minY + maxY) / 2;

  const out = new Float32Array(flat.length);
  for (let i = 0; i < n; i++) {
    out[i * 2]     = (flat[i * 2] - cx) / span + 0.5;
    // Flip Y: data-space up should be screen-space up.
    out[i * 2 + 1] = 0.5 - (flat[i * 2 + 1] - cy) / span;
  }
  return out;
}

/**
 * Normalize a flat [x0,y0,z0,…] array into a [-0.5,0.5] cube, preserving
 * aspect across all three axes (so PC variance ratios survive).
 */
function normalize3(flat) {
  const n = flat.length / 3;
  const lo = [Infinity, Infinity, Infinity];
  const hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < n; i++) {
    for (let d = 0; d < 3; d++) {
      const v = flat[i * 3 + d];
      if (v < lo[d]) lo[d] = v;
      if (v > hi[d]) hi[d] = v;
    }
  }
  const span = Math.max(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) || 1;
  const c = [(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2];

  const out = new Float32Array(flat.length);
  for (let i = 0; i < n; i++) {
    out[i * 3]     =  (flat[i * 3]     - c[0]) / span;
    out[i * 3 + 1] = -(flat[i * 3 + 1] - c[1]) / span;   // screen Y is down
    out[i * 3 + 2] =  (flat[i * 3 + 2] - c[2]) / span;
  }
  return out;
}

/**
 * Uniform spatial grid in CSR form, over normalized [0,1] coordinates.
 * Built once per settled layout; queried on every mousemove.
 */
function buildGrid(pos) {
  const n = pos.length >> 1;
  const cellOf = (v) => {
    const c = Math.floor(v * GRID);
    return c < 0 ? 0 : c >= GRID ? GRID - 1 : c;
  };

  const starts = new Int32Array(GRID * GRID + 1);
  const cellIdx = new Int32Array(n);

  for (let i = 0; i < n; i++) {
    const c = cellOf(pos[i * 2 + 1]) * GRID + cellOf(pos[i * 2]);
    cellIdx[i] = c;
    starts[c + 1]++;
  }
  for (let c = 0; c < GRID * GRID; c++) starts[c + 1] += starts[c];

  const cursor = starts.slice(0, GRID * GRID);
  const items = new Int32Array(n);
  for (let i = 0; i < n; i++) items[cursor[cellIdx[i]]++] = i;

  return {
    /** Nearest point to (nx, ny) within maxDist (normalized units), or -1. */
    query(nx, ny, maxDist) {
      const cx = cellOf(nx), cy = cellOf(ny);
      const reach = Math.max(1, Math.ceil(maxDist * GRID));
      let best = -1, bestD2 = maxDist * maxDist;

      for (let r = 0; r <= reach; r++) {
        // Once the ring's nearest possible edge is further than the best hit
        // so far, no further ring can improve on it.
        if (best >= 0 && ((r - 1) / GRID) * ((r - 1) / GRID) > bestD2) break;

        const y0 = Math.max(0, cy - r), y1 = Math.min(GRID - 1, cy + r);
        const x0 = Math.max(0, cx - r), x1 = Math.min(GRID - 1, cx + r);

        for (let gy = y0; gy <= y1; gy++) {
          const onYEdge = gy === cy - r || gy === cy + r;
          for (let gx = x0; gx <= x1; gx++) {
            // Only walk the ring's perimeter; the interior was covered already.
            if (r > 0 && !onYEdge && gx !== cx - r && gx !== cx + r) continue;
            const c = gy * GRID + gx;
            for (let k = starts[c]; k < starts[c + 1]; k++) {
              const i = items[k];
              const dx = pos[i * 2] - nx, dy = pos[i * 2 + 1] - ny;
              const d2 = dx * dx + dy * dy;
              if (d2 < bestD2) { bestD2 = d2; best = i; }
            }
          }
        }
      }
      return best;
    },
  };
}

/* --- renderer ------------------------------------------------------------ */

/**
 * @param {HTMLCanvasElement} canvas
 * @param {object} opts
 *   interactive {boolean}  enable hover hit-testing (default true)
 *   drift       {boolean}  gentle idle motion once at rest
 *   padding     {number}   inset in CSS px (default 18)
 *   onHover     {function} (index|null, {x,y}) => void
 */
export function createScatter(canvas, opts = {}) {
  const ctx = canvas.getContext('2d', { alpha: true });
  const interactive = opts.interactive !== false;
  const drift = !!opts.drift;
  const padding = opts.padding ?? 18;
  const onHover = opts.onHover || (() => {});

  let data = null;
  let posFrom = null, posTo = null, posNow = null;
  let grid = null;
  let morphStart = 0, morphDur = 900, morphing = false;

  let colorBy = null;       // { kind:'categorical'|'continuous', … }
  let active = null;        // Set of visible level indices, or null = all
  let phase = null;         // per-point drift phase

  // 3D intro state
  let intro = null;         // { start, rotateMs, collapseMs, target }
  let depth = null;         // per-point camera depth during the intro
  let dbStarts = null, dbItems = null, dbCursor = null;   // depth×color CSR

  let cssW = 0, cssH = 0, dpr = 1;
  let side = 0, offX = 0, offY = 0;
  let raf = 0, needsDraw = false, running = true;
  let hoverIdx = -1;
  let lastDriftFrame = 0;

  const reduced = prefersReducedMotion();

  /* -- geometry -- */

  function measure() {
    const rect = canvas.getBoundingClientRect();
    if (!rect.width || !rect.height) return false;

    cssW = rect.width;
    cssH = rect.height;
    dpr = Math.min(window.devicePixelRatio || 1, MAX_DPR);

    const bw = Math.round(cssW * dpr), bh = Math.round(cssH * dpr);
    if (canvas.width !== bw || canvas.height !== bh) {
      canvas.width = bw;
      canvas.height = bh;
    }
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    side = Math.max(10, Math.min(cssW, cssH) - padding * 2);
    offX = (cssW - side) / 2;
    offY = (cssH - side) / 2;
    return true;
  }

  function pointRadius() {
    if (!data) return 2;
    // Keep total ink roughly constant as n grows, within sane bounds.
    const r = Math.sqrt((side * side) / data.n) * 0.34;
    return Math.max(1.05, Math.min(r, 3.4));
  }

  /* -- color assignment -- */

  /** Bucket point indices by color so each color is one batched path. */
  function buildBatches() {
    if (!data || !colorBy) return null;
    const n = data.n;

    if (colorBy.kind === 'categorical') {
      const set = colorBy.set;
      const buckets = Array.from({ length: set.levels.length }, () => []);
      const dimmed = [];
      for (let i = 0; i < n; i++) {
        const v = set.values[i];
        if (active && !active.has(v)) dimmed.push(i);
        else if (buckets[v]) buckets[v].push(i);
      }
      return {
        dimmed,
        groups: buckets.map((idx, li) => ({ idx, color: categoricalColor(li) })),
      };
    }

    // Continuous: quantize into CONT_BINS so we still batch draw calls.
    const sc = colorBy.score;
    const [lo, hi] = sc.range;
    const span = (hi - lo) || 1;
    const buckets = Array.from({ length: CONT_BINS }, () => []);
    for (let i = 0; i < n; i++) {
      let t = (sc.values[i] - lo) / span;
      t = t < 0 ? 0 : t > 1 ? 1 : t;
      buckets[Math.min(CONT_BINS - 1, Math.floor(t * CONT_BINS))].push(i);
    }
    return {
      dimmed: [],
      groups: buckets.map((idx, b) => ({
        idx,
        color: viridis((b + 0.5) / CONT_BINS),
      })),
    };
  }

  let batches = null;

  /* -- 3D intro ------------------------------------------------------------
     Rotate the PC1/PC2/PC3 cloud, then collapse it onto the 2D embedding.
     Points are drawn back-to-front in DEPTH_SLICES bands, with size and
     alpha falling off with distance, so the cloud reads as having volume
     without needing a per-frame sort.
     --------------------------------------------------------------------- */

  function projectIntro(now) {
    const c3 = data.coords3d;
    const n = data.n;
    const elapsed = now - intro.start;

    // Rotation eases to a stop as the collapse takes over.
    const spin = Math.min(elapsed, intro.rotateMs + intro.collapseMs);
    const yaw = (spin / 1000) * 0.55;
    const tilt = 0.32;

    const cy = Math.cos(yaw), sy = Math.sin(yaw);
    const cx = Math.cos(tilt), sx = Math.sin(tilt);

    const collapseT = elapsed <= intro.rotateMs ? 0
      : Math.min(1, (elapsed - intro.rotateMs) / intro.collapseMs);
    const e = easeInOutCubic(collapseT);

    const focal = 2.6;
    let dLo = Infinity, dHi = -Infinity;

    for (let i = 0; i < n; i++) {
      const x = c3[i * 3], y = c3[i * 3 + 1], z = c3[i * 3 + 2];

      // yaw about Y, then tilt about X
      const x1 = x * cy + z * sy;
      const z1 = -x * sy + z * cy;
      const y2 = y * cx - z1 * sx;
      const z2 = y * sx + z1 * cx;

      const p = focal / (focal + z2);
      // Perspective flattens out as the cloud collapses into the plane.
      const pp = p + (1 - p) * e;

      const px = 0.5 + x1 * pp;
      const py = 0.5 + y2 * pp;

      posNow[i * 2]     = px + (posTo[i * 2] - px) * e;
      posNow[i * 2 + 1] = py + (posTo[i * 2 + 1] - py) * e;

      const d = z2 * (1 - e);
      depth[i] = d;
      if (d < dLo) dLo = d;
      if (d > dHi) dHi = d;
    }

    intro.depthLo = dLo;
    intro.depthHi = dHi;
    intro.t = e;
    return collapseT >= 1;
  }

  /** Refill the depth × color CSR buckets for this frame. */
  function binByDepth(groupOf, nGroups) {
    const n = data.n;
    const bins = DEPTH_SLICES * nGroups;
    if (!dbStarts || dbStarts.length !== bins + 1) {
      dbStarts = new Int32Array(bins + 1);
      dbItems = new Int32Array(n);
      dbCursor = new Int32Array(bins);
    }
    dbStarts.fill(0);

    const lo = intro.depthLo, span = (intro.depthHi - lo) || 1;
    const slice = new Int32Array(n);
    for (let i = 0; i < n; i++) {
      let s = Math.floor(((depth[i] - lo) / span) * DEPTH_SLICES);
      s = s < 0 ? 0 : s >= DEPTH_SLICES ? DEPTH_SLICES - 1 : s;
      // Larger z is further away, so draw high slices first.
      const b = (DEPTH_SLICES - 1 - s) * nGroups + groupOf(i);
      slice[i] = b;
      dbStarts[b + 1]++;
    }
    for (let b = 0; b < bins; b++) dbStarts[b + 1] += dbStarts[b];
    dbCursor.set(dbStarts.subarray(0, bins));
    for (let i = 0; i < n; i++) dbItems[dbCursor[slice[i]]++] = i;
    return bins;
  }

  function drawIntro() {
    const isCat = colorBy && colorBy.kind === 'categorical';
    const groups = isCat ? colorBy.set.levels.length : 1;
    const groupOf = isCat ? (i) => colorBy.set.values[i] : () => 0;

    binByDepth(groupOf, groups);

    const r = pointRadius();
    const flat = intro.t;                     // 0 = full 3D, 1 = flat

    for (let s = 0; s < DEPTH_SLICES; s++) {
      // Back slices are smaller, dimmer, and cooler — depth cues without a sort.
      const f = s / (DEPTH_SLICES - 1);       // 0 = furthest, 1 = nearest
      const depthAlpha = 0.34 + 0.56 * f;
      const depthScale = 0.72 + 0.46 * f;

      ctx.globalAlpha = depthAlpha + (0.82 - depthAlpha) * flat;
      const rr = r * (depthScale + (1 - depthScale) * flat);

      for (let g = 0; g < groups; g++) {
        const b = s * groups + g;
        const from = dbStarts[b], to = dbStarts[b + 1];
        if (to <= from) continue;
        ctx.fillStyle = isCat ? categoricalColor(g) : '#35B779';
        ctx.beginPath();
        for (let k = from; k < to; k++) {
          const i = dbItems[k];
          const px = offX + posNow[i * 2] * side;
          const py = offY + posNow[i * 2 + 1] * side;
          ctx.moveTo(px + rr, py);
          ctx.arc(px, py, rr, 0, TAU);
        }
        ctx.fill();
      }
    }
    ctx.globalAlpha = 1;
  }

  /* -- drawing -- */

  function strokeGroup(idx, r, driftAmp, t) {
    if (!idx.length) return;
    ctx.beginPath();
    for (let k = 0; k < idx.length; k++) {
      const i = idx[k];
      let x = posNow[i * 2], y = posNow[i * 2 + 1];
      if (driftAmp) {
        const p = phase[i];
        x += Math.sin(t * 0.00035 + p) * driftAmp;
        y += Math.cos(t * 0.00029 + p * 1.7) * driftAmp;
      }
      const px = offX + x * side, py = offY + y * side;
      ctx.moveTo(px + r, py);     // moveTo first, or arcs chain together
      ctx.arc(px, py, r, 0, TAU);
    }
    ctx.fill();
  }

  function draw(now) {
    if (!data || !posNow || !side) return;
    ctx.clearRect(0, 0, cssW, cssH);

    if (intro) { drawIntro(); return; }

    const r = pointRadius();
    const amp = drift && !reduced ? 0.0022 : 0;
    const t = now || 0;

    if (!batches) batches = buildBatches();
    if (!batches) return;

    // Dimmed (deselected) points sit behind, at low alpha.
    if (batches.dimmed.length) {
      ctx.globalAlpha = 0.12;
      ctx.fillStyle = '#8B98A9';
      strokeGroup(batches.dimmed, Math.max(1, r * 0.82), amp, t);
    }

    ctx.globalAlpha = drift ? 0.78 : 0.88;
    for (const g of batches.groups) {
      if (!g.idx.length) continue;
      ctx.fillStyle = g.color;
      strokeGroup(g.idx, r, amp, t);
    }

    // Hovered point: ring it so the tooltip has an obvious anchor.
    if (hoverIdx >= 0) {
      const x = offX + posNow[hoverIdx * 2] * side;
      const y = offY + posNow[hoverIdx * 2 + 1] * side;
      ctx.globalAlpha = 1;
      ctx.beginPath();
      ctx.arc(x, y, r + 3.5, 0, TAU);
      ctx.strokeStyle = '#FDE725';
      ctx.lineWidth = 1.5;
      ctx.stroke();
    }

    ctx.globalAlpha = 1;
  }

  /* -- frame loop -- */

  function frame(now) {
    raf = 0;
    if (!running) return;

    let again = false;

    if (intro) {
      const done = projectIntro(now);
      needsDraw = true;
      if (done) {
        intro = null;
        depth = null;
        posNow.set(posTo);
        grid = interactive ? buildGrid(posNow) : null;
        batches = null;
      }
      again = true;
    } else if (morphing) {
      const p = Math.min(1, (now - morphStart) / morphDur);
      const e = easeInOutCubic(p);
      for (let i = 0; i < posNow.length; i++) {
        posNow[i] = posFrom[i] + (posTo[i] - posFrom[i]) * e;
      }
      if (p >= 1) {
        morphing = false;
        posNow.set(posTo);
        grid = interactive ? buildGrid(posNow) : null;
      } else {
        again = true;
      }
      needsDraw = true;
    }

    if (!intro && drift && !reduced) {
      // 30fps is plenty for idle motion and halves the CPU cost.
      if (now - lastDriftFrame > 33) { lastDriftFrame = now; needsDraw = true; }
      again = true;
    }

    if (needsDraw) { needsDraw = false; draw(now); }
    if (again) raf = requestAnimationFrame(frame);
  }

  function schedule() {
    needsDraw = true;
    if (!raf && running) raf = requestAnimationFrame(frame);
  }

  /* -- hover -- */

  function pick(clientX, clientY) {
    if (!interactive || !grid || !side || morphing || intro) return -1;
    const rect = canvas.getBoundingClientRect();
    const nx = (clientX - rect.left - offX) / side;
    const ny = (clientY - rect.top - offY) / side;
    if (nx < -0.05 || nx > 1.05 || ny < -0.05 || ny > 1.05) return -1;
    // Hit radius in normalized units — a little larger than the dot itself.
    return grid.query(nx, ny, (pointRadius() + 5) / side);
  }

  function handleMove(ev) {
    const idx = pick(ev.clientX, ev.clientY);
    if (idx === hoverIdx) return;
    hoverIdx = idx;
    schedule();
    if (idx < 0) { onHover(null, null); return; }
    onHover(idx, {
      x: offX + posNow[idx * 2] * side,
      y: offY + posNow[idx * 2 + 1] * side,
    });
  }

  function handleLeave() {
    if (hoverIdx === -1) return;
    hoverIdx = -1;
    onHover(null, null);
    schedule();
  }

  if (interactive) {
    canvas.addEventListener('mousemove', handleMove);
    canvas.addEventListener('mouseleave', handleLeave);
    canvas.addEventListener('pointerdown', (ev) => {
      if (ev.pointerType === 'touch') handleMove(ev);
    });
  }

  /* -- resize -- */

  const ro = new ResizeObserver(() => {
    if (measure()) schedule();
  });
  ro.observe(canvas);

  /* -- pause when off-screen: an idle hero canvas should cost nothing -- */

  const io = new IntersectionObserver((entries) => {
    running = entries[0].isIntersecting;
    if (running) { lastDriftFrame = 0; schedule(); }
    else if (raf) { cancelAnimationFrame(raf); raf = 0; }
  }, { threshold: 0 });
  io.observe(canvas);

  /* -- public API -- */

  const api = {
    /** @param {object} parsed  output of parseEmbedding() */
    setData(parsed, initialEmbedding) {
      data = parsed;
      phase = new Float32Array(data.n);
      // Deterministic pseudo-phase — no Math.random, so frames are stable.
      for (let i = 0; i < data.n; i++) phase[i] = (i * 2.399963229728653) % TAU;

      const name = initialEmbedding || data.embeddingNames[0];
      posTo = data.embeddings[name];
      posNow = Float32Array.from(posTo);
      posFrom = Float32Array.from(posTo);
      data.current = name;
      measure();
      grid = interactive ? buildGrid(posNow) : null;
      batches = null;
      schedule();
      return api;
    },

    /**
     * Rotating 3D cloud that collapses onto `target`. No-ops (landing
     * directly on the 2D target) when coords3d is absent or the reader
     * has asked for reduced motion.
     */
    startIntro3D(target, { rotateMs = 1500, collapseMs = 1600 } = {}) {
      if (!data || !data.coords3d) return api;
      const name = data.embeddings[target] ? target : data.embeddingNames[0];
      posTo = data.embeddings[name];
      data.current = name;

      if (reduced) {
        posNow.set(posTo);
        grid = interactive ? buildGrid(posNow) : null;
        schedule();
        return api;
      }

      depth = new Float32Array(data.n);
      intro = { start: performance.now(), rotateMs, collapseMs, t: 0,
                depthLo: -0.5, depthHi: 0.5 };
      grid = null;                 // meaningless mid-flight; rebuilt on settle
      schedule();
      return api;
    },

    setEmbedding(name, { animate = true } = {}) {
      if (!data || !data.embeddings[name] || data.current === name) return api;
      data.current = name;
      posTo = data.embeddings[name];
      if (!animate || reduced || intro) {
        posNow.set(posTo);
        grid = interactive ? buildGrid(posNow) : null;
        schedule();
      } else {
        posFrom.set(posNow);
        morphStart = performance.now();
        morphing = true;
        grid = null;              // stale during flight; rebuilt on settle
        schedule();
      }
      return api;
    },

    setColorBy(spec) { colorBy = spec; batches = null; schedule(); return api; },
    setActiveLevels(set) { active = set; batches = null; schedule(); return api; },
    getHovered() { return hoverIdx; },
    redraw() { measure(); schedule(); return api; },

    destroy() {
      running = false;
      if (raf) cancelAnimationFrame(raf);
      ro.disconnect();
      io.disconnect();
      canvas.removeEventListener('mousemove', handleMove);
      canvas.removeEventListener('mouseleave', handleLeave);
    },
  };

  return api;
}

/* --- data loading -------------------------------------------------------- */

/**
 * Validate and normalize an embedding.json payload.
 * Throws with a readable message rather than failing silently — a malformed
 * export should be obvious during development, not a blank box in production.
 */
export function parseEmbedding(raw) {
  if (!raw || typeof raw !== 'object') throw new Error('payload is not an object');
  if (!raw.coords || typeof raw.coords !== 'object') throw new Error('missing "coords"');

  const names = Object.keys(raw.coords);
  if (!names.length) throw new Error('"coords" contains no embeddings');

  const first = raw.coords[names[0]];
  if (!Array.isArray(first) || first.length < 2 || first.length % 2 !== 0) {
    throw new Error(`"coords.${names[0]}" must be a flat [x,y,…] array of even length`);
  }
  const n = first.length >> 1;

  const embeddings = {};
  for (const name of names) {
    const flat = raw.coords[name];
    if (!Array.isArray(flat) || flat.length !== n * 2) {
      throw new Error(`"coords.${name}" has ${(flat || []).length >> 1} points, expected ${n}`);
    }
    embeddings[name] = normalize2(flat);
  }

  // Optional 3-column cloud used only by the hero intro. Its absence is not
  // an error — the hero simply starts flat.
  let coords3d = null;
  if (Array.isArray(raw.coords3d)) {
    if (raw.coords3d.length !== n * 3) {
      console.warn(`[scatter] coords3d has ${raw.coords3d.length / 3} points, ` +
                   `expected ${n}; skipping the 3D intro.`);
    } else {
      coords3d = normalize3(raw.coords3d);
    }
  }

  const labelSets = {};
  for (const [key, set] of Object.entries(raw.labelSets || {})) {
    if (!Array.isArray(set.levels) || !Array.isArray(set.values)) continue;
    if (set.values.length !== n) {
      throw new Error(`labelSets.${key} has ${set.values.length} values, expected ${n}`);
    }
    const counts = new Array(set.levels.length).fill(0);
    for (const v of set.values) if (counts[v] !== undefined) counts[v]++;
    labelSets[key] = {
      key,
      label: set.label || key,
      levels: set.levels,
      values: Int16Array.from(set.values),
      counts,
    };
  }

  const scores = {};
  for (const [key, sc] of Object.entries(raw.scores || {})) {
    if (!Array.isArray(sc.values)) continue;
    if (sc.values.length !== n) {
      throw new Error(`scores.${key} has ${sc.values.length} values, expected ${n}`);
    }
    let lo = Infinity, hi = -Infinity;
    for (const v of sc.values) { if (v < lo) lo = v; if (v > hi) hi = v; }
    scores[key] = {
      key,
      label: sc.label || key,
      values: Float32Array.from(sc.values),
      range: Array.isArray(sc.range) && sc.range.length === 2 ? sc.range : [lo, hi],
    };
  }

  return {
    n,
    meta: raw.meta || {},
    embeddings,
    embeddingNames: names,
    coords3d,
    labelSets,
    scores,
    current: names[0],
  };
}

/** Fetch + parse, with the HTTP failure surfaced as a readable error. */
export async function loadEmbedding(url) {
  const res = await fetch(url, { cache: 'no-cache' });
  if (!res.ok) throw new Error(`${res.status} ${res.statusText} fetching ${url}`);
  let raw;
  try {
    raw = await res.json();
  } catch (e) {
    throw new Error(`${url} is not valid JSON (${e.message})`);
  }
  return parseEmbedding(raw);
}
