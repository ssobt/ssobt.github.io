/* ============================================================================
   explorer.js — controls for the single-cell explorer
   ----------------------------------------------------------------------------
   Loads assets/data/embedding.json once and drives two mounts from it:
   the decorative hero canvas and the interactive figure in §02.

   If the data fails to load, the figure degrades to an explanatory panel
   rather than an empty black box — this page gets read by people deciding
   whether to interview someone, and a broken widget is worse than none.
   ========================================================================= */

import {
  createScatter, loadEmbedding,
  categoricalColor, viridisGradient,
} from './scatter.js';
import { mountCohort, loadCohort } from './cohort.js';

const DATA_URL = new URL('../data/embedding.json', import.meta.url).href;
const COHORT_URL = new URL('../data/cohort.json', import.meta.url).href;

/* --- small DOM helpers --------------------------------------------------- */

const $ = (sel, root = document) => root.querySelector(sel);
const el = (tag, cls, text) => {
  const n = document.createElement(tag);
  if (cls) n.className = cls;
  if (text != null) n.textContent = text;
  return n;
};
const fmt = (n) => n.toLocaleString('en-US');

/* --- hero ---------------------------------------------------------------- */

function mountHero(data) {
  const canvas = $('#hero-canvas');
  if (!canvas) return;

  const names = data.embeddingNames;
  const rest = names.includes('umap') ? 'umap' : names[names.length - 1];

  const scatter = createScatter(canvas, {
    interactive: false,
    drift: true,
    padding: 26,
  });

  // Color by the first categorical label set, so the hero is already
  // showing real structure rather than an undifferentiated cloud.
  const firstSet = Object.values(data.labelSets)[0];
  scatter.setData(data, rest);
  if (firstSet) scatter.setColorBy({ kind: 'categorical', set: firstSet });

  // The page's opening gesture: a rotating PC1/PC2/PC3 cloud that collapses
  // into the 2D embedding — dimensionality reduction, performed. Falls back
  // to a static 2D figure when coords3d is absent or motion is reduced.
  scatter.startIntro3D(rest, { rotateMs: 1500, collapseMs: 1700 });

  const cap = $('#hero-figcap-text');
  if (cap) {
    cap.textContent = data.meta.heroCaption
      || `${fmt(data.n)} cells · ${rest.toUpperCase()}`;
  }
}

/* --- explorer ------------------------------------------------------------ */

function mountExplorer(data) {
  const root = $('#explorer');
  if (!root) return;

  const canvas = $('#explorer-canvas', root);
  const tip = $('#explorer-tip', root);
  const legendEl = $('#explorer-legend', root);
  const colorbar = $('#explorer-colorbar', root);
  const embedSeg = $('#explorer-embed', root);
  const colorSeg = $('#explorer-color', root);
  const tableWrap = $('#explorer-table');        // lives outside #explorer
  const capEl = $('#explorer-caption', root);

  /* ---- color modes: every label set, then every continuous score ---- */

  const modes = [
    ...Object.values(data.labelSets).map((set) => ({
      key: set.key, label: set.label, kind: 'categorical', set,
    })),
    ...Object.values(data.scores).map((sc) => ({
      key: sc.key, label: sc.label, kind: 'continuous', score: sc,
    })),
  ];
  if (!modes.length) {
    fail(root, 'The embedding contains no label sets or scores to color by.');
    return;
  }

  let mode = modes[0];
  let active = null;          // Set of visible level indices, or null = all

  const scatter = createScatter(canvas, {
    interactive: true,
    padding: 22,
    onHover: showTip,
  });

  scatter.setData(data, data.embeddingNames[0]);

  /* ---- tooltip ---- */

  function showTip(idx, pos) {
    if (idx == null || !pos) { tip.dataset.on = 'false'; return; }

    tip.replaceChildren();
    for (const set of Object.values(data.labelSets)) {
      const row = el('div');
      row.append(el('dt', null, `${set.label}:`));
      row.append(el('dd', null, set.levels[set.values[idx]] ?? '—'));
      tip.append(row);
    }
    for (const sc of Object.values(data.scores)) {
      const row = el('div');
      row.append(el('dt', null, `${sc.label}:`));
      row.append(el('dd', null, sc.values[idx].toFixed(3)));
      tip.append(row);
    }

    tip.style.left = `${pos.x}px`;
    tip.style.top = `${pos.y}px`;
    tip.dataset.on = 'true';
  }

  /* ---- legend ---- */

  function renderLegend() {
    legendEl.replaceChildren();

    if (mode.kind !== 'categorical') {
      colorbar.dataset.on = 'true';
      $('.colorbar__ramp', colorbar).style.background = viridisGradient();
      const [lo, hi] = mode.score.range;
      const ticks = $('.colorbar__ticks', colorbar);
      ticks.replaceChildren(
        el('span', null, lo.toFixed(2)),
        el('span', null, mode.score.label),
        el('span', null, hi.toFixed(2)),
      );
      return;
    }

    colorbar.dataset.on = 'false';
    const { set } = mode;

    set.levels.forEach((name, i) => {
      const li = el('li');
      const btn = el('button', 'legend__item');
      btn.type = 'button';
      btn.setAttribute('aria-pressed', String(!active || active.has(i)));

      const sw = el('span', 'legend__swatch');
      sw.style.background = categoricalColor(i);
      btn.append(sw, el('span', 'legend__name', name),
                 el('span', 'legend__n', fmt(set.counts[i])));

      btn.addEventListener('click', () => { toggleLevel(i); });
      li.append(btn);
      legendEl.append(li);
    });
  }

  function toggleLevel(i) {
    if (active === null) {
      active = new Set([i]);                      // first click isolates
    } else if (active.has(i)) {
      if (active.size === 1) active = null;       // clicking the lone one resets
      else { active.delete(i); if (!active.size) active = null; }
    } else {
      active.add(i);
    }
    scatter.setActiveLevels(active);
    renderLegend();
    renderTable();
  }

  /* ---- accessible mirror of the canvas ---- */

  function renderTable() {
    tableWrap.replaceChildren();
    if (mode.kind !== 'categorical') {
      const [lo, hi] = mode.score.range;
      tableWrap.append(el('p', 'mono',
        `${mode.label}: ${fmt(data.n)} cells, range ${lo.toFixed(3)} – ${hi.toFixed(3)}, ` +
        `colored on a viridis scale from low (dark purple) to high (yellow).`));
      return;
    }

    const { set } = mode;
    const table = el('table', 'tbl');
    const thead = el('thead');
    const hr = el('tr');
    ['Group', 'Cells', '% of total', 'Shown'].forEach((h, i) => {
      const th = el('th', i > 0 ? 'num' : null, h);
      th.scope = 'col';
      hr.append(th);
    });
    thead.append(hr);

    const tbody = el('tbody');
    set.levels.forEach((name, i) => {
      const tr = el('tr');
      const shown = !active || active.has(i);
      tr.append(el('td', null, name));
      tr.append(el('td', 'num', fmt(set.counts[i])));
      tr.append(el('td', 'num', `${((set.counts[i] / data.n) * 100).toFixed(1)}%`));
      tr.append(el('td', 'num', shown ? 'yes' : 'no'));
      tbody.append(tr);
    });

    table.append(thead, tbody);
    const caption = el('caption', 'visually-hidden',
      `Cell counts per ${set.label} group, and whether each is currently shown in the plot.`);
    table.prepend(caption);

    const wrap = el('div', 'tbl-wrap');
    wrap.append(table);
    tableWrap.append(wrap);
  }

  /* ---- segmented controls ---- */

  function buildSeg(container, items, current, onPick) {
    container.replaceChildren();
    items.forEach((it) => {
      const b = el('button', null, it.label);
      b.type = 'button';
      b.dataset.key = it.key;
      b.setAttribute('aria-pressed', String(it.key === current));
      b.addEventListener('click', () => onPick(it));
      container.append(b);
    });
  }


  const embedItems = data.embeddingNames.map((name) => ({
    key: name,
    label: (data.meta.embeddingLabels && data.meta.embeddingLabels[name]) || name.toUpperCase(),
  }));

  function selectEmbedding(it) {
    scatter.setEmbedding(it.key, { animate: true });
    tip.dataset.on = 'false';
    buildSeg(embedSeg, embedItems, it.key, selectEmbedding);
  }

  function selectMode(it) {
    mode = it;
    active = null;
    scatter.setActiveLevels(null);
    scatter.setColorBy(
      it.kind === 'categorical'
        ? { kind: 'categorical', set: it.set }
        : { kind: 'continuous', score: it.score },
    );
    buildSeg(colorSeg, modes, it.key, selectMode);
    renderLegend();
    renderTable();
  }

  buildSeg(embedSeg, embedItems, data.embeddingNames[0], selectEmbedding);
  selectMode(mode);

  /* ---- caption ---- */

  if (capEl && data.meta.caption) {
    capEl.replaceChildren(
      el('b', null, 'Figure 1. '),
      document.createTextNode(data.meta.caption),
    );
  }
}

/* --- failure state ------------------------------------------------------- */

function fail(root, message) {
  const status = $('#explorer-status', root) || $('#explorer-status');
  if (!status) return;
  status.dataset.on = 'true';
  status.replaceChildren(el('p', null, message));
}

/* --- boot ---------------------------------------------------------------- */

/* The two figures are independent: one failing must not blank the other. */
(async function bootCohort() {
  const root = $('#cohort');
  if (!root) return;
  const status = $('#cohort-status', root);
  try {
    const data = await loadCohort(COHORT_URL);
    if (status) status.dataset.on = 'false';
    mountCohort(data, root);
  } catch (err) {
    console.error('[cohort] could not load:', err);
    if (status) {
      status.dataset.on = 'true';
      status.replaceChildren(el('p', null,
        'The cohort figure could not load its data. The analysis it shows is ' +
        'described in the text above and in the linked preprint.'));
    }
  }
})();

(async function boot() {
  const status = $('#explorer-status');

  let data;
  try {
    data = await loadEmbedding(DATA_URL);
  } catch (err) {
    console.error('[explorer] could not load embedding:', err);
    fail(document,
      'The interactive figure could not load its data. The analysis it renders is ' +
      'described in the caption below, and the code that produced it is linked from ' +
      'the project cards above.');
    const cap = $('#hero-figcap-text');
    if (cap) cap.textContent = 'figure unavailable';
    return;
  }

  if (status) status.dataset.on = 'false';

  // Mounted independently: a fault in one must not take out the other.
  try { mountHero(data); }
  catch (err) { console.error('[explorer] hero failed:', err); }

  try { mountExplorer(data); }
  catch (err) {
    console.error('[explorer] figure failed:', err);
    fail(document, 'The interactive figure could not be drawn. Its caption and the ' +
                   'equivalent data table below still describe the analysis.');
  }
})();
