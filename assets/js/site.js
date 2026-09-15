/* ============================================================================
   site.js — page chrome. Scrollspy, sticky-bar state, reveal-on-scroll,
   copy-to-clipboard for the email address.

   Everything here is an enhancement: with JS disabled the page is fully
   readable, every link works, and nothing is hidden.
   ========================================================================= */

/* --- sticky bar gets a border once you've left the top ------------------- */

const bar = document.querySelector('.topbar');
if (bar) {
  const onScroll = () => {
    bar.dataset.scrolled = String(window.scrollY > 8);
  };
  addEventListener('scroll', onScroll, { passive: true });
  onScroll();
}

/* --- scrollspy ----------------------------------------------------------- */

const navLinks = [...document.querySelectorAll('.topbar__nav a[href^="#"]')];
const sections = navLinks
  .map((a) => document.querySelector(a.getAttribute('href')))
  .filter(Boolean);

if (sections.length) {
  const setCurrent = (id) => {
    navLinks.forEach((a) => {
      a.setAttribute('aria-current', String(a.getAttribute('href') === `#${id}`));
    });
  };

  const spy = new IntersectionObserver((entries) => {
    // Pick the entry nearest the top of the viewport that's currently visible.
    const visible = entries
      .filter((e) => e.isIntersecting)
      .sort((a, b) => a.boundingClientRect.top - b.boundingClientRect.top);
    if (visible.length) setCurrent(visible[0].target.id);
  }, { rootMargin: '-15% 0px -70% 0px', threshold: 0 });

  sections.forEach((s) => spy.observe(s));
}

/* --- reveal on scroll ---------------------------------------------------
   The reveal is decorative, but it works by hiding content until an observer
   fires — so every failure mode has to end with the content visible. Hence
   the immediate first pass and the unconditional timeout below.
   ------------------------------------------------------------------------ */

const revealables = [...document.querySelectorAll('.reveal')];
const reduced = matchMedia('(prefers-reduced-motion: reduce)').matches;
const showAll = () => revealables.forEach((n) => { n.dataset.shown = 'true'; });

if (reduced || !('IntersectionObserver' in window)) {
  showAll();
} else {
  const io = new IntersectionObserver((entries, obs) => {
    entries.forEach((e) => {
      if (!e.isIntersecting) return;
      e.target.dataset.shown = 'true';
      obs.unobserve(e.target);
    });
  }, { rootMargin: '0px 0px -8% 0px', threshold: 0.05 });

  revealables.forEach((n) => {
    // Anything already on screen at load should not wait for a callback.
    const r = n.getBoundingClientRect();
    if (r.top < innerHeight && r.bottom > 0) n.dataset.shown = 'true';
    else io.observe(n);
  });

  // Failsafe: if the observer never delivers — odd embedding contexts,
  // prerender, a throttled tab — show everything anyway. Invisible content
  // is a far worse outcome than a missed animation.
  setTimeout(showAll, 2500);
}

/* --- copy email ---------------------------------------------------------- */

document.querySelectorAll('[data-copy]').forEach((btn) => {
  btn.addEventListener('click', async (ev) => {
    const value = btn.dataset.copy;
    if (!navigator.clipboard) return;          // let the mailto: link do its job
    ev.preventDefault();
    try {
      await navigator.clipboard.writeText(value);
      const label = btn.querySelector('.v') || btn;
      const original = label.textContent;
      label.textContent = 'copied ✓';
      setTimeout(() => { label.textContent = original; }, 1400);
    } catch {
      location.href = `mailto:${value}`;
    }
  });
});

/* --- current year in the footer ------------------------------------------ */

const yearEl = document.querySelector('[data-year]');
if (yearEl) yearEl.textContent = String(new Date().getFullYear());
