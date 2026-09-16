#!/usr/bin/env python3
"""
make_og_image.py — build the 1200x630 social preview card.

Draws the card as an SVG, using the real UMAP coordinates from
assets/data/embedding.json so the preview shows the same structure the site
does. Rasterize the result with headless Chrome:

    python3 tools/make_og_image.py
    python3 -m http.server 8765 &
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \\
        --headless --disable-gpu --hide-scrollbars \\
        --screenshot=assets/img/og.png --window-size=1200,630 \\
        http://localhost:8765/assets/img/og.svg

(`qlmanage -t` also rasterizes SVG but letterboxes into a square, which
crops the plot panel — use Chrome.)

Re-run both steps whenever the headline copy or the embedding changes.
"""

from __future__ import annotations

import argparse
import json
import os

PALETTE = ["#56B4E9", "#E69F00", "#009E73"]

W, H = 1200, 630
PANEL_X, PANEL_Y, PANEL_SIZE = 700, 80, 470
SUBSAMPLE = 3          # every Nth cell — ~2,000 dots keeps the SVG small

HEADLINE = "Simon Sobti"
EYEBROW = "COMPUTATIONAL BIOLOGIST · PhD, UCSF"
ROLE = "multi-omics · single-cell · machine learning"
PITCH = ["Turning high-dimensional genomics", "into therapeutic targets."]
URL = "ssobt.github.io"

MONO = "ui-monospace,SF Mono,Menlo,monospace"
SANS = "-apple-system,Helvetica Neue,Arial,sans-serif"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", default="assets/data/embedding.json")
    ap.add_argument("--output", default="assets/img/og.svg")
    ap.add_argument("--label", default="perturbation",
                    help="labelSet used to color the dots")
    args = ap.parse_args()

    with open(args.data) as fh:
        d = json.load(fh)

    um = d["coords"]["umap"]
    groups_of = d["labelSets"][args.label]["values"]
    n = len(groups_of)

    xs, ys = um[0::2], um[1::2]
    span = max(max(xs) - min(xs), max(ys) - min(ys))
    cx = (min(xs) + max(xs)) / 2
    cy = (min(ys) + max(ys)) / 2

    buckets: dict[int, list[tuple[float, float]]] = {}
    for i in range(0, n, SUBSAMPLE):
        x = (um[i * 2] - cx) / span + 0.5
        y = 0.5 - (um[i * 2 + 1] - cy) / span
        buckets.setdefault(groups_of[i], []).append(
            (PANEL_X + x * PANEL_SIZE, PANEL_Y + y * PANEL_SIZE))

    layers = []
    for gi in sorted(buckets):
        dots = "".join(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.6"/>'
                       for x, y in buckets[gi])
        layers.append(f'<g fill="{PALETTE[gi % len(PALETTE)]}" opacity="0.85">{dots}</g>')

    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">
<rect width="{W}" height="{H}" fill="#0A0C10"/>
<rect x="0" y="0" width="{W}" height="4" fill="#35B779"/>
{"".join(layers)}
<text x="80" y="214" fill="#35B779" font-family="{MONO}" font-size="21" letter-spacing="3.4">{EYEBROW}</text>
<text x="76" y="316" fill="#E6EDF3" font-family="{SANS}" font-size="82" font-weight="650" letter-spacing="-2.6">{HEADLINE}</text>
<text x="80" y="372" fill="#8B98A9" font-family="{MONO}" font-size="23">{ROLE}</text>
<text x="80" y="446" fill="#8B98A9" font-family="{SANS}" font-size="27">{PITCH[0]}</text>
<text x="80" y="482" fill="#8B98A9" font-family="{SANS}" font-size="27">{PITCH[1]}</text>
<text x="80" y="556" fill="#5B6675" font-family="{MONO}" font-size="19">{URL}</text>
</svg>'''

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as fh:
        fh.write(svg)
    print(f"Wrote {args.output} ({round(len(svg) / 1024)} KB). "
          f"Rasterize it with headless Chrome — see the module docstring.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
