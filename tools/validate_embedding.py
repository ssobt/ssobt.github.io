#!/usr/bin/env python3
"""
validate_embedding.py — check an embedding.json against assets/data/SCHEMA.md

Mirrors the checks parseEmbedding() performs in assets/js/scatter.js, so a bad
export fails here at the command line instead of in front of a hiring manager.

    python3 tools/validate_embedding.py assets/data/embedding.json

Exit status 0 = valid, 1 = invalid.
"""

from __future__ import annotations

import json
import os
import sys

MAX_RECOMMENDED_KB = 1024
MAX_RECOMMENDED_LEVELS = 24


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 1
    path = argv[1]

    errors: list[str] = []
    warnings: list[str] = []

    try:
        with open(path) as fh:
            d = json.load(fh)
    except FileNotFoundError:
        print(f"FAIL  no such file: {path}")
        return 1
    except json.JSONDecodeError as e:
        print(f"FAIL  not valid JSON: {e}")
        return 1

    if not isinstance(d, dict):
        print("FAIL  top level must be an object")
        return 1

    # --- coords -------------------------------------------------------------
    coords = d.get("coords")
    if not isinstance(coords, dict) or not coords:
        print('FAIL  missing or empty "coords"')
        return 1

    n = None
    for name, flat in coords.items():
        if not isinstance(flat, list) or len(flat) < 2 or len(flat) % 2:
            errors.append(f'coords.{name}: must be a flat [x,y,...] list of even length')
            continue
        if any(not isinstance(v, (int, float)) for v in flat):
            errors.append(f"coords.{name}: contains non-numeric values")
            continue
        this_n = len(flat) // 2
        if n is None:
            n = this_n
        elif this_n != n:
            errors.append(f"coords.{name}: {this_n} points, expected {n}")

    if n is None:
        print("FAIL  no usable embedding in coords")
        return 1

    if len(coords) < 2:
        warnings.append("only one embedding present — the PCA->UMAP morph needs two")
    if n < 3000:
        warnings.append(f"n = {n}; below ~3,000 the structure tends to look thin")
    if n > 20000:
        warnings.append(f"n = {n}; above ~20,000 the file starts to hurt first paint")

    # --- labelSets ----------------------------------------------------------
    label_sets = d.get("labelSets") or {}
    if not isinstance(label_sets, dict):
        errors.append('"labelSets" must be an object')
        label_sets = {}

    for key, s in label_sets.items():
        levels, values = s.get("levels"), s.get("values")
        if not isinstance(levels, list) or not isinstance(values, list):
            errors.append(f"labelSets.{key}: needs both 'levels' and 'values' lists")
            continue
        if len(values) != n:
            errors.append(f"labelSets.{key}: {len(values)} values, expected {n}")
        bad = [v for v in values if not isinstance(v, int) or v < 0 or v >= len(levels)]
        if bad:
            errors.append(
                f"labelSets.{key}: {len(bad)} value(s) outside 0..{len(levels) - 1} "
                f"(first offender: {bad[0]}) — values must be integer indices into 'levels'")
        if len(levels) > MAX_RECOMMENDED_LEVELS:
            warnings.append(f"labelSets.{key}: {len(levels)} levels makes an unwieldy legend")
        used = set(values)
        unused = [lv for i, lv in enumerate(levels) if i not in used]
        if unused:
            warnings.append(f"labelSets.{key}: unused level(s) {unused} will show n = 0")

    # --- scores -------------------------------------------------------------
    scores = d.get("scores") or {}
    if not isinstance(scores, dict):
        errors.append('"scores" must be an object')
        scores = {}

    for key, s in scores.items():
        values = s.get("values")
        if not isinstance(values, list):
            errors.append(f"scores.{key}: needs a 'values' list")
            continue
        if len(values) != n:
            errors.append(f"scores.{key}: {len(values)} values, expected {n}")
        if any(not isinstance(v, (int, float)) for v in values):
            errors.append(f"scores.{key}: contains non-numeric values")
            continue
        rng = s.get("range")
        if rng is not None:
            if not (isinstance(rng, list) and len(rng) == 2):
                errors.append(f"scores.{key}: 'range' must be [min, max]")
            elif values and (min(values) < rng[0] or max(values) > rng[1]):
                warnings.append(
                    f"scores.{key}: values fall outside the declared range {rng}; "
                    f"points will clamp at the ends of the colourbar")

    if not label_sets and not scores:
        errors.append("no labelSets and no scores — the explorer has nothing to colour by")

    # --- meta ---------------------------------------------------------------
    meta = d.get("meta") or {}
    if not meta.get("caption"):
        warnings.append('meta.caption is empty — the figure will render without a legend')
    if meta.get("n") not in (None, n):
        warnings.append(f"meta.n = {meta.get('n')} disagrees with the actual {n} points")

    # --- size ---------------------------------------------------------------
    kb = round(os.path.getsize(path) / 1024)
    if kb > MAX_RECOMMENDED_KB:
        warnings.append(f"{kb} KB is over the {MAX_RECOMMENDED_KB} KB guideline — lower --n")

    # --- report -------------------------------------------------------------
    for w in warnings:
        print(f"WARN  {w}")
    for e in errors:
        print(f"FAIL  {e}")

    if errors:
        print(f"\n{len(errors)} error(s). See assets/data/SCHEMA.md")
        return 1

    print(f"\nOK    {path}")
    print(f"      {n:,} cells · {kb} KB · embeddings: {', '.join(coords)}")
    if label_sets:
        print(f"      labelSets: {', '.join(f'{k} ({len(v[chr(108)+chr(101)+chr(118)+chr(101)+chr(108)+chr(115)])})' for k, v in label_sets.items())}")
    if scores:
        print(f"      scores: {', '.join(scores)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
