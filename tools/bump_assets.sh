#!/usr/bin/env bash
# Re-stamp the ?v= query on every stylesheet and entry-point script.
#
# index.html and site.css are cached independently by GitHub Pages
# (cache-control: max-age=600). Without this stamp a deploy that changes both
# can serve new HTML against a visitor's 10-minute-old CSS, which renders the
# page wrong until their cache expires. Run this after any CSS/JS change.
set -euo pipefail
cd "$(dirname "$0")/.."
V=$(date +%Y%m%d%H%M)
for f in index.html projects/*.html; do
  sed -i '' -E "s|(href=\"[^\"]*assets/css/site\.css)(\?v=[0-9]+)?\"|\1?v=$V\"|g" "$f"
  sed -i '' -E "s|(src=\"[^\"]*assets/js/[a-z]+\.js)(\?v=[0-9]+)?\"|\1?v=$V\"|g" "$f"
done

# Entry-point scripts alone are not enough: a module imported by URL is cached
# under that URL, so bare './cohort.js' would still come from cache. Stamp the
# import specifiers too.
for f in assets/js/*.js; do
  sed -i '' -E "s|(from '\./[a-z]+\.js)(\?v=[0-9]+)?'|\1?v=$V'|g" "$f"
done
echo "stamped assets ?v=$V"
