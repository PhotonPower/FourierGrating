#!/usr/bin/env bash
# Builds docs/implementation_plan/main.pdf (pdflatex, three passes for references and TOC).
set -euo pipefail
cd "$(dirname "$0")/../docs/implementation_plan"
for i in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error main.tex > build.log || { tail -40 build.log; exit 1; }
done
rm -f main.aux main.log main.out main.toc build.log
echo "written: docs/implementation_plan/main.pdf"
