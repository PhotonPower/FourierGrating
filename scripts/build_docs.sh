#!/usr/bin/env bash
# Builds docs/implementation_plan/main.pdf and docs/handbuch/handbuch.pdf (pdflatex, three passes each;
# makeindex after the first pass of the handbook).
set -euo pipefail
cd "$(dirname "$0")/../docs/implementation_plan"
for i in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error main.tex > build.log || { tail -40 build.log; exit 1; }
done
rm -f main.aux main.log main.out main.toc build.log
echo "written: docs/implementation_plan/main.pdf"

cd ../handbuch
for i in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error handbuch.tex > build.log || { tail -40 build.log; exit 1; }
  if [ "$i" = 1 ]; then makeindex -q -s handbuch.ist handbuch.idx; fi   # Stichwortverzeichnis
done
rm -f handbuch.aux handbuch.log handbuch.out handbuch.toc handbuch.idx handbuch.ilg handbuch.ind build.log
echo "written: docs/handbuch/handbuch.pdf"
