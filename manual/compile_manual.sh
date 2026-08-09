#!/usr/bin/env bash

set -e

cd "$(dirname "$0")"

case "${1:-build}" in
    build)
        echo "Compiling JETSPIN manual..."

        # First LaTeX pass: creates .aux with citation information
        pdflatex -interaction=nonstopmode -halt-on-error manual.tex

        # Build bibliography
        bibtex manual

        # Two further passes resolve bibliography, references and TOC
        pdflatex -interaction=nonstopmode -halt-on-error manual.tex
        pdflatex -interaction=nonstopmode -halt-on-error manual.tex

        echo
        echo "Created: $(pwd)/manual.pdf"
        ;;

    clean)
        rm -f \
            manual.aux \
            manual.bbl \
            manual.blg \
            manual.brf \
            manual.log \
            manual.out \
            manual.toc
        echo "Auxiliary files removed."
        ;;

    distclean)
        rm -f \
            manual.aux \
            manual.bbl \
            manual.blg \
            manual.brf \
            manual.log \
            manual.out \
            manual.toc \
            manual.pdf
        echo "Auxiliary files and manual.pdf removed."
        ;;

    *)
        echo "Usage: $0 [build|clean|distclean]"
        exit 1
        ;;
esac
