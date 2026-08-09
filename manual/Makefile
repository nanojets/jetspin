LATEXMK ?= latexmk
PDFLATEX ?= pdflatex
BIBTEX ?= bibtex

.PHONY: all clean

all: manual.pdf

SECTION_TEX := $(filter-out manual.tex,$(wildcard *.tex))

manual.pdf: manual.tex $(SECTION_TEX) bibliography.bib elsarticle-num.bst figure01.png vmd1.png vmd2.png
	@if command -v $(LATEXMK) >/dev/null 2>&1; then \
		$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error manual.tex; \
	else \
		$(PDFLATEX) -interaction=nonstopmode -halt-on-error manual.tex && \
		$(BIBTEX) manual && \
		$(PDFLATEX) -interaction=nonstopmode -halt-on-error manual.tex && \
		$(PDFLATEX) -interaction=nonstopmode -halt-on-error manual.tex; \
	fi

clean:
	@if command -v $(LATEXMK) >/dev/null 2>&1; then \
		$(LATEXMK) -C manual.tex; \
	else \
		rm -f manual.aux manual.bbl manual.blg manual.log manual.out manual.toc; \
	fi
