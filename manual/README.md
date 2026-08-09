# JETSPIN user manual

The PDF is kept in the repository for immediate access. Its LaTeX source,
figures, bibliography, and bibliography style are versioned alongside it.

To rebuild the manual locally, run:

```sh
make -C manual
```

The build uses `latexmk` when available. Otherwise it falls back to
`pdflatex` and `bibtex`. Documentation changes are also compiled by the
`Build manual` GitHub Actions workflow, which uploads the resulting PDF as
an artifact.

Run `make -C manual clean` to remove intermediate LaTeX files.

## Alternative build script

The manual can also be compiled with the bundled shell script:

```sh
manual/compile_manual.sh
```

The script runs `pdflatex`, `bibtex`, and the two additional `pdflatex`
passes required to resolve the bibliography, cross-references, and table of
contents. Unlike the Makefile, this method does not attempt to use
`latexmk`, so `pdflatex` and `bibtex` must be available in `PATH`.

To remove auxiliary LaTeX files, run:

```sh
manual/compile_manual.sh clean
```

To remove both auxiliary files and the generated PDF, run:

```sh
manual/compile_manual.sh distclean
```
