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
