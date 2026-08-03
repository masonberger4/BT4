"""Packaged codon-usage and tRNA data (``*.tsv`` + ``*.provenance.json``).

This package is intentionally code-free: it exists so the data directory is a
*regular* package rather than an implicit namespace package. Regular packages
resolve reliably under ``importlib.resources.files(...)`` in frozen builds
(PyInstaller one-file **and** macOS ``.app`` bundles), which is how the codon and
tAI loaders read these files at runtime.
"""
