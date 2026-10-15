"""Exploratory figure sketches for collaborator discussion (branch ``figure-ideas``).

Sketches reuse the publication adapters' validated inputs, style, notation and
shared artists, but are deliberately absent from ``registry.py``: they have no
family, are never released, and cannot touch production. Builds write to
``publication/output/sketches/<name>/`` (Git-ignored) with an ``index.json``
describing each sketch. A sketch adopted for the manuscript becomes a family
figure; it is not promoted from here.

    python -m publication.sketches --output-dir publication/output/sketches/NAME [--only KEY ...]
"""
