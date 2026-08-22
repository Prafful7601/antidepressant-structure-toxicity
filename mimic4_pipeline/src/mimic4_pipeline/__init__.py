"""MIMIC-IV antidepressant exposure & organ-toxicity reproducibility pipeline.

Tests whether antidepressant organ-system toxicity signals derived from FDA
FAERS spontaneous reports are reproducible in measured MIMIC-IV clinical
data. Built and validated against the open MIMIC-IV Clinical Database Demo
v2.2; designed so that only the database path (see :mod:`mimic4_pipeline.config`)
changes once full, credentialed MIMIC-IV access is available.
"""
