"""Free-text dose parsing and the low-dose stratification rule for task 3.

``dose_val_rx`` in ``hosp.prescriptions`` is free text, not a clean number
(plain values like ``"50"``, ranges like ``"25-50"``, and unparseable
strings all occur). This is fuzzy-input logic in the same spirit as the
drug-name matcher, so it gets its own tested module rather than living
inline in the cohort SQL.
"""

from __future__ import annotations

import re

# Daily-dose thresholds (mg/day) below which amitriptyline, doxepin, and
# trazodone are conventionally non-antidepressant (analgesic/hypnotic) use
# rather than depression treatment (brief, task 3). Antidepressant dosing
# for the two TCAs conventionally starts at 75 mg/day; trazodone's hypnotic
# range (25-100 mg) is well below its antidepressant range (150-400 mg).
# These are literature-convention cutoffs, not derived from this dataset --
# adjust here if a different convention is preferred; every downstream user
# reads the threshold from this one dict.
DOSE_THRESHOLDS_MG_PER_DAY: dict[str, float] = {
    "amitriptyline": 75.0,
    "doxepin": 75.0,
    "trazodone": 150.0,
}

_NUM_RE = re.compile(r"[0-9]*\.?[0-9]+")


def parse_dose_mg(dose_val_rx: str | float | None) -> float | None:
    """Parse one ``dose_val_rx`` value into a numeric mg amount.

    - Plain numbers ("50", "7.5") parse directly.
    - Hyphenated ranges ("25-50") take the **upper** bound: the conservative
      choice for a low-dose screen, since it never under-counts a dose that
      might actually reach the standard/antidepressant range.
    - Anything with no parseable number (blank, "QID", None, NaN) -> None.
    """
    if dose_val_rx is None:
        return None
    s = str(dose_val_rx).strip()
    if not s or s.lower() == "nan":
        return None
    nums = [float(n) for n in _NUM_RE.findall(s)]
    if not nums:
        return None
    return max(nums)


def estimate_daily_dose_mg(dose_val_rx, doses_per_24_hrs) -> float | None:
    """Estimate total daily dose in mg from one prescription row.

    ``doses_per_24_hrs`` is frequently null in MIMIC-IV; missing or
    non-positive values are treated as a single administration (the
    conservative floor -- this can only under-, never over-, estimate daily
    dose, consistent with parse_dose_mg's own conservative range handling).
    """
    mg = parse_dose_mg(dose_val_rx)
    if mg is None:
        return None
    n = doses_per_24_hrs
    if n is None or (isinstance(n, float) and n != n) or n <= 0:  # n != n catches NaN
        n = 1.0
    return mg * float(n)


def class_is_dose_ambiguous_low(agents_in_class: list[str], dose_stratum_by_agent: dict[str, str]) -> bool:
    """True iff a drug class's presence in one admission is fully explained
    by low-dose, presumptively non-psychiatric use of dose-ambiguous agents
    -- i.e. every agent contributing to that class in this admission is one
    of the three dose-ambiguous agents (amitriptyline/doxepin/trazodone)
    AND is classified ``"low_dose"``.

    This is the test behind the primary multi-class resolution rule (see
    :mod:`mimic4_pipeline.cohort`): a class failing this test (because it
    has a non-dose-ambiguous agent like bupropion, or a dose-ambiguous
    agent at standard or unknown dose) is never dropped -- only low-dose
    dose-ambiguous-only classes are reclassified away.
    """
    if not agents_in_class:
        return False
    for agent in agents_in_class:
        if agent not in DOSE_THRESHOLDS_MG_PER_DAY:
            return False
        if dose_stratum_by_agent.get(agent) != "low_dose":
            return False
    return True


def classify_dose_stratum(agent: str, daily_doses_mg: list[float | None]) -> str:
    """Classify one admission's exposure to a dose-ambiguous agent.

    ``daily_doses_mg`` is every estimated daily dose observed for this
    agent during the admission (dose can change over a stay -- the
    admission-level classification uses the *maximum*, i.e. the strongest
    evidence of standard/antidepressant-range dosing).

    Returns one of ``"not_applicable"`` (agent isn't dose-ambiguous),
    ``"low_dose"``, ``"standard_dose"``, or ``"unknown_dose"`` (no dose
    could be parsed for any administration).
    """
    threshold = DOSE_THRESHOLDS_MG_PER_DAY.get(agent)
    if threshold is None:
        return "not_applicable"
    # A value that failed to parse comes back as None from
    # estimate_daily_dose_mg, but pandas silently turns None into float
    # NaN inside a float column (e.g. after a groupby/.tolist() round
    # trip) -- `nan is not None` is True, so filtering on identity alone
    # would let an unparseable dose through as if it were real, and
    # `max()` with a NaN in the list is likewise unreliable. Both must be
    # excluded explicitly.
    resolved = [d for d in daily_doses_mg if d is not None and not (isinstance(d, float) and d != d)]
    if not resolved:
        return "unknown_dose"
    return "standard_dose" if max(resolved) >= threshold else "low_dose"
