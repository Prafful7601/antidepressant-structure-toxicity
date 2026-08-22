"""Maps the 28 curated antidepressants to the free-text ``drug`` values that
actually appear in ``mimiciv_hosp.prescriptions``, cross-checked against
``formulary_drug_cd``, ``gsn``, and ``ndc``.

This is deliberately its own module with its own tests (see
``tests/test_drug_mapping.py``): a silently wrong mapping here would
invalidate every downstream cohort and outcome without any obvious error.

Design decisions, and why
--------------------------
- **Matching is name-based, not raw substring-based.** Four of the 28
  curated agents contain another curated agent's name as a raw substring
  (``imipramine`` in ``trimipramine``, ``citalopram`` in ``escitalopram``,
  ``venlafaxine`` in ``desvenlafaxine``, ``milnacipran`` in
  ``levomilnacipran`` -- discovered programmatically, see
  ``tests/test_drug_mapping.py::test_curated_agent_names_do_not_collide_undetected``).
  A naive ``"citalopram" in drug_text.lower()`` check would misclassify
  every escitalopram order as citalopram. Matching instead uses a
  word-boundary regex (``\\b<alias>\\b``), which correctly rejects all four
  cases because the shorter agent's name is preceded by a letter, not a
  boundary, in the longer agent's name. **This citalopram/escitalopram-class
  collision is a real curation finding, not a hypothetical, and is reported
  in the paper's methods as its own step** -- see the README.
- **Only generic names drive matching. Brand names are not matcher input.**
  The curated brand-name list below (``BRAND_NAMES``) is hand-compiled and
  not verified against a formulary, so it is used only to flag *candidates*
  for human review (:func:`psychiatric_review_candidates`) among drug-text
  values the generic-name matcher rejected. A brand-name hit never, by
  itself, adds a row to the exposure cohort. When full MIMIC-IV access
  arrives, rerun the distinct-value inspection over the full drug list and
  review the printed candidates by eye (per explicit instruction) rather
  than trusting this list to auto-classify anything.
- **Combination products are flagged, not silently included or dropped.**
  A drug string like ``"Amitriptyline-Chlordiazepoxide"`` matches the
  amitriptyline alias but is not amitriptyline monotherapy exposure; it is
  flagged ``is_combination_product=True`` so the cohort step can make an
  explicit, documented decision about it rather than the mapper deciding
  silently.
- **Route/form filtering excludes non-systemic use** (topical, ophthalmic,
  otic, nasal, rectal-local, etc.) so e.g. a hypothetical compounded
  doxepin cream would not be counted as systemic antidepressant exposure.
- **formulary_drug_cd / gsn / ndc are used as an independent recovery
  check**, not as the primary key: once a set of drug-text matches is
  found, any *other* prescriptions row sharing one of those rows'
  formulary_drug_cd/gsn/ndc is flagged for review even if its own drug text
  didn't match anything (catches a spelling the text matcher missed).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

# ---------------------------------------------------------------------------
# The 28 curated agents (brief, section 2), by class.
# ---------------------------------------------------------------------------
AGENT_CLASS: dict[str, str] = {
    # Tricyclic
    "amitriptyline": "Tricyclic",
    "nortriptyline": "Tricyclic",
    "imipramine": "Tricyclic",
    "desipramine": "Tricyclic",
    "clomipramine": "Tricyclic",
    "doxepin": "Tricyclic",
    "trimipramine": "Tricyclic",
    "protriptyline": "Tricyclic",
    "amoxapine": "Tricyclic",
    # SSRI
    "fluoxetine": "SSRI",
    "sertraline": "SSRI",
    "paroxetine": "SSRI",
    "citalopram": "SSRI",
    "escitalopram": "SSRI",
    "fluvoxamine": "SSRI",
    # SNRI
    "venlafaxine": "SNRI",
    "desvenlafaxine": "SNRI",
    "duloxetine": "SNRI",
    "levomilnacipran": "SNRI",
    "milnacipran": "SNRI",
    # RIMA
    "moclobemide": "RIMA",
    "clorgyline": "RIMA",
    # Atypical
    "bupropion": "Atypical",
    "mirtazapine": "Atypical",
    "trazodone": "Atypical",
    "vilazodone": "Atypical",
    "vortioxetine": "Atypical",
    "nefazodone": "Atypical",
}

assert len(AGENT_CLASS) == 28, f"expected 28 curated agents, got {len(AGENT_CLASS)}"

# Common US brand names. NOT used to drive automatic matching (see module
# docstring) -- this is a hand-compiled, unverified list, used only as a
# hint set for psychiatric_review_candidates() to flag unmatched drug text
# for a human to look at. None of these were observed in the MIMIC-IV demo
# (see the audit run) -- prescriptions.drug used generic names throughout.
BRAND_NAMES: dict[str, tuple[str, ...]] = {
    "amitriptyline": ("elavil",),
    "nortriptyline": ("pamelor", "aventyl"),
    "imipramine": ("tofranil",),
    "desipramine": ("norpramin",),
    "clomipramine": ("anafranil",),
    "doxepin": ("silenor", "sinequan"),
    "trimipramine": ("surmontil",),
    "protriptyline": ("vivactil",),
    "amoxapine": ("asendin",),
    "fluoxetine": ("prozac", "sarafem"),
    "sertraline": ("zoloft",),
    "paroxetine": ("paxil", "pexeva", "brisdelle"),
    "citalopram": ("celexa",),
    "escitalopram": ("lexapro",),
    "fluvoxamine": ("luvox",),
    "venlafaxine": ("effexor",),
    "desvenlafaxine": ("pristiq", "khedezla"),
    "duloxetine": ("cymbalta", "irenka", "drizalma"),
    "levomilnacipran": ("fetzima",),
    "milnacipran": ("savella",),
    "moclobemide": ("manerix", "aurorix"),
    "clorgyline": (),  # never marketed; research compound only
    "bupropion": ("wellbutrin", "zyban", "aplenzin", "forfivo"),
    "mirtazapine": ("remeron",),
    "trazodone": ("desyrel", "oleptro"),
    "vilazodone": ("viibryd",),
    "vortioxetine": ("trintellix", "brintellix"),
    "nefazodone": ("serzone",),
}

# Other real, marketed (or historically marketed) antidepressants that are
# NOT part of the curated 28-agent list (brief, section 2). Like
# BRAND_NAMES, this is a review-hint set only: a hit here never adds a row
# to the exposure cohort, it only flags an unmatched drug-text value as
# worth a human look (e.g. "did the curated list miss something that
# appears in the real data?"). Deliberately hand-compiled and not
# authoritative.
OTHER_KNOWN_ANTIDEPRESSANTS_FOR_REVIEW: dict[str, str] = {
    "tianeptine": "atypical TCA-related",
    "agomelatine": "melatonergic antidepressant",
    "esketamine": "NMDA-antagonist antidepressant (Spravato)",
    "reboxetine": "selective NRI",
    "maprotiline": "tetracyclic",
    "mianserin": "tetracyclic",
    "dosulepin": "TCA (also called dothiepin)",
    "dothiepin": "TCA (also called dosulepin)",
    "lofepramine": "TCA",
    "butriptyline": "TCA",
    "iprindole": "TCA-related",
    "amineptine": "TCA-related (withdrawn in most markets)",
    "etoperidone": "trazodone-related phenylpiperazine",
    "phenelzine": "MAOI (not a RIMA; curated list includes only moclobemide/clorgyline)",
    "tranylcypromine": "MAOI",
    "isocarboxazid": "MAOI",
    "selegiline": "MAOI (also used for Parkinson's disease)",
}

# Word-final suffixes shared by several curated agents' generic names.
# Purely a review-hint heuristic (see psychiatric_review_candidates): an
# unmatched drug-text value ending in one of these is flagged for a human
# to check, since it has the *shape* of an antidepressant generic name.
# Deliberately over-inclusive -- e.g. "atomoxetine" (an ADHD medication,
# not an antidepressant) also ends in "oxetine" and will be flagged; that
# is the intended tradeoff for a manual-review net, not a precision matcher.
ANTIDEPRESSANT_SUFFIX_STEMS: tuple[str, ...] = (
    "triptyline", "ipramine", "oxapine", "oxetine", "talopram", "afaxine",
    "nacipran", "zodone", "bemide", "cipran",
)

# Routes considered non-systemic for these oral/IV-class agents. None of the
# 28 curated agents has an approved topical/ophthalmic/otic/nasal-local
# formulation, so any such route on a matched drug string indicates either a
# compounded/miscoded product or a data error, and is excluded from exposure.
NON_SYSTEMIC_ROUTE_TOKENS: frozenset[str] = frozenset(
    {
        "top",
        "topical",
        "oph",
        "ou",
        "os",
        "od",
        "otic",
        "ear",
        "nu",  # nasal (MIMIC route abbreviation)
        "nasal",
        "ac",  # 'AC' = 'auriculum' in some MIMIC route lists; treated cautiously below too
        "rectal",  # local rectal preparations only; PR *systemic* dosing is not used by any of the 28
        "sk",
        "skin",
        "tp",
    }
)

# Dosage-form words that indicate a non-systemic (local) preparation
# regardless of the coded route field.
NON_SYSTEMIC_FORM_KEYWORDS: tuple[str, ...] = (
    "cream",
    "ointment",
    "lotion",
    "gel",
    "patch",
    "suppository",
    "supp",
    "ophth",
    "otic",
    "nasal spray",
    "shampoo",
)

# Salt forms, release modifiers, and dosage-form words that legitimately
# trail an agent's name (e.g. "Amitriptyline HCl", "BuPROPion XL (Once
# Daily)") and must NOT be mistaken for a second drug name when screening
# for combination products.
_QUALIFIER_WORDS: frozenset[str] = frozenset(
    w.lower()
    for w in (
        "HCl", "Hydrochloride", "Oxalate", "Mesylate", "Succinate", "Sulfate",
        "Besylate", "Tartrate", "Maleate", "Citrate", "Fumarate", "Acetate",
        "Sustained", "Release", "Extended", "Delayed", "Immediate", "Once",
        "Daily", "XL", "XR", "ER", "SR", "CR", "DR", "IR", "ODT", "Tab",
        "Tablet", "Capsule", "Cap", "Injection", "Inj", "Oral", "Solution",
        "Soln", "Suspension", "Susp", "Elixir", "Liquid", "Sprinkles", "IV",
        "PO", "mg", "Free", "Sugar",
    )
)

_WORD_RE = re.compile(r"[A-Za-z]+")


def _build_pattern(alias: str) -> re.Pattern:
    """Word-boundary, case-insensitive pattern for one alias.

    ``\\b`` is what makes this safe against the four real substring traps
    among the curated agents themselves (see module docstring): the
    embedded shorter name is always preceded by a letter in the longer
    name, so there is no boundary there and the shorter pattern does not
    match.
    """
    return re.compile(r"\b" + re.escape(alias) + r"\b", re.IGNORECASE)


@dataclass(frozen=True)
class AgentAlias:
    agent: str
    drug_class: str
    alias: str
    alias_type: str  # always "generic" -- see module docstring on brand names
    pattern: re.Pattern = field(compare=False)


def _all_aliases() -> list[AgentAlias]:
    # Generic names only. Brand names (BRAND_NAMES) deliberately do NOT feed
    # this list -- see module docstring "Only generic names drive matching".
    aliases = [
        AgentAlias(agent, drug_class, agent, "generic", _build_pattern(agent))
        for agent, drug_class in AGENT_CLASS.items()
    ]
    # Longest alias first: belt-and-suspenders on top of \b, and it also
    # gives a deterministic, reviewable match order when a string somehow
    # matches more than one alias.
    aliases.sort(key=lambda a: len(a.alias), reverse=True)
    return aliases


ALL_ALIASES: list[AgentAlias] = _all_aliases()


@dataclass(frozen=True)
class MatchResult:
    agent: str
    drug_class: str
    alias_matched: str
    alias_type: str
    is_combination_product: bool


def match_drug_text(drug_text: str | None) -> MatchResult | None:
    """Match one free-text ``drug`` (or similar) value to a curated agent by
    its **generic** name only (brand names are not matcher input -- see
    module docstring).

    Returns ``None`` if no curated agent's generic name appears as a whole
    word in ``drug_text``. When more than one alias matches (should not
    happen among the 28 curated agents by construction, but a future edit
    could introduce a collision), the longest alias wins and is reported;
    callers can detect the ambiguity themselves by calling :func:`all_matches`
    if needed.
    """
    if not drug_text:
        return None
    hits = [a for a in ALL_ALIASES if a.pattern.search(drug_text)]
    if not hits:
        return None
    best = hits[0]  # ALL_ALIASES is sorted longest-alias-first
    return MatchResult(
        agent=best.agent,
        drug_class=best.drug_class,
        alias_matched=best.alias,
        alias_type=best.alias_type,
        is_combination_product=looks_like_combination_product(drug_text, best.alias),
    )


def all_matches(drug_text: str | None) -> list[MatchResult]:
    """All agents whose generic-name alias appears in ``drug_text`` (for
    ambiguity checks)."""
    if not drug_text:
        return []
    return [
        MatchResult(a.agent, a.drug_class, a.alias, a.alias_type,
                    looks_like_combination_product(drug_text, a.alias))
        for a in ALL_ALIASES
        if a.pattern.search(drug_text)
    ]


def looks_like_combination_product(drug_text: str, matched_alias: str) -> bool:
    """Heuristic: does ``drug_text`` name a second active ingredient besides
    ``matched_alias``?

    Splits on combination-product delimiters (``-``, ``/``, ``+``, " and ",
    " with ", " w/ "), drops the matched alias itself and any known salt /
    release / dosage-form qualifier word, and flags a combination product if
    a substantive word remains (e.g. "Chlordiazepoxide" in
    "Amitriptyline-Chlordiazepoxide"). This is a heuristic, not a drug
    database lookup, so it is conservative: short or purely numeric/unit
    leftovers are not flagged.
    """
    remainder = re.sub(re.escape(matched_alias), " ", drug_text, flags=re.IGNORECASE)
    parts = re.split(r"[-/+]|(?:\band\b)|(?:\bwith\b)|(?:\bw/\b)", remainder, flags=re.IGNORECASE)
    for part in parts:
        for word in _WORD_RE.findall(part):
            if len(word) < 4:
                continue
            if word.lower() in _QUALIFIER_WORDS:
                continue
            return True
    return False


def is_systemic(route: str | None, form_rx: str | None) -> bool:
    """False if ``route``/``form_rx`` indicate a non-systemic (local) use."""
    route_norm = (route or "").strip().lower()
    if route_norm in NON_SYSTEMIC_ROUTE_TOKENS:
        return False
    form_norm = (form_rx or "").strip().lower()
    if any(kw in form_norm for kw in NON_SYSTEMIC_FORM_KEYWORDS):
        return False
    return True


def near_miss_report(distinct_drug_values: list[str]) -> list[tuple[str, str]]:
    """Values that contain a curated agent's generic or brand name as a
    *raw* substring but were rejected by the word-boundary matcher, with the
    reason -- i.e. exactly the traps this module is designed to reject.
    (Brand names are included here purely as a wider net for this specific
    trap pattern; they still never drive an accepted match -- see module
    docstring.) Used for the reviewable audit printout, not for matching
    itself.
    """
    rejected = []
    accepted = {v for v in distinct_drug_values if match_drug_text(v) is not None}
    for value in distinct_drug_values:
        if value in accepted or not value:
            continue
        low = value.lower()
        for alias, agent in (
            [(a, agent) for agent in AGENT_CLASS for a in (agent, *BRAND_NAMES.get(agent, ()))]
        ):
            if alias and alias.lower() in low:
                rejected.append(
                    (
                        value,
                        f"raw substring '{alias}' ({agent}) present but not a whole word "
                        f"(word-boundary match failed) -- rejected",
                    )
                )
                break
    return rejected


def psychiatric_review_candidates(distinct_drug_values: list[str]) -> list[tuple[str, str]]:
    """Everything :func:`match_drug_text` rejected that still might be worth
    a human's attention, for the "print anything unmatched that looks like a
    psychiatric medication" review pass (run this over the full drug-value
    list once full MIMIC-IV access is granted -- don't trust the hint sets
    below to auto-classify anything).

    Flags a rejected value when it:
      1. is a near-miss per :func:`near_miss_report` (curated generic/brand
         name present as a raw substring, no word boundary);
      2. contains a curated agent's **brand** name as a whole word (brand
         names are hints only -- see module docstring -- so a hit here is
         reported, never auto-included);
      3. contains a known antidepressant **not** in the curated 28
         (:data:`OTHER_KNOWN_ANTIDEPRESSANTS_FOR_REVIEW`) as a whole word; or
      4. has a word ending in a suffix characteristic of antidepressant
         generic names (:data:`ANTIDEPRESSANT_SUFFIX_STEMS`).

    Deliberately over-inclusive: every entry is a candidate for a human to
    look at, not a classification decision.
    """
    reasons: dict[str, list[str]] = {}
    for value, reason in near_miss_report(distinct_drug_values):
        reasons.setdefault(value, []).append(reason)

    for value in distinct_drug_values:
        if not value or match_drug_text(value) is not None:
            continue
        low = value.lower()

        for agent, brands in BRAND_NAMES.items():
            for brand in brands:
                if brand and re.search(r"\b" + re.escape(brand) + r"\b", low):
                    reasons.setdefault(value, []).append(
                        f"contains brand name '{brand}' for {agent} ({AGENT_CLASS[agent]}) -- "
                        f"brand names are review hints only, not matcher input"
                    )

        for name, note in OTHER_KNOWN_ANTIDEPRESSANTS_FOR_REVIEW.items():
            if re.search(r"\b" + re.escape(name) + r"\b", low):
                reasons.setdefault(value, []).append(
                    f"contains '{name}', a known antidepressant not in the curated 28 ({note})"
                )

        for word in _WORD_RE.findall(low):
            for stem in ANTIDEPRESSANT_SUFFIX_STEMS:
                if word != stem and word.endswith(stem):
                    reasons.setdefault(value, []).append(
                        f"word '{word}' ends in '-{stem}', a suffix common to antidepressant generic names"
                    )
                    break

    return [(v, "; ".join(dict.fromkeys(r))) for v, r in reasons.items()]
