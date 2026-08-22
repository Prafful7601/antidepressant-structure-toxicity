"""Unit tests for the antidepressant free-text drug matcher.

Includes deliberate traps: strings where a curated agent's name is a raw
substring of something that is NOT that agent, to prove the matcher does
not rely on naive substring matching.
"""

from __future__ import annotations

import itertools

import pytest

from mimic4_pipeline.drug_mapping import (
    AGENT_CLASS,
    BRAND_NAMES,
    is_systemic,
    looks_like_combination_product,
    match_drug_text,
    near_miss_report,
    psychiatric_review_candidates,
)


def test_curated_agent_list_has_28_agents_correct_classes():
    assert len(AGENT_CLASS) == 28
    counts = {}
    for cls in AGENT_CLASS.values():
        counts[cls] = counts.get(cls, 0) + 1
    assert counts == {"Tricyclic": 9, "SSRI": 6, "SNRI": 5, "RIMA": 2, "Atypical": 6}


def test_curated_agent_names_do_not_collide_undetected():
    """Documents the 4 real substring traps found among the 28 curated
    agents themselves, and proves the matcher resolves all of them
    correctly via word-boundary matching."""
    traps = []
    for a, b in itertools.permutations(AGENT_CLASS, 2):
        if a in b:
            traps.append((a, b))
    # These four are the actual, known collisions -- if this assertion ever
    # fails because the set changed, the new pair needs its own explicit
    # match_drug_text assertion added below.
    assert set(traps) == {
        ("imipramine", "trimipramine"),
        ("citalopram", "escitalopram"),
        ("venlafaxine", "desvenlafaxine"),
        ("milnacipran", "levomilnacipran"),
    }


@pytest.mark.parametrize(
    "text,expected_agent",
    [
        ("Trimipramine", "trimipramine"),
        ("TRIMIPRAMINE 25mg Tab", "trimipramine"),
        ("Escitalopram Oxalate", "escitalopram"),
        ("Desvenlafaxine Succinate ER", "desvenlafaxine"),
        ("Levomilnacipran ER", "levomilnacipran"),
    ],
)
def test_longer_agent_name_not_misclassified_as_embedded_shorter_agent(text, expected_agent):
    """The core trap: does 'citalopram' embedded in 'escitalopram' get
    misclassified as citalopram instead of escitalopram? (And likewise for
    the other three embedded pairs.)"""
    result = match_drug_text(text)
    assert result is not None
    assert result.agent == expected_agent


@pytest.mark.parametrize(
    "text",
    [
        "Imipramine",
        "IMIPRAMINE HCl 25mg",
    ],
)
def test_shorter_agent_still_matches_on_its_own(text):
    result = match_drug_text(text)
    assert result is not None
    assert result.agent == "imipramine"


def test_case_insensitive_matching():
    for text in ["amitriptyline", "AMITRIPTYLINE", "AmiTriptyline", "BuPROPion"]:
        result = match_drug_text(text)
        assert result is not None


def test_unrelated_drug_containing_agent_name_as_raw_substring_is_rejected():
    """Deliberate trap: a fictional compound whose name contains
    'amitriptyline' as a raw substring, immediately followed by another
    letter (no word boundary). A naive `"amitriptyline" in text.lower()`
    check would misclassify this as amitriptyline exposure."""
    result = match_drug_text("Pseudamitriptylinex Complex")
    assert result is None
    # sanity: the real agent still matches on its own
    assert match_drug_text("Amitriptyline").agent == "amitriptyline"


def test_unrelated_word_containing_agent_stem_is_rejected():
    """Deliberate trap: 'trazodone' embedded with a letter directly on
    each side, so it is a raw substring but not a whole word."""
    assert match_drug_text("Xtrazodonex Complex") is None
    # sanity: the real agent still matches on its own
    assert match_drug_text("Trazodone").agent == "trazodone"


def test_no_match_for_unrelated_drugs():
    for text in ["Metoprolol Tartrate", "Sodium Chloride 0.9%", "Insulin Human Regular", "Tricor"]:
        assert match_drug_text(text) is None


def test_salt_form_suffix_still_matches_and_is_not_flagged_as_combination():
    result = match_drug_text("Amitriptyline HCl")
    assert result is not None and result.agent == "amitriptyline"
    assert not looks_like_combination_product("Amitriptyline HCl", "amitriptyline")


def test_release_modifier_suffix_still_matches_and_is_not_combination():
    result = match_drug_text("BuPROPion XL (Once Daily)")
    assert result is not None and result.agent == "bupropion"
    assert not result.is_combination_product


def test_combination_product_is_flagged():
    result = match_drug_text("Amitriptyline-Chlordiazepoxide")
    assert result is not None and result.agent == "amitriptyline"
    assert result.is_combination_product is True


def test_olanzapine_fluoxetine_combination_is_flagged():
    result = match_drug_text("OLANZapine-Fluoxetine (Symbyax)")
    assert result is not None and result.agent == "fluoxetine"
    assert result.is_combination_product is True


def test_brand_names_do_not_drive_automatic_matching():
    """Brand names are review hints only (per explicit instruction): a pure
    brand-name string must NOT be auto-classified as exposure."""
    assert match_drug_text("Zoloft 50mg Tablet") is None
    assert match_drug_text("Paxil CR") is None


def test_brand_name_hit_surfaces_as_review_candidate_not_a_match():
    values = ["Zoloft 50mg Tablet"]
    assert match_drug_text(values[0]) is None
    report = dict(psychiatric_review_candidates(values))
    assert "Zoloft 50mg Tablet" in report
    assert "sertraline" in report["Zoloft 50mg Tablet"]


def test_brand_name_substring_trap_does_not_spuriously_match():
    """'Paxil' brand for paroxetine should not spuriously match as exposure
    on unrelated text with no word boundary -- but it IS a legitimate raw
    substring, so it should still surface as a review candidate."""
    assert match_drug_text("Paxillin Injection") is None
    report = dict(psychiatric_review_candidates(["Paxillin Injection"]))
    assert "Paxillin Injection" in report
    assert "paxil" in report["Paxillin Injection"]


def test_other_known_antidepressant_not_curated_is_flagged_for_review_only():
    """A real antidepressant outside the curated 28 (e.g. agomelatine) must
    never be auto-classified, but should be flagged for a human to see."""
    assert match_drug_text("Agomelatine 25mg") is None
    report = dict(psychiatric_review_candidates(["Agomelatine 25mg"]))
    assert "Agomelatine 25mg" in report
    assert "agomelatine" in report["Agomelatine 25mg"]


def test_suffix_stem_heuristic_flags_unmatched_novel_name():
    """A fictional/misspelled generic-sounding name should be flagged by
    the suffix heuristic even though it matches nothing else."""
    assert match_drug_text("Zoravoxetine") is None
    report = dict(psychiatric_review_candidates(["Zoravoxetine"]))
    assert "Zoravoxetine" in report
    assert "oxetine" in report["Zoravoxetine"]


def test_unrelated_drug_is_not_flagged_by_review_pass():
    report = dict(psychiatric_review_candidates(["Metoprolol Tartrate", "Insulin Human Regular"]))
    assert report == {}


def test_non_systemic_route_excluded():
    assert is_systemic("PO", None) is True
    assert is_systemic("PO/NG", None) is True
    assert is_systemic("IV", None) is True
    assert is_systemic("TOPICAL", None) is False
    assert is_systemic("OU", None) is False  # both eyes, ophthalmic
    assert is_systemic(None, "Doxepin Cream 5%") is False
    assert is_systemic(None, "Trazodone Patch") is False


def test_near_miss_report_flags_rejected_traps_only():
    values = [
        "Amitriptyline",  # accepted -> not a near-miss
        "Pseudamitriptylinex Complex",  # rejected trap -> should be reported
        "Metoprolol Tartrate",  # no agent substring at all -> not reported
    ]
    report = dict(near_miss_report(values))
    assert "Pseudamitriptylinex Complex" in report
    assert "Amitriptyline" not in report
    assert "Metoprolol Tartrate" not in report


def test_every_agent_has_at_least_a_generic_alias_pattern():
    for agent in AGENT_CLASS:
        result = match_drug_text(agent.capitalize())
        assert result is not None and result.agent == agent


def test_clorgyline_never_marketed_has_no_brand_names():
    # Documented in-code: clorgyline was never marketed. Confirm the data
    # structure reflects that rather than a placeholder brand name sneaking in.
    assert BRAND_NAMES["clorgyline"] == ()
