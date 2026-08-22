"""Unit tests for the dose-based multi-class reclassification rule (task 3,
primary rule as of 2026-08-23) -- built directly against synthetic
prescription rows rather than the database, so these run without a build.
"""

from __future__ import annotations

import pandas as pd

from mimic4_pipeline.cohort import _exposure_per_admission


def _rx(hadm_id, drug_class, agent, dose_val_rx="100", doses_per_24_hrs=1.0, is_systemic=True):
    return {
        "hadm_id": hadm_id,
        "drug_class": drug_class,
        "agent": agent,
        "dose_val_rx": dose_val_rx,
        "doses_per_24_hrs": doses_per_24_hrs,
        "is_systemic": is_systemic,
    }


def _by_hadm(df: pd.DataFrame, hadm_id: int) -> dict:
    row = df[df["hadm_id"] == hadm_id]
    assert len(row) == 1
    return row.iloc[0].to_dict()


def test_single_class_admission_is_unaffected():
    rx = pd.DataFrame([_rx(1, "SSRI", "sertraline", "100")])
    result = _by_hadm(_exposure_per_admission(rx), 1)
    assert result["exposure_class"] == "SSRI"
    assert result["exposure_assignment_method"] == "single_class"
    assert result["reclassified_away_classes"] == ""


def test_low_dose_amitriptyline_plus_ssri_reclassifies_to_ssri():
    rx = pd.DataFrame(
        [
            _rx(1, "Tricyclic", "amitriptyline", "50"),  # below 75mg threshold -> low_dose
            _rx(1, "SSRI", "sertraline", "100"),
        ]
    )
    result = _by_hadm(_exposure_per_admission(rx), 1)
    assert result["exposure_class"] == "SSRI"
    assert result["exposure_assignment_method"] == "dose_reclassified"
    assert result["reclassified_away_classes"] == "Tricyclic"
    assert result["reclassified_away_agents"] == "amitriptyline"
    # dose_stratum reflects the FINAL assigned class (SSRI has no dose
    # ambiguity), not the dropped Tricyclic class's low_dose status.
    assert result["dose_stratum"] == "not_applicable"


def test_low_dose_trazodone_plus_ssri_reclassifies_to_ssri():
    rx = pd.DataFrame(
        [
            _rx(1, "Atypical", "trazodone", "25"),  # below 150mg threshold -> low_dose
            _rx(1, "SSRI", "sertraline", "100"),
        ]
    )
    result = _by_hadm(_exposure_per_admission(rx), 1)
    assert result["exposure_class"] == "SSRI"
    assert result["exposure_assignment_method"] == "dose_reclassified"
    assert result["reclassified_away_classes"] == "Atypical"


def test_standard_dose_amitriptyline_plus_ssri_remains_genuine_multiclass():
    rx = pd.DataFrame(
        [
            _rx(1, "Tricyclic", "amitriptyline", "150"),  # >=75mg -> standard_dose: real TCA use
            _rx(1, "SSRI", "sertraline", "100"),
        ]
    )
    result = _by_hadm(_exposure_per_admission(rx), 1)
    assert result["exposure_class"] is None
    assert result["exposure_assignment_method"] == "excluded_multiclass"


def test_unknown_dose_amitriptyline_plus_ssri_remains_excluded_not_reclassified():
    # Unparseable dose is not evidence of low dose -- conservative: don't
    # reclassify away a class we can't actually confirm is low-dose.
    rx = pd.DataFrame(
        [
            _rx(1, "Tricyclic", "amitriptyline", "QID"),
            _rx(1, "SSRI", "sertraline", "100"),
        ]
    )
    result = _by_hadm(_exposure_per_admission(rx), 1)
    assert result["exposure_class"] is None
    assert result["exposure_assignment_method"] == "excluded_multiclass"


def test_two_non_dose_ambiguous_classes_remains_excluded():
    rx = pd.DataFrame(
        [
            _rx(1, "SSRI", "sertraline", "100"),
            _rx(1, "SNRI", "duloxetine", "60"),
        ]
    )
    result = _by_hadm(_exposure_per_admission(rx), 1)
    assert result["exposure_class"] is None
    assert result["exposure_assignment_method"] == "excluded_multiclass"


def test_bupropion_presence_blocks_reclassification_of_its_own_atypical_class():
    # Atypical here is backed by BOTH low-dose trazodone AND bupropion --
    # bupropion has no dose ambiguity, so its presence is real
    # antidepressant-dose evidence and the whole Atypical class must not be
    # dropped just because trazodone alone was low-dose.
    rx = pd.DataFrame(
        [
            _rx(1, "Atypical", "trazodone", "25"),
            _rx(1, "Atypical", "bupropion", "150"),
            _rx(1, "SSRI", "sertraline", "100"),
        ]
    )
    result = _by_hadm(_exposure_per_admission(rx), 1)
    assert result["exposure_class"] is None
    assert result["exposure_assignment_method"] == "excluded_multiclass"


def test_all_classes_low_dose_ambiguous_only_has_nothing_to_assign():
    rx = pd.DataFrame(
        [
            _rx(1, "Tricyclic", "amitriptyline", "50"),
            _rx(1, "Atypical", "trazodone", "25"),
        ]
    )
    result = _by_hadm(_exposure_per_admission(rx), 1)
    assert result["exposure_class"] is None
    assert result["exposure_assignment_method"] == "excluded_ambiguous_low_dose_only"
    assert set(result["reclassified_away_classes"].split("|")) == {"Tricyclic", "Atypical"}


def test_n_classes_raw_reflects_pre_reclassification_count_for_sensitivity_analysis():
    rx = pd.DataFrame(
        [
            _rx(1, "Tricyclic", "amitriptyline", "50"),
            _rx(1, "SSRI", "sertraline", "100"),
        ]
    )
    result = _by_hadm(_exposure_per_admission(rx), 1)
    # exposure_class is reclassified to SSRI, but n_classes_raw still shows
    # 2 -- the original pure-exclusion sensitivity cohort is reproduced by
    # filtering n_classes_raw == 1, independent of this reclassification.
    assert result["n_classes_raw"] == 2
    assert result["exposure_class"] == "SSRI"


def test_non_systemic_rows_are_excluded_from_exposure_before_reclassification():
    rx = pd.DataFrame(
        [
            _rx(1, "Tricyclic", "amitriptyline", "150", is_systemic=False),  # topical/non-systemic, ignored
            _rx(1, "SSRI", "sertraline", "100"),
        ]
    )
    result = _by_hadm(_exposure_per_admission(rx), 1)
    assert result["exposure_class"] == "SSRI"
    assert result["exposure_assignment_method"] == "single_class"
