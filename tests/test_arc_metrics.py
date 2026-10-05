"""Edge-case coverage for ARC-style evaluation metrics."""

from math import inf

import numpy as np

from vcm.evaluation.arc_metrics import ARCEvaluator


def test_des_uses_available_gene_count_when_less_than_top_n():
    evaluator = ARCEvaluator()
    assert evaluator.calculate_des({"A": 2.0}, {"A": 3.0}, top_n=100) == 1.0


def test_des_partial_overlap_uses_effective_denominator():
    evaluator = ARCEvaluator()
    score = evaluator.calculate_des({"A": 3.0, "B": 2.0}, {"A": 4.0, "C": 1.0}, top_n=100)
    assert score == 0.5


def test_des_empty_inputs_return_zero():
    evaluator = ARCEvaluator()
    assert evaluator.calculate_des({}, {}) == 0.0
    assert evaluator.calculate_des({"A": 1.0}, {"A": 1.0}, top_n=0) == 0.0


def test_des_ties_are_deterministic():
    evaluator = ARCEvaluator()
    assert evaluator.calculate_des({"B": 1.0, "A": 1.0}, {"A": 1.0}, top_n=1) == 1.0


def test_pds_identifies_correct_perturbation():
    evaluator = ARCEvaluator()
    effects = {"sample": {"drug_a": 0.2, "drug_b": 0.9}}
    assert evaluator.calculate_pds(effects, {"sample": "drug_b"}) == 1.0


def test_pds_skips_missing_sample_ids():
    evaluator = ARCEvaluator()
    assert evaluator.calculate_pds({"other": {"drug_a": 1.0}}, {"sample": "drug_a"}) == 0.0


def test_pds_empty_inputs_return_zero():
    evaluator = ARCEvaluator()
    assert evaluator.calculate_pds({}, {}) == 0.0


def test_pds_ties_use_lexicographic_perturbation_id():
    evaluator = ARCEvaluator()
    effects = {"sample": {"zeta": 1.0, "alpha": 1.0}}
    assert evaluator.calculate_pds(effects, {"sample": "alpha"}) == 1.0


def test_mae_aligns_only_common_genes():
    evaluator = ARCEvaluator()
    assert evaluator.calculate_mae({"A": 3.0, "missing": 99.0}, {"A": 1.0}) == 2.0


def test_mae_weights_viral_genes_more_heavily():
    evaluator = ARCEvaluator()
    score = evaluator.calculate_mae({"viral_LT": 3.0, "house": 3.0}, {"viral_LT": 1.0, "house": 1.0})
    assert score == 2.0


def test_mae_custom_weights_override_defaults():
    evaluator = ARCEvaluator()
    score = evaluator.calculate_mae({"A": 3.0}, {"A": 1.0}, {"housekeeping": 2.0})
    assert score == 2.0


def test_mae_supports_explicit_gene_weight_override():
    evaluator = ARCEvaluator()
    score = evaluator.calculate_mae(
        {"A": 3.0, "B": 3.0},
        {"A": 1.0, "B": 1.0},
        {"A": 3.0, "housekeeping": 1.0},
    )
    assert score == 2.0


def test_mae_non_finite_values_return_infinity():
    evaluator = ARCEvaluator()
    assert evaluator.calculate_mae({"A": np.nan}, {"A": 1.0}) == inf
    assert evaluator.calculate_mae({"A": np.inf}, {"A": 1.0}) == inf


def test_mae_no_common_genes_returns_infinity():
    evaluator = ARCEvaluator()
    assert evaluator.calculate_mae({"A": 1.0}, {"B": 1.0}) == inf


def test_mae_none_inputs_return_infinity():
    evaluator = ARCEvaluator()
    assert evaluator.calculate_mae(None, None) == inf


def test_cell_cycle_scores_handle_missing_markers():
    evaluator = ARCEvaluator()
    s_score, g2m_score, phase = evaluator.calculate_cell_cycle_scores({}, ["MCM2"], ["TOP2A"])
    assert (s_score, g2m_score, phase) == (0.0, 0.0, "G1")


def test_evaluate_full_returns_all_metrics():
    evaluator = ARCEvaluator()
    result = evaluator.evaluate_full(
        {"A": 2.0}, {"A": 2.0},
        {"sample": {"drug": 1.0}}, {"sample": "drug"},
        {"A": 2.0}, {"A": 2.0},
    )
    assert (result.des, result.pds, result.mae) == (1.0, 1.0, 0.0)
