from grimperium_mini.metrics import absolute_deviation, relative_deviation_percent
from grimperium_mini.reactions import calculate_reactions, reaction_enthalpy


def test_reaction_enthalpy_formula():
    assert reaction_enthalpy({"A": 2 * -100.0, "B": -50.0}, {"C": -400.0}) == -150.0


def test_calculate_reactions_and_metrics():
    formation = [
        {"mol_id": "A", "method": "PM7", "hof_selected_kJmol": -100.0},
        {"mol_id": "B", "method": "PM7", "hof_selected_kJmol": -50.0},
        {"mol_id": "C", "method": "PM7", "hof_selected_kJmol": -400.0},
    ]
    reactions = [
        {
            "reaction_id": "r1",
            "reaction_type": "test",
            "compound_A": "A",
            "compound_B": "B",
            "compound_C": "C",
            "compound_D": "",
            "coef_a": "2",
            "coef_b": "1",
            "coef_c": "1",
            "coef_d": "0",
            "Hr_exp_kJmol": "-140",
        }
    ]
    out = calculate_reactions(reactions, formation)[0]
    assert out["Hr_PM7_crest_mopac_kJmol"] == -150.0
    assert out["AD_PM7_kJmol"] == 10.0


def test_ad_rd_metrics():
    assert absolute_deviation(-150.0, -140.0) == 10.0
    assert round(relative_deviation_percent(-150.0, -140.0), 6) == round(
        10 / 140 * 100, 6
    )
