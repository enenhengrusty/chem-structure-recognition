"""Unit tests for chem_cli.py."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

from PIL import Image

SCRIPTS_DIR = (
    Path(__file__).resolve().parent.parent
    / "skills"
    / "chem-structure-recognition"
    / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIR))

import chem_cli  # noqa: E402


class TestChemCli(unittest.TestCase):
    def test_sanitize_valid_and_mw_match(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            res = chem_cli.sanitize_structure(
                "CC(=O)Oc1ccccc1C(=O)O",
                note_mw=180.16,
                out_dir=tmp,
                name="aspirin",
            )
            self.assertTrue(res["valid"])
            self.assertEqual(res["formula"], "C9H8O4")
            self.assertEqual(res["mw_audit"]["status"], "MATCH")
            self.assertIn("carboxylic_acid", res["functional_motifs"])
            self.assertIn("ester", res["functional_motifs"])
            self.assertTrue(Path(res["preview_png"]).is_file())
            self.assertTrue(Path(res["preview_svg"]).is_file())

    def test_sanitize_invalid_smiles(self) -> None:
        res = chem_cli.sanitize_structure("C1(C)(C)(C)(C)C")
        self.assertFalse(res["valid"])
        self.assertIn("error", res)

    def test_audit_mw_hydrate_match(self) -> None:
        # Aspirin (180.159) + H2O (18.0153) ~= 198.174
        res = chem_cli.audit_mw_against_note(180.159, 198.17)
        self.assertEqual(res["status"], "ADDUCT_OR_SALT_MATCH")
        adducts = [a["adduct"] for a in res["adduct_matches"]]
        self.assertIn("·H2O", adducts)

    def test_audit_mw_structural_delta(self) -> None:
        # Missing one CH2 (~14.0266)
        res = chem_cli.audit_mw_against_note(180.159, 194.186)
        self.assertEqual(res["status"], "MISMATCH_REVIEW_REQUIRED")
        hyps = [h["hypothesis"] for h in res["structural_hypotheses"]]
        self.assertTrue(any("+CH2" in h for h in hyps))

    def test_render_reaction(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            res = chem_cli.render_reaction(
                reactants=["Oc1ccccc1C(=O)O", "CC(=O)OC(C)=O"],
                products=["CC(=O)Oc1ccccc1C(=O)O"],
                out_dir=tmp,
                name="aspirin_rxn",
            )
            self.assertTrue(res["valid"])
            self.assertTrue(Path(res["preview_png"]).is_file())
            self.assertTrue(Path(res["preview_svg"]).is_file())
            # Product minus reactants lacks acetic acid (C2H4O2)
            self.assertEqual(
                res["element_delta_product_minus_reactant"],
                {"C": -2, "H": -4, "O": -2},
            )

    def test_audit_table(self) -> None:
        payload = {
            "ref_name": "salicylic_acid",
            "solvent_vol_ml": 10.0,
            "actual_product_g": 1.45,
            "rows": [
                {
                    "role": "reactant",
                    "name": "salicylic_acid",
                    "smiles": "Oc1ccccc1C(=O)O",
                    "note_mw": 138.12,
                    "eq": 1.0,
                    "mass_g": 1.3812,
                },
                {
                    "role": "reactant",
                    "name": "acetic_anhydride",
                    "smiles": "CC(=O)OC(C)=O",
                    "note_mw": 102.09,
                    "eq": 1.5,
                    "mass_g": 1.531,
                },
                {
                    "role": "product",
                    "name": "aspirin",
                    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                    "eq": 1.0,
                },
            ],
        }
        res = chem_cli.audit_table(payload)
        self.assertEqual(res["ref_name"], "salicylic_acid")
        self.assertAlmostEqual(res["base_1eq_mmol_by_note_mw"], 10.0, places=3)
        self.assertAlmostEqual(res["substrate_molarity_M"], 1.0, places=2)
        prod = res["theoretical_products"][0]
        self.assertAlmostEqual(prod["theoretical_100pct_g"], 1.8016, places=2)
        self.assertGreater(prod["yield_pct"], 75.0)

    def test_stoichiometry_anomaly_warning(self) -> None:
        res = chem_cli.calculate_stoichiometry(mw=180.16, mmol=10.0, actual_g=2.1)
        self.assertGreater(res["yield_percentage"], 105.0)
        self.assertIn("anomaly_warning", res)

    def test_prep_image(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            src = tmp_path / "sample.png"
            Image.new("RGB", (400, 300), color=(245, 245, 245)).save(src)
            out_dir = tmp_path / "prep"
            res = chem_cli.prep_image(str(src), str(out_dir), enhance=True)
            self.assertEqual(res["width"], 400)
            self.assertEqual(res["height"], 300)
            self.assertTrue(Path(res["upright"]).is_file())
            self.assertTrue(Path(res["enhanced"]).is_file())
            self.assertEqual(len(res["crops"]), 6)


if __name__ == "__main__":
    unittest.main()
