"""Numerical invariants; run with python -m unittest discover -s analysis -v."""
import math
import json
from pathlib import Path
import unittest
from build_data import POINTWISE_Z, classify, number, overlaps, valid_estimate


class IntervalTests(unittest.TestCase):
    def test_missing_is_not_zero(self):
        for value in [None, "NA", "nan", "inf", "", "-Inf"]:
            self.assertIsNone(number(value))
        self.assertEqual(number("0"), 0.0)

    def test_invalid_estimates(self):
        for beta, se in [(None, 1), (1, None), (1, 0), (1, -1), (math.nan, 1), (1, math.inf)]:
            self.assertFalse(valid_estimate(beta, se))
            self.assertEqual(classify(beta, se), "unavailable")

    def test_classification(self):
        self.assertEqual(classify(0.01, 0.02), "bounded-small")
        self.assertEqual(classify(0.5, 0.01), "above-threshold")
        self.assertEqual(classify(-0.5, 0.01), "above-threshold")
        self.assertEqual(classify(0.0, 0.1), "unresolved")
        self.assertEqual(classify(0.1, 0.01), "unresolved")

    def test_boundaries(self):
        self.assertEqual(classify(0, 0.1 / POINTWISE_Z), "bounded-small")
        self.assertEqual(classify(0.15, 0.025, delta=0.1, z=2), "unresolved")
        self.assertEqual(classify(0.151, 0.025, delta=0.1, z=2), "above-threshold")

    def test_sign_symmetry_and_threshold_monotonicity(self):
        for beta in [0, 0.01, 0.08, 0.2, 0.5]:
            for se in [0.005, 0.025, 0.1, 0.5]:
                for delta in [0.02, 0.1, 0.5]:
                    self.assertEqual(classify(beta, se, delta), classify(-beta, se, delta))
                if classify(beta, se, 0.1) == "bounded-small":
                    self.assertEqual(classify(beta, se, 0.5), "bounded-small")

    def test_conservative_interval_does_not_gain_precision(self):
        for beta in [0, 0.01, 0.08, 0.2, 0.5]:
            for se in [0.005, 0.025, 0.1, 0.5]:
                if classify(beta, se, z=4.7) == "bounded-small":
                    self.assertEqual(classify(beta, se), "bounded-small")

    def test_closed_interval_overlap(self):
        self.assertTrue(overlaps(1, 10, 10, 20))
        self.assertTrue(overlaps(5, 25, 10, 20))
        self.assertTrue(overlaps(12, 18, 10, 20))
        self.assertFalse(overlaps(1, 9, 10, 20))
        self.assertFalse(overlaps(21, 30, 10, 20))

    def test_invalid_controls(self):
        for delta, z in [(0, 2), (-0.1, 2), (0.1, 0), (math.nan, 2)]:
            with self.assertRaises(ValueError):
                classify(0, 1, delta, z)


class ReleasedDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((Path(__file__).resolve().parents[1] / "data.json").read_text())

    def test_coverage_regression(self):
        self.assertEqual([t["nAssayed"] for t in self.data["traits"]], [18505, 18524, 18514])
        self.assertEqual([len(t["loci"]) for t in self.data["traits"]], [50, 50, 50])
        self.assertEqual([sum(len(l["genes"]) for l in t["loci"]) for t in self.data["traits"]], [2963, 5247, 3948])

    def test_defaults_resolve_and_source_intervals(self):
        meta = self.data["meta"]
        t = next(t for t in self.data["traits"] if t["id"] == meta["defaultTrait"])
        locus = next(l for l in t["loci"] if l["id"] == meta["defaultLocus"])
        selected = next(g for g in locus["genes"] if g["id"] == meta["defaultGene"])
        self.assertEqual(selected["symbol"], "SORT1")
        self.assertEqual(locus["defaultCounts"], {"bounded-small": 4, "above-threshold": 0, "unresolved": 44, "unavailable": 22})
        self.assertEqual(locus["familyCounts"], {"bounded-small": 1, "above-threshold": 0, "unresolved": 47, "unavailable": 22})
        self.assertAlmostEqual(selected["beta"] - POINTWISE_Z * selected["se"], -0.24865317389761266)

    def test_all_intervals_and_counts(self):
        for trait in self.data["traits"]:
            for locus in trait["loci"]:
                self.assertEqual(len(locus["genes"]), len({g["id"] for g in locus["genes"]}))
                self.assertEqual(sum(locus["defaultCounts"].values()), len(locus["genes"]))
                self.assertEqual(sum(locus["familyCounts"].values()), len(locus["genes"]))
                self.assertLessEqual(locus["familyCounts"]["bounded-small"], locus["defaultCounts"]["bounded-small"])
                for gene in locus["genes"]:
                    self.assertEqual(gene["chr"], locus["chr"])
                    self.assertTrue(overlaps(gene["start"], gene["end"], locus["start"], locus["end"]))

    def test_unavailable_values_are_null(self):
        unavailable = 0
        for trait in self.data["traits"]:
            for locus in trait["loci"]:
                for gene in locus["genes"]:
                    if gene["availability"] != "available":
                        unavailable += 1
                        self.assertIsNone(gene["beta"])
                        self.assertIsNone(gene["se"])
                        self.assertIsNone(gene["burdenP"])
                    else:
                        self.assertTrue(valid_estimate(gene["beta"], gene["se"]))
                        self.assertGreaterEqual(gene["burdenP"], 0)
                        self.assertLessEqual(gene["burdenP"], 1)
        self.assertGreater(unavailable, 0)

    def test_static_javascript_matches_json(self):
        text = (Path(__file__).resolve().parents[1] / "data.js").read_text()
        self.assertTrue(text.startswith("window.SPENCE_DATA = "))
        self.assertEqual(json.loads(text.removeprefix("window.SPENCE_DATA = ").rstrip(";\n")), self.data)


if __name__ == "__main__":
    unittest.main()
