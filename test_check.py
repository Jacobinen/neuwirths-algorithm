"""Tests for the command line (check.py) and the Heegaard diagrams (heegaard.py).

Run with ``uv run python -m unittest`` (about half a minute: the diagrams are laid out
with a force-directed relaxation).
"""
import os
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

import check
import heegaard

HERE = os.path.dirname(os.path.abspath(__file__))


def run_cli(*args, cwd=None):
    return subprocess.run([sys.executable, os.path.join(HERE, "check.py")] + list(args),
                          capture_output=True, text=True, cwd=cwd or HERE)


class ParseTests(unittest.TestCase):

    def test_generators_are_inferred(self):
        gens, rels = check.parse_presentation(["BBa", "Ab"])
        self.assertEqual(gens, [["a", "A"], ["b", "B"]])
        self.assertEqual(rels, ["BBa", "Ab"])
        gens, _ = check.parse_presentation(["c", "a", "b"])
        self.assertEqual([g[0] for g in gens], ["a", "b", "c"])

    def test_explicit_generators(self):
        gens, _ = check.parse_presentation(["ab", "ba"], generators="ab")
        self.assertEqual([g[0] for g in gens], ["a", "b"])
        with self.assertRaises(ValueError):
            check.parse_presentation(["ac", "ca"], generators="ab")      # c not declared
        with self.assertRaises(ValueError):
            check.parse_presentation(["aa", "aA"], generators="ab")      # b never occurs

    def test_rejects_bad_input(self):
        with self.assertRaises(ValueError):
            check.parse_presentation(["BBa"])                  # not balanced
        with self.assertRaises(ValueError):
            check.parse_presentation(["a1", "b"])              # not a letter
        with self.assertRaises(ValueError):
            check.parse_presentation([])

    def test_presentation_text(self):
        gens, rels = check.parse_presentation(["BBa", "Ab"])
        self.assertEqual(check.presentation_text(gens, rels), "<a, b | b^-2 a, a^-1 b>")


class CliTests(unittest.TestCase):

    def test_spine(self):
        result = run_cli("BBa", "Ab")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("<a, b | b^-2 a, a^-1 b>", result.stdout)
        self.assertIn("IS the spine", result.stdout)
        self.assertIn("A = (1,4)(2,5)(3,9)(6,8)(7,10)", result.stdout)

    def test_not_a_spine(self):
        result = run_cli("aaBBB", "ABAbab", "-q")
        self.assertEqual(result.returncode, 1)
        self.assertIn("NOT the spine", result.stdout)

    def test_bad_input(self):
        result = run_cli("BBa")
        self.assertEqual(result.returncode, 2)
        self.assertIn("balanced", result.stderr)

    def test_writes_diagram_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = run_cli("BBa", "Ab", "--heegaard", "--no-open",
                             "--svg", "d.svg", "--tikz", "d.tex", cwd=tmp)
            self.assertEqual(result.returncode, 0, result.stderr)
            for name in ("heegaard-BBa-Ab.html", "d.svg", "d.tex"):
                self.assertTrue(os.path.exists(os.path.join(tmp, name)), name)
            ET.parse(os.path.join(tmp, "d.svg"))       # well-formed XML
            with open(os.path.join(tmp, "d.tex")) as fh:
                self.assertIn("\\begin{tikzpicture}", fh.read())

    def test_no_diagram_without_spine(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = run_cli("aaBBB", "ABAbab", "--heegaard", "--no-open", cwd=tmp)
            self.assertEqual(result.returncode, 1)
            self.assertIn("No Heegaard diagram", result.stdout)
            self.assertEqual(os.listdir(tmp), [])


class DiagramTests(unittest.TestCase):

    CASES = [
        # relators, L, number of regions
        (["BBa", "Ab"], 1, 3),
        (["B", "Ab"], 1, 1),
        (["B", "A"], 2, 2),
        (["BBa", "AbAbb"], 1, 6),
        (["a", "b", "c"], 3, 3),
    ]

    def _diagram(self, rels):
        gens, rels = check.parse_presentation(rels)
        return heegaard.diagram_for_presentation(gens, rels, relax_iterations=120)

    def test_combinatorics_and_planarity(self):
        for rels, L, regions in self.CASES:
            with self.subTest(rels=rels):
                d = self._diagram(rels)
                self.assertIsNotNone(d)
                self.assertEqual(d.L, L)
                self.assertEqual(len(d.faces), regions)
                self.assertEqual(len(d.arcs), sum(len(w) for w in rels))
                # every sphere: V - E + F = 2
                self.assertEqual(d.euler_characteristics(), [2] * L)
                # the drawing is a planar embedding with the prescribed rotation system
                self.assertEqual(d.crossings(), 0)
                self.assertTrue(d.rotations_ok())
                # the arcs follow the relators: exit dart of a letter, entry dart of the next
                for arc in d.arcs:
                    self.assertEqual(d.Amap[arc["exit"]], arc["entry"])
                # disks are visible
                for layout in d.layouts:
                    for r in layout["radii"].values():
                        self.assertGreater(r, 0.03)
                if L > 1:
                    self.assertEqual(len(d.tube_faces), L)

    def test_not_a_spine_gives_none(self):
        gens, rels = check.parse_presentation(["aaBBB", "ABAbab"])
        self.assertIsNone(heegaard.diagram_for_presentation(gens, rels))

    def test_bundles_partition_the_arcs(self):
        d = self._diagram(["BBa", "AbAbb"])
        comp = d.components[0]
        bundles = d._bundles(comp)
        darts = sorted(e for b in bundles for e in b)
        self.assertEqual(len(darts), len(d.arcs))
        self.assertEqual(sorted({d.edge_of[e] for e in darts}), list(range(len(d.arcs))))
        # the a-b and abar-bbar arcs form 3-strand bundles; the two b-bbar arcs enclose
        # the flowers at a and abar, so they are not parallel
        self.assertEqual(sorted(len(b) for b in bundles), [1, 1, 3, 3])

    def test_outputs(self):
        d = self._diagram(["BBa", "Ab"])
        svg = d.to_svg()
        root = ET.fromstring(svg)
        self.assertTrue(root.tag.endswith("svg"))
        self.assertEqual(len([el for el in root.iter() if el.tag.endswith("circle") and
                              el.get("stroke") == "#222"]), 4)          # the four disks
        html = d.to_html()
        self.assertIn("&lt;a, b | b^-2 a, a^-1 b&gt;", html)
        tikz = d.to_tikz()
        self.assertIn("\\end{tikzpicture}", tikz)
        self.assertEqual(tikz.count("\\draw[disk]"), 4)
        self.assertEqual(tikz.count("\\draw[beta"), 5)
        self.assertEqual(d.presentation_text(latex=True), r"\langle a, b \mid b^{-2}a, a^{-1}b \rangle")


if __name__ == "__main__":
    unittest.main()
