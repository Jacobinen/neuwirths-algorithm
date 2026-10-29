"""Tests for the command line (check.py) and the Heegaard diagrams (heegaard.py).

Run with ``uv run python -m unittest`` (about half a minute: the diagrams are laid out
with a force-directed relaxation).
"""
import math
import os
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

import check
import heegaard
import pretzel
import surface

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
        self.assertIn(r"\mathcal{P} = \langle a, b \mid b^{-2}a, a^{-1}b \rangle", html)   # the heading, for KaTeX
        tikz = d.to_tikz()
        self.assertIn("\\end{tikzpicture}", tikz)
        self.assertEqual(tikz.count("\\draw[disk]"), 4)
        self.assertEqual(tikz.count("\\draw[beta"), 5)
        self.assertEqual(d.presentation_text(latex=True), r"\langle a, b \mid b^{-2}a, a^{-1}b \rangle")


class SurfaceTests(unittest.TestCase):

    def _drawing(self, rels):
        gens, rels = check.parse_presentation(rels)
        d = heegaard.diagram_for_presentation(gens, rels, relax_iterations=120)
        return d, surface.SurfaceDrawing(d)

    def test_strands_match_the_planar_arcs(self):
        for rels in (["BBa", "Ab"], ["BBa", "AbAbb"], ["B", "A"], ["a", "b", "c"]):
            with self.subTest(rels=rels):
                d, sd = self._drawing(rels)
                self.assertEqual(len(sd.handles), d.g)
                self.assertEqual(sum(len(h["strands"]) for h in sd.handles), sum(len(w) for w in rels))
                self.assertEqual(len(sd.joins), d.L - 1)
                for i, h in enumerate(sd.handles):
                    tube = h["tube"]
                    self.assertEqual(sorted(s["sheet"] for s in h["strands"]), list(range(1, d.prev[i] + 1)))
                    # a strand starts at the point of the foot circle where the planar arc ends,
                    # at both feet, and its twist is a lift of the difference of the two angles
                    for s in h["strands"]:
                        e = [e for e in d.vertex_darts[2 * i] if d.dart_info[e][1] == s["sheet"]][0]
                        for foot, dart, pos in ((2 * i, e, 0.0), (2 * i + 1, d.B(e), math.pi)):
                            phi = s["lift_p"] if pos == 0.0 else s["lift_p"] + s["twist"]
                            x, y = tube.point(pos, phi)
                            cx, cy = sd.feet[foot]
                            ang = sd.dart_angle[dart]
                            self.assertAlmostEqual(x, cx + sd.r * math.cos(ang), places=9)
                            self.assertAlmostEqual(y, cy + sd.r * math.sin(ang), places=9)
                        self.assertTrue(s["runs"])          # (a strand may run entirely behind the tube)
                    # strands keep their cyclic order along the tube: lifted angles increase
                    lifts = [s["lift_p"] for s in h["strands"]]
                    self.assertEqual(lifts, sorted(lifts))
                    ends = [s["lift_p"] + s["twist"] for s in h["strands"]]
                    self.assertEqual(ends, sorted(ends))

    def test_outputs(self):
        d, sd = self._drawing(["BBa", "Ab"])
        root = ET.fromstring(sd.to_svg())
        self.assertTrue(root.tag.endswith("svg"))
        texts = [el.text for el in root.iter() if el.tag.endswith("text") and el.text]
        self.assertEqual(sum(1 for t in texts if t.startswith("\u03b1")), 2)       # alpha_a, alpha_b
        tikz = sd.to_tikz()
        self.assertIn("\\end{tikzpicture}", tikz)
        self.assertEqual(tikz.count("alpha_{"), 2)
        page = d.to_html()
        self.assertIn("cut open along the", page)
        self.assertEqual(page.count("<svg"), 2)


class PretzelTests(unittest.TestCase):

    def test_hoop_geometry(self):
        P = (math.cos(0.4), math.sin(0.4), 0.0)
        Q = (math.cos(2.1), math.sin(2.1), 0.0)
        h = pretzel.Hoop(P, Q, 0.2)
        self.assertLess(math.dist(h.centre(0.0), P), 1e-9)
        self.assertLess(math.dist(h.centre(h.length), Q), 1e-9)
        self.assertAlmostEqual(h.radius(0.0), 0.2)
        self.assertAlmostEqual(h.radius(h.length / 2), h.r)
        for k in range(1, 10):
            s = h.length * k / 10
            for phi in (0.3, 2.0, 4.5):
                p = h.point(s, phi)
                dist, s2 = h.distance(p)
                self.assertAlmostEqual(dist, h.radius(s), places=9)
                self.assertAlmostEqual(s2, s, places=9)
                direction = pretzel._add(p, h.centre(s), -1.0)
                self.assertAlmostEqual((h.angle_of(s, direction) - phi + math.pi) % (2 * math.pi) - math.pi, 0.0, places=9)
        # the hoop lies in the plane of the equator: all of it has z = 0
        self.assertLess(max(abs(h.centre(h.length * k / 20)[2]) for k in range(21)), 1e-9)

    def test_pictures(self):
        for rels, flat in ((["BBa", "Ab"], True), (["B", "Ab"], True), (["aab", "bbA"], None)):
            with self.subTest(rels=rels):
                gens, rels = check.parse_presentation(rels)
                d = heegaard.diagram_for_presentation(gens, rels, relax_iterations=120)
                picture = pretzel.Pretzel(d)
                if flat is not None:
                    self.assertEqual(picture.flat, flat)
                self.assertEqual(len(picture.hoops), d.g)
                self.assertEqual(len(picture.alpha), d.g)
                self.assertEqual(sum(len(st) for st in picture.strands), sum(len(w) for w in rels))
                # the beta-curves are closed chains of arcs and strands, one strand per letter,
                # lying on the surface
                for r, pieces in enumerate(picture.beta):
                    self.assertEqual(len(pieces), 2 * len(rels[r]))
                    for a, b in zip(pieces, pieces[1:] + pieces[:1]):
                        self.assertLess(math.dist(a["points"][-1], b["points"][0]), 1e-9)
                    for piece in pieces:
                        for p in piece["points"]:
                            self.assertLess(abs(picture.field(p)), 0.1)
                # each strand twists by at most the allowed amount around its hoop
                for st in picture.strands:
                    for s in st.values():
                        self.assertLessEqual(abs(s["twist"]), pretzel.MAX_TWIST + 0.4)
                # the surface has the right genus
                self.assertEqual(picture.genus(h=0.16), d.g)

    def test_page(self):
        gens, rels = check.parse_presentation(["BBa", "Ab"])
        d = heegaard.diagram_for_presentation(gens, rels, relax_iterations=120)
        page = d.to_html()
        self.assertIn("pretzel-viewer", page)
        self.assertIn("three.module.js", page)
        self.assertIn("cut open along the", page)
        self.assertEqual(page.count("<svg"), 2)          # the fallback drawing and the planar diagram
        # a diagram with several components keeps the drawing of a ball with handles
        gens, rels = check.parse_presentation(["a", "b", "c"])
        d = heegaard.diagram_for_presentation(gens, rels, relax_iterations=120)
        self.assertIsNone(d.pretzel())
        page = d.to_html()
        self.assertNotIn("pretzel-viewer", page)
        self.assertEqual(page.count("<svg"), 2)


if __name__ == "__main__":
    unittest.main()
