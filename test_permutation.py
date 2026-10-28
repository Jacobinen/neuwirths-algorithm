"""Regression tests for permutation.py and the Sage-free Neuwirth check.

Run with ``python -m unittest`` (a few seconds).  Set ``NEUWIRTH_SLOW_TESTS=1``
to also run the presentations with large search spaces (several minutes).

The expected values for the permutation primitives are the examples from the
SageMath documentation of the functions they replace.  The expected verdicts of
the presentations come from the original SageMath version of the code, which
was run in full on every presentation except the two with a 4.8 million
configuration search space (compared with this version on a 30 000
configuration prefix instead).
"""
import io
import os
import unittest
from contextlib import redirect_stdout

from permutation import Permutation, PermutationGroup
import neuwirthcheck


class PermutationTests(unittest.TestCase):

    def test_constructors(self):
        self.assertEqual(list(Permutation([2, 1, 3])), [2, 1, 3])
        self.assertEqual(list(Permutation('(1,2)')), [2, 1])
        self.assertEqual(list(Permutation('(1,2)(3,4,5)')), [2, 1, 4, 5, 3])
        # size is the largest element mentioned, 1-cycles are allowed
        self.assertEqual(list(Permutation('(1)(3)(2)(4)')), [1, 2, 3, 4])
        self.assertEqual(list(Permutation('(5,10,16)')),
                         [1, 2, 3, 4, 10, 6, 7, 8, 9, 16, 11, 12, 13, 14, 15, 5])
        self.assertEqual(list(Permutation('()')), [])
        with self.assertRaises(ValueError):
            Permutation('(1,2)(2,4)')
        with self.assertRaises(ValueError):
            Permutation([1, 2, 2, 4])
        with self.assertRaises(ValueError):
            Permutation([1, 2, 4, 5])

    def test_repr_and_access(self):
        p = Permutation([2, 1, 3])
        self.assertEqual(repr(p), '[2, 1, 3]')
        self.assertEqual(str(p), '[2, 1, 3]')
        self.assertEqual(len(p), 3)
        self.assertEqual(p[0], 2)
        self.assertEqual(p(1), 2)
        self.assertEqual(p[:], [2, 1, 3])

    def test_products(self):
        # sage.combinat.permutation.Permutation.right_action_product examples
        p, q = Permutation([2, 1, 3]), Permutation([3, 1, 2])
        self.assertEqual(list(p.right_action_product(q)), [1, 3, 2])
        self.assertEqual(list(q.right_action_product(p)), [3, 2, 1])
        self.assertEqual(list(p.left_action_product(q)), [3, 2, 1])
        self.assertEqual(list(q.left_action_product(p)), [1, 3, 2])
        # different sizes are padded with fixed points
        p4 = Permutation([2, 1, 3, 4])
        self.assertEqual(list(p4.right_action_product(q)), [1, 3, 2, 4])
        self.assertEqual(list(q.right_action_product(p4)), [3, 2, 1, 4])

    def test_inverse(self):
        self.assertEqual(list(Permutation([3, 8, 5, 10, 9, 4, 6, 1, 7, 2]).inverse()),
                         [8, 10, 1, 6, 3, 7, 9, 2, 5, 4])
        self.assertEqual(list(Permutation([2, 4, 1, 5, 3]).inverse()), [3, 1, 5, 2, 4])

    def test_cycles(self):
        self.assertEqual(Permutation([2, 1, 3, 4]).to_cycles(), [(1, 2), (3,), (4,)])
        self.assertEqual(Permutation([2, 1, 3, 4]).to_cycles(singletons=False), [(1, 2)])
        self.assertEqual(Permutation([2, 1, 3, 4]).to_cycles(use_min=False), [(4,), (3,), (2, 1)])
        self.assertEqual(Permutation([4, 1, 5, 2, 6, 3]).to_cycles(), [(1, 4, 2), (3, 5, 6)])
        self.assertEqual(Permutation([4, 1, 5, 2, 6, 3]).to_cycles(use_min=False), [(6, 3, 5), (4, 2, 1)])
        self.assertEqual(Permutation([3, 1, 2, 4]).cycle_type(), [3, 1])
        self.assertEqual(Permutation([2, 1, 3]).cycle_string(), '(1,2)')
        self.assertEqual(Permutation([2, 1, 3]).cycle_string(singletons=True), '(1,2)(3)')

    def test_orbits(self):
        # sage.groups.perm_gps.permgroup.PermutationGroup.orbits examples
        G = PermutationGroup([Permutation('(3,4)'), Permutation('(1,3)')])
        self.assertEqual(G.orbits(), ((1, 3, 4), (2,)))
        G = PermutationGroup([Permutation('(1,2)(3,4)'), Permutation('(1,2,3,4,10)')])
        self.assertEqual(G.orbits(), ((1, 2, 3, 4, 10), (5,), (6,), (7,), (8,), (9,)))
        # generators are sorted lexicographically (Sage canonicalisation) and the
        # points of an orbit come in GAP's breadth-first discovery order
        G = PermutationGroup([Permutation('(1,4)(2,6,5)'), Permutation('(3,4)(5,6)')])
        self.assertEqual([list(g) for g in G.gens()], [[1, 2, 4, 3, 6, 5], [4, 6, 3, 1, 2, 5]])
        self.assertEqual(G.orbits(), ((1, 4, 3), (2, 6, 5)))


class NeuwirthTests(unittest.TestCase):
    """Verdicts of the presentations listed in neuwirthcheck.main()."""

    GEN = [['a', 'A'], ['b', 'B']]
    FAST = [
        (['aaBBB', 'ABAbab'], False),
        (['aBBBa', 'BAbabA'], False),
        (['aBBBa', 'BAbaBBa'], False),
        (['aBBBa', 'AbaBBaB'], False),
        (['BBa', 'AbAbb'], True),
        (['BBa', 'Ab'], True),
        (['B', 'Ab'], True),
        (['B', 'A'], True),
    ]
    SLOW = [
        (['aBBaBBaB', 'AbaBBaB'], False),
        (['aBBaBBaB', 'BabAbbA'], False),
        (['aBBaBBaB', 'BaBaB'], True),
        (['BBaBBaBa', 'AbAbb'], True),
    ]

    def _check(self, cases):
        for rel, expected in cases:
            with self.subTest(rel=rel):
                with redirect_stdout(io.StringIO()):
                    result = neuwirthcheck.check_presentation(self.GEN, rel)
                self.assertEqual(result, expected)

    def test_fast_presentations(self):
        self._check(self.FAST)

    @unittest.skipUnless(os.environ.get('NEUWIRTH_SLOW_TESTS'), 'set NEUWIRTH_SLOW_TESTS=1 to run')
    def test_slow_presentations(self):
        self._check(self.SLOW)

    def test_spine_output(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            result = neuwirthcheck.check_presentation(self.GEN, ['BBa', 'Ab'])
        self.assertTrue(result)
        self.assertEqual(buf.getvalue(),
                         'THE FOLLOWING CONFIGURATION YIELDS A SPINE:\n\n'
                         'Reordering permutation:  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]\n'
                         '(1)(2)(3)(4)(5)(6)(7)(8)(9)(10)\n'
                         'Involution A:  [4, 5, 9, 1, 2, 8, 10, 6, 3, 7]\n'
                         '(1,4)(2,5)(3,9)(6,8)(7,10)\n')


if __name__ == '__main__':
    unittest.main()
