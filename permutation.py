"""Minimal, dependency-free stand-ins for the SageMath objects used by Neuwirth's algorithm.

Only the handful of methods that ``neuwirthcheck.py`` calls are provided, and each
one reproduces the behaviour of the corresponding SageMath / GAP routine exactly
(same conventions, same orderings), so that the algorithm file could be ported
unchanged apart from its imports:

``Permutation``
    replaces ``sage.combinat.permutation.Permutation``.  Construct it from a
    cycle string ``'(1,2)(3,4,5)'`` (size = largest element mentioned, as in Sage)
    or from a one-line value table ``[2, 1, 5, 3, 4]``.

``PermutationGroup``
    replaces ``sage.groups.perm_gps.permgroup.PermutationGroup``; only ``orbits()``
    is provided.  It reproduces (a) Sage's domain ``{1, ..., n}`` with ``n`` the
    size of the largest generator, (b) Sage's canonicalisation of the generator
    list (duplicates removed, generators sorted lexicographically by one-line
    notation) and (c) the order in which GAP's ``OrbitsPerms`` lists orbits and
    the points inside each orbit.  (c) matters: ``construct_F`` in the algorithm
    takes the *last* point of each orbit.
"""

import itertools


class Permutation:
    """A permutation of ``{1, ..., n}`` stored in one-line notation.

    ``p[i]`` (0-based, like Sage) is the image of ``i + 1``; ``p(i)`` (1-based)
    is the image of ``i``.
    """

    __slots__ = ("_list",)

    def __init__(self, data):
        if isinstance(data, Permutation):
            self._list = list(data._list)
        elif isinstance(data, str):
            self._list = self._from_cycle_string(data)
        else:
            lst = [int(x) for x in data]
            self._check_one_line(lst)
            self._list = lst

    # -- construction helpers -------------------------------------------------

    @staticmethod
    def _from_cycle_string(s):
        # Same parsing as sage.combinat.permutation.Permutation.__classcall_private__:
        # split on ")(", strip the outer parentheses, size = largest entry.
        if s == "()" or s == "":
            return []
        cycles = s.split(")(")
        cycles[0] = cycles[0][1:]
        cycles[-1] = cycles[-1][:-1]
        cycle_list = [[int(k) for k in c.split(",")] for c in cycles]
        n = max(max(c) for c in cycle_list)
        return Permutation._from_cycles(n, cycle_list)

    @staticmethod
    def _from_cycles(n, cycles):
        # Same checks as sage.combinat.permutation.from_cycles.
        p = n * [None]
        for cycle in cycles:
            cycle_length = len(cycle)
            for i in range(cycle_length):
                k = cycle[i]
                pk = cycle[(i + 1) % cycle_length]
                if k < 1 or pk < 1:
                    raise ValueError("all elements should be strictly positive "
                                     "integers, but I found %s" % min(k, pk))
                if k > n or pk > n:
                    raise ValueError("you claimed that this is a permutation on "
                                     "1...%s, but it contains %s" % (n, max(k, pk)))
                if p[k - 1] is not None:
                    raise ValueError("the element %s appears more than once "
                                     "in the input" % k)
                p[k - 1] = pk
        for i in range(n):
            if p[i] is None:
                p[i] = i + 1
        return p

    @staticmethod
    def _check_one_line(lst):
        # Same checks as sage.combinat.permutation.Permutation.__init__.
        if not lst:
            return
        for i in lst:
            if i < 1:
                raise ValueError("the elements must be strictly positive integers")
        srt = sorted(lst)
        if srt[-1] != len(srt):
            raise ValueError("the permutation has length %s but its maximal "
                             "element is %s. Some element may be repeated, or "
                             "an element is missing, but there is something "
                             "wrong with its length." % (len(srt), srt[-1]))
        previous = srt[0] - 1
        for i in srt:
            if i == previous:
                raise ValueError("an element appears twice in the input")
            previous = i

    # -- list-like interface (as Sage's CombinatorialElement) ------------------

    def __len__(self):
        return len(self._list)

    def __iter__(self):
        return iter(self._list)

    def __getitem__(self, i):
        return self._list[i]

    def __call__(self, i):
        return self._list[i - 1]

    def __eq__(self, other):
        return isinstance(other, Permutation) and self._list == other._list

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(tuple(self._list))

    def __repr__(self):
        # Sage's default display option 'list'.
        return repr(self._list)

    __str__ = __repr__

    # -- group operations ------------------------------------------------------

    def inverse(self):
        w = list(range(len(self._list)))
        for i, j in enumerate(self._list, start=1):
            w[j - 1] = i
        return Permutation(w)

    __invert__ = inverse

    def right_action_product(self, rp):
        """``self`` applied first, ``rp`` afterwards (``i -> rp(self(i))``).

        Permutations of different sizes are padded with fixed points, exactly
        as ``sage.combinat.permutation_cython.right_action_product`` does.
        """
        s = list(self._list)
        r = list(rp)
        for i in range(len(s) + 1, len(r) + 1):
            s.append(i)
        for i in range(len(r) + 1, len(s) + 1):
            r.append(i)
        return Permutation([r[i - 1] for i in s])

    def left_action_product(self, lp):
        """``lp`` applied first, ``self`` afterwards (``i -> self(lp(i))``)."""
        return Permutation(lp).right_action_product(self)

    # -- cycle structure -------------------------------------------------------

    def to_cycles(self, singletons=True, use_min=True):
        """Disjoint cycles as tuples, in Sage's order.

        With ``use_min=True`` the cycles are listed by increasing smallest
        element and each cycle starts with its smallest element; with
        ``use_min=False`` by increasing largest element, each starting with its
        largest element.  Fixed points are included unless ``singletons=False``.
        """
        cycles = []
        l = list(self._list)
        groundset = range(len(l)) if use_min else reversed(range(len(l)))
        for i in groundset:
            if not l[i]:
                continue
            cycle_first = i + 1
            cycle = [cycle_first]
            l[i], nxt = False, l[i]
            while nxt != cycle_first:
                cycle.append(nxt)
                l[nxt - 1], nxt = False, l[nxt - 1]
            if singletons or len(cycle) > 1:
                cycles.append(tuple(cycle))
        return cycles

    cycle_tuples = to_cycles

    def cycle_string(self, singletons=False):
        cycles = self.to_cycles(singletons=singletons)
        if not cycles:
            return "()"
        return "".join("(" + ",".join(str(x) for x in c) + ")" for c in cycles)

    def cycle_type(self):
        """Non-increasing list of cycle lengths (a partition of ``len(self)``)."""
        return sorted((len(c) for c in self.to_cycles()), reverse=True)

    def show(self, representation="cycles", **kwargs):
        """Print the cycle decomposition (singletons included).

        Sage's ``show`` draws the permutation as a directed graph; without a
        plotting library the same information is printed as text instead.
        """
        print(self.cycle_string(singletons=True))


class PermutationGroup:
    """The group generated by some ``Permutation`` objects, acting on ``{1, ..., n}``.

    Only the pieces of Sage's ``PermutationGroup`` needed for orbit counting are
    reproduced, see the module docstring.
    """

    def __init__(self, gens):
        gens = [Permutation(g) for g in gens]
        # Sage: the domain is {1, ..., n} with n the largest size among the
        # generators (every point of a generator counts, fixed points included),
        # and at least {1}.
        n = max([1] + [len(g) for g in gens])
        padded = []
        for g in gens:
            lst = list(g) + list(range(len(g) + 1, n + 1))
            padded.append(tuple(lst))
        if not padded:
            padded = [tuple(range(1, n + 1))]
        # Sage (canonicalize=True): remove duplicates and sort the generators
        # lexicographically by their one-line notation.
        self._degree = n
        self._gens = tuple(Permutation(list(g)) for g in sorted(set(padded)))

    def gens(self):
        return self._gens

    def degree(self):
        return self._degree

    def orbits(self):
        """Orbits of the group on ``{1, ..., n}``, in GAP's ``OrbitsPerms`` order.

        Orbits are listed by increasing smallest element.  Inside an orbit the
        points appear in the order in which GAP discovers them: a breadth-first
        search from the smallest point, applying the generators in the group's
        (sorted) generator order.  Fixed points of the whole group are orbits of
        length one.  Returns a tuple of tuples, as Sage 10.5 does.

        Known exception: for a group given by a *single* non-identity generator
        that has a fixed point in ``{1, ..., n}``, GAP walks the orbits with a
        polycyclic generating sequence instead of the generator, so the order
        of the points inside an orbit can differ from the one returned here.
        The orbits themselves (and hence their number) are the same.  Neuwirth's
        algorithm only reads the point order in ``construct_F``, whose generators
        always include a fixed-point-free involution, so this case never arises
        there.
        """
        n = self._degree
        images = [tuple(g) for g in self._gens]
        seen = [False] * (n + 1)
        orbits = []
        for start in range(1, n + 1):
            if seen[start]:
                continue
            orb = [start]
            seen[start] = True
            for pnt in orb:  # orb grows while we iterate: breadth-first search
                for g in images:
                    img = g[pnt - 1]
                    if not seen[img]:
                        orb.append(img)
                        seen[img] = True
            orbits.append(tuple(orb))
        return tuple(orbits)

    def __repr__(self):
        return "Permutation Group with generators %s" % (
            [g.cycle_string() for g in self._gens],)


__all__ = ["Permutation", "PermutationGroup"]
