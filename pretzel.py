"""Heegaard diagrams on a three-dimensional genus-g surface.

The handlebody is drawn as a round body with the g handles as fat hoops around it,
blended into one smooth surface (a pretzel for g = 2, 3 when the hoops lie in the
plane of the body's equator); the alpha-curves are the belt circles of the hoops (the
meridians half way round) and the beta-curves run over the body and along the hoops as
the relators say.

The curves come from the planar diagram (``heegaard.HeegaardDiagram``).  When the disks
of the diagram sit around its middle in convex position, a homeomorphism of the plane
moves them onto a circle, the rim, and the rim becomes the equator of the body
(stereographic projection), so that the hoops lie flat; otherwise the plane is mapped
to the body as it is and the hoops rise wherever their feet are.  The planar arcs are
carried along, the strands over the hoops are computed on the hoops themselves, with
the twists that the cyclic orders of the sheets force, and everything is pushed onto
the blended surface.  The web page renders the surface with three.js, where it can be
turned around.
"""

import json
import math

import surface

FLAT = 0.72           # the picture is squashed vertically by this factor
BLEND = 0.09          # width of the transition between the body and the hoops
TUBE = 0.022          # radius of the drawn curves
LIFT = 0.012          # curves float this far above the surface
SAMPLES = 40          # points per strand over a hoop
FAT = 0.38            # radius of the hoops away from the feet
MAX_TWIST = 1.2       # strands twist at most this much (radians) around a hoop
ALPHA = "#1f2430"     # colour of the alpha-curves (the beta-curves use the palette of the planar view)


# --------------------------------------------------------------------------- vectors

def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _unit(v):
    n = math.sqrt(_dot(v, v))
    return (v[0] / n, v[1] / n, v[2] / n)


def _add(a, b, s=1.0):
    return (a[0] + s * b[0], a[1] + s * b[1], a[2] + s * b[2])


def _scale(v, s):
    return (v[0] * s, v[1] * s, v[2] * s)


# --------------------------------------------------------------------------- the plane

def _periodic_spline(thetas, values):
    """Periodic cubic spline on the circle through (thetas[k], values[k]); the angles
    are increasing in [0, 2 pi).  Returns a function of the angle."""
    n = len(thetas)
    if n == 1:
        return lambda t: values[0]
    two_pi = 2 * math.pi
    h = [((thetas[(k + 1) % n] - thetas[k]) % two_pi) or two_pi for k in range(n)]
    # second derivatives M_k from the cyclic tridiagonal system
    A = [[0.0] * n for _ in range(n)]
    b = [0.0] * n
    for k in range(n):
        A[k][(k - 1) % n] += h[k - 1]
        A[k][k] += 2 * (h[k - 1] + h[k])
        A[k][(k + 1) % n] += h[k]
        b[k] = 6 * ((values[(k + 1) % n] - values[k]) / h[k] - (values[k] - values[k - 1]) / h[k - 1])
    for i in range(n):            # Gaussian elimination, the system is tiny
        piv = max(range(i, n), key=lambda r: abs(A[r][i]))
        A[i], A[piv] = A[piv], A[i]
        b[i], b[piv] = b[piv], b[i]
        for r in range(n):
            if r != i and A[r][i] != 0.0:
                f = A[r][i] / A[i][i]
                A[r] = [A[r][j] - f * A[i][j] for j in range(n)]
                b[r] -= f * b[i]
    M = [b[i] / A[i][i] for i in range(n)]

    def f(t):
        u = (t - thetas[0]) % two_pi
        acc = 0.0
        k = 0
        while k < n - 1 and acc + h[k] <= u:
            acc += h[k]
            k += 1
        x = u - acc
        hk = h[k]
        y0, y1 = values[k], values[(k + 1) % n]
        m0, m1 = M[k], M[(k + 1) % n]
        return (m0 * (hk - x) ** 3 + m1 * x ** 3) / (6 * hk) + (y0 / hk - m0 * hk / 6) * (hk - x) + \
            (y1 / hk - m1 * hk / 6) * x
    return f


def _refine(points, max_length):
    """Subdivide the segments of a polyline until none is longer than max_length (the
    maps applied to the curves are not affine, so the points must be dense)."""
    out = [points[0]]
    for p, q in zip(points, points[1:]):
        pieces = max(1, int(math.ceil(math.hypot(q[0] - p[0], q[1] - p[1]) / max_length)))
        for k in range(1, pieces + 1):
            out.append((p[0] + (q[0] - p[0]) * k / pieces, p[1] + (q[1] - p[1]) * k / pieces))
    return out


def _untwist(plane, shifts):
    """Turn the ends of the arcs around their disks: ``shifts`` maps a dart to the angle
    its port is to be moved by.  Around every disk the plane is sheared by an angle that
    interpolates the shifts of its ports at the disk circle and fades out with the
    distance, which is a homeomorphism as long as the cyclic order of the ports is
    kept; shifts that would break it are scaled down.  Updates ``plane.curves`` and
    ``plane.port_angle`` in place."""
    d = plane.d
    two_pi = 2 * math.pi
    warps = {}
    for v in plane.centres:
        ports = sorted((plane.port_angle[e] % two_pi, shifts.get(e, 0.0)) for e in d.vertex_darts[v])
        if not ports or all(abs(s) < 1e-9 for _, s in ports):
            continue
        n = len(ports)
        factor = 1.0
        if n > 1:
            for k in range(n):
                a0, s0 = ports[k]
                a1, s1 = ports[(k + 1) % n]
                gap = (a1 - a0) % two_pi or two_pi
                slope = (s1 - s0) / gap
                if slope < -0.85:
                    factor = min(factor, 0.85 / -slope)
        support = 3.0 * plane.radii[v]
        qx, qy = plane.centres[v]
        for w in plane.centres:
            if w != v:
                support = min(support, 0.9 * (math.hypot(plane.centres[w][0] - qx, plane.centres[w][1] - qy)
                                              - plane.radii[w]))
        warps[v] = ([(a, s * factor) for a, s in ports], max(support, 1.2 * plane.radii[v]))

    def delta(ports, a):
        n = len(ports)
        if n == 1:
            return ports[0][1]
        a %= two_pi
        for k in range(n):
            a0, s0 = ports[k]
            a1, s1 = ports[(k + 1) % n]
            gap = (a1 - a0) % two_pi or two_pi
            t = (a - a0) % two_pi
            if t <= gap:
                return s0 + (s1 - s0) * t / gap
        return ports[0][1]

    def warp(p):
        for v, (ports, S) in warps.items():
            qx, qy = plane.centres[v]
            dx, dy = p[0] - qx, p[1] - qy
            rho = math.hypot(dx, dy)
            if rho >= S or rho < 1e-12:
                continue
            r = plane.radii[v]
            t = 0.0 if rho <= r else (rho - r) / (S - r)
            w = 1.0 - t * t * (3 - 2 * t)
            a = math.atan2(dy, dx) + w * delta(ports, math.atan2(dy, dx))
            p = (qx + rho * math.cos(a), qy + rho * math.sin(a))
        return p

    plane.curves = {idx: [warp(p) for p in pts] for idx, pts in plane.curves.items()}
    plane.port_angle = dict(plane.port_angle)
    for idx, pts in plane.curves.items():
        arc = d.arcs[idx]
        for end, e in ((0, arc["exit"]), (-1, arc["entry"])):
            qx, qy = plane.centres[d.vertex_of[e]]
            plane.port_angle[e] = math.atan2(pts[end][1] - qy, pts[end][0] - qx)
    plane.crossings = d._count_crossings(plane.curves)


class PlaneLayout:
    """The planar diagram as the three-dimensional picture needs it: ``curves`` (arc
    index -> polyline from exit to entry, ending on the disk circles), ``centres``,
    ``radii``, ``port_angle`` (dart -> angle of its end point around its disk), and
    ``centre``, ``scale``: the point and the length that become the pole and the
    radius of the equator of the body.  This is the diagram's own layout; ``RimLayout``
    moves the disks onto a circle first."""

    def __init__(self, diagram, layout):
        self.d = diagram
        vertices = layout["vertices"]
        self.centres = dict(layout["centres"])
        self.radii = dict(layout["radii"])
        n = len(vertices)
        self.centre = (sum(self.centres[v][0] for v in vertices) / n, sum(self.centres[v][1] for v in vertices) / n)
        self.scale = math.sqrt(sum(math.dist(self.centres[v], self.centre) ** 2 for v in vertices) / n)
        fine = max(0.15 * min(self.radii.values()), 0.008)
        self.curves = {idx: _refine(pts, fine) for idx, pts in diagram._curves(layout).items()}
        self.port_angle = {}
        for idx, pts in self.curves.items():
            arc = diagram.arcs[idx]
            for end, e in ((0, arc["exit"]), (-1, arc["entry"])):
                qx, qy = self.centres[diagram.vertex_of[e]]
                self.port_angle[e] = math.atan2(pts[end][1] - qy, pts[end][0] - qx)
        self.crossings = diagram._count_crossings(self.curves)
        self.flat = False


class RimLayout:
    """The planar diagram carried by a homeomorphism of the plane onto a picture in which
    the disks sit on a circle (the rim) with round boundaries.  Same attributes as
    ``PlaneLayout``; ``order`` lists the disks counter-clockwise around the rim and
    ``theta`` gives their angles; ``flat`` says whether the picture is a good one (the
    disks were seen at similar distances and well separated angles from the centre, and
    the two feet of every handle are neighbours on the rim)."""

    def __init__(self, diagram, layout):
        self.d = diagram
        vertices = layout["vertices"]
        centres = layout["centres"]
        radii = layout["radii"]
        self.centre = self._choose_centre(vertices, centres)
        cx, cy = self.centre

        def angle(v):
            return math.atan2(centres[v][1] - cy, centres[v][0] - cx) % (2 * math.pi)

        pairs = sorted((angle(v), v) for v in vertices)
        thetas = [t for t, v in pairs]
        self.order = [v for t, v in pairs]
        dist = [math.hypot(centres[v][0] - cx, centres[v][1] - cy) for v in self.order]
        self.scale = math.exp(sum(math.log(x) for x in dist) / len(dist))
        self.rho = _periodic_spline(thetas, [math.log(x) for x in dist])   # log radius of the rim curve
        self.theta = {v: t for v, t in zip(self.order, thetas)}
        n = len(vertices)
        gaps = [(thetas[(k + 1) % n] - thetas[k]) % (2 * math.pi) for k in range(n)]
        self.flat = (min(gaps) >= 0.5 * 2 * math.pi / n and max(dist) <= 2.6 * min(dist) and
                     all(abs(self.order.index(2 * i) - self.order.index(2 * i + 1)) in (1, n - 1)
                         for i in range(diagram.g) if 2 * i in self.theta and 2 * i + 1 in self.theta))

        # 1. the radial map that straightens the rim curve into a circle
        fine = max(0.15 * min(radii.values()), 0.008)
        curves = {idx: _refine(pts, fine) for idx, pts in diagram._curves(layout).items()}
        curves = {idx: [self._radial(p) for p in pts] for idx, pts in curves.items()}
        self.centres = {v: self._radial(centres[v]) for v in vertices}
        # the image of every disk boundary, as a radius function around the new centre
        self.images = {}
        for v in vertices:
            qx, qy = self.centres[v]
            samples = []
            for k in range(72):
                a = 2 * math.pi * k / 72
                p = self._radial((centres[v][0] + radii[v] * math.cos(a), centres[v][1] + radii[v] * math.sin(a)))
                samples.append((math.atan2(p[1] - qy, p[0] - qx) % (2 * math.pi), math.hypot(p[0] - qx, p[1] - qy)))
            samples.sort()
            self.images[v] = samples
        # 2. every hoop gets the radius of its smaller foot; the warp shrinks the disks
        # (radially, around their centres) to exactly that circle
        self.radii = {}
        for v in vertices:
            partner = v ^ 1
            feet = [v, partner] if partner in self.images else [v]
            self.radii[v] = min(r for w in feet for _, r in self.images[w])
        self.support = {}
        for v in vertices:
            qx, qy = self.centres[v]
            s = 4.0 * self.radii[v]
            for w in vertices:
                if w != v:
                    gap = math.hypot(self.centres[w][0] - qx, self.centres[w][1] - qy) - max(r for _, r in self.images[w])
                    s = min(s, 0.9 * gap)
            self.support[v] = max(s, 1.3 * max(r for _, r in self.images[v]))
        self.curves = {}
        self.port_angle = {}
        for idx, pts in curves.items():
            arc = diagram.arcs[idx]
            out = [self._warp(p) for p in pts]
            for end, e in ((0, arc["exit"]), (-1, arc["entry"])):
                v = diagram.vertex_of[e]
                qx, qy = self.centres[v]
                a = math.atan2(out[end][1] - qy, out[end][0] - qx)
                out[end] = (qx + self.radii[v] * math.cos(a), qy + self.radii[v] * math.sin(a))
                self.port_angle[e] = a
            self.curves[idx] = out
        self.crossings = diagram._count_crossings(self.curves)

    def _choose_centre(self, vertices, centres):
        """A point from which the disks are seen in well separated directions and at
        similar distances: the best of a grid around the centroid."""
        n = len(vertices)
        mx = sum(centres[v][0] for v in vertices) / n
        my = sum(centres[v][1] for v in vertices) / n
        if n == 1:
            return (mx, my)
        if n == 2:
            # a torus: the two feet 100 degrees apart on the rim
            (ax, ay), (bx, by) = (centres[v] for v in vertices)
            d = math.hypot(bx - ax, by - ay)
            nx, ny = -(by - ay) / d, (bx - ax) / d
            off = 0.5 * d / math.tan(math.radians(50))
            return (mx + off * nx, my + off * ny)
        spread = max(math.hypot(centres[v][0] - mx, centres[v][1] - my) for v in vertices)

        def score(c):
            angles = sorted(math.atan2(centres[v][1] - c[1], centres[v][0] - c[0]) for v in vertices)
            gaps = [(angles[(k + 1) % n] - angles[k]) % (2 * math.pi) for k in range(n)]
            dist = [math.hypot(centres[v][0] - c[0], centres[v][1] - c[1]) for v in vertices]
            logs = [math.log(max(x, 1e-9)) for x in dist]
            mean = sum(logs) / n
            return (max(abs(l - mean) for l in logs) + 2.0 * max(0.0, 0.5 * 2 * math.pi / n - min(gaps))
                    + 0.3 * math.hypot(c[0] - mx, c[1] - my) / spread)

        best = ((mx, my), score((mx, my)))
        for i in range(-6, 7):
            for j in range(-6, 7):
                c = (mx + 0.5 * spread * i / 6, my + 0.5 * spread * j / 6)
                s = score(c)
                if s < best[1]:
                    best = (c, s)
        return best[0]

    def _radial(self, p):
        cx, cy = self.centre
        dx, dy = p[0] - cx, p[1] - cy
        r = math.hypot(dx, dy)
        if r < 1e-12:
            return (cx, cy)
        t = math.atan2(dy, dx)
        scale = self.scale / math.exp(self.rho(t))
        return (cx + dx * scale, cy + dy * scale)

    def _image_radius(self, v, a):
        samples = self.images[v]
        a %= 2 * math.pi
        n = len(samples)
        for k in range(n):
            a0, r0 = samples[k]
            a1, r1 = samples[(k + 1) % n]
            span = (a1 - a0) % (2 * math.pi)
            t = (a - a0) % (2 * math.pi)
            if t <= span and span > 0:
                return r0 + (r1 - r0) * t / span
        return samples[0][1]

    def _warp(self, p):
        """Shrink every disk image to its circle; the plane around it follows (less and
        less, up to the support radius), so the map is a homeomorphism."""
        for v, (qx, qy) in self.centres.items():
            dx, dy = p[0] - qx, p[1] - qy
            r = math.hypot(dx, dy)
            S = self.support[v]
            if r >= S or r < 1e-12:
                continue
            a = math.atan2(dy, dx)
            r_img = self._image_radius(v, a)
            lam = self.radii[v] / r_img
            if r <= r_img:
                scale = lam
            else:
                t = (r - r_img) / (S - r_img)
                scale = lam ** (1.0 - t)
            p = (qx + dx * scale, qy + dy * scale)
        return p


# --------------------------------------------------------------------------- the surface

class Hoop:
    """A fat tube along a circular arch from the foot P to the foot Q (unit vectors, on
    the body): the arc of the circle through P and Q whose centre lies on their
    bisector, in the plane through the centre of the body (or that plane turned about
    the chord PQ by ``tilt``), reaching out to the distance ``far``.  The tube has the
    radius of the feet at the feet and widens to ``r`` along a short neck."""

    def __init__(self, P, Q, r_foot, r=None, far=None, tilt=0.0):
        self.P, self.Q = P, Q
        self.r_foot = r_foot
        self.tilt = tilt
        cosb = max(-1.0, min(1.0, _dot(P, Q)))
        beta = math.acos(cosb) / 2                       # half the angle between the feet
        self.beta = beta
        M = _scale(_add(P, Q), 0.5)
        b = _unit(M) if _dot(M, M) > 1e-12 else _unit(_cross(P, (0.0, 0.0, 1.0)))
        s = _unit(_add(P, Q, -1.0))                       # from the midpoint towards P
        m = _cross(b, s)                                  # normal of the plane of the hoop
        if tilt:
            c, t = math.cos(tilt), math.sin(tilt)
            b, m = _add(_scale(b, c), m, t), _add(_scale(m, c), b, -t)
        F = far if far is not None else 1.0 + 1.05 * 2 * math.sin(beta)
        D = (F * F - 1.0) / (2 * (F - math.cos(beta)))
        rho = F - D
        self.r = max(r_foot, min(r if r is not None else FAT, 0.55 * rho))
        psi = math.atan2(math.sin(beta), math.cos(beta) - D)
        self.C = _add(M, b, D - math.cos(beta))
        self.e1, self.e2, self.m = b, s, m
        self.rho, self.psi = rho, psi
        self.length = rho * 2 * psi
        self.neck = min(0.5, 0.3 * self.length)

    def radius(self, s):
        """Radius of the tube at arc length s: the feet's radius at both ends, r in
        the middle."""
        t = min(s, self.length - s) / self.neck
        if t >= 1.0:
            return self.r
        t = max(t, 0.0)
        return self.r_foot + (self.r - self.r_foot) * t * t * (3 - 2 * t)

    def _at(self, s):
        """Centre line point, unit tangent and the two cross-section directions at arc
        length s (u1 = -m, u2 = m x t)."""
        s = max(0.0, min(self.length, s))
        a = self.psi - s / self.rho
        ca, sa = math.cos(a), math.sin(a)
        p = _add(_add(self.C, self.e1, self.rho * ca), self.e2, self.rho * sa)
        t = _add(_scale(self.e1, sa), self.e2, -ca)
        return p, t, _scale(self.m, -1.0), _cross(self.m, t)

    def centre(self, s):
        return self._at(s)[0]

    def point(self, s, phi):
        """Point of the tube at arc length s and angle phi around it."""
        p, t, u1, u2 = self._at(s)
        r = self.radius(s)
        return _add(_add(p, u1, r * math.cos(phi)), u2, r * math.sin(phi))

    def angle_of(self, s, direction):
        """The angle phi of a direction (projected onto the cross-section at s)."""
        p, t, u1, u2 = self._at(s)
        return math.atan2(_dot(direction, u2), _dot(direction, u1))

    def distance(self, p):
        """(distance from p to the centre line of the hoop, arc length of the nearest
        point)."""
        d = _add(p, self.C, -1.0)
        x, y, z = _dot(d, self.e1), _dot(d, self.e2), _dot(d, self.m)
        a = math.atan2(y, x)
        if -self.psi <= a <= self.psi:
            return math.hypot(math.hypot(x, y) - self.rho, z), (self.psi - a) * self.rho
        best = None
        for sgn, s in ((1.0, 0.0), (-1.0, self.length)):
            dd = math.sqrt((x - self.rho * math.cos(sgn * self.psi)) ** 2 +
                           (y - self.rho * math.sin(sgn * self.psi)) ** 2 + z * z)
            if best is None or dd < best[0]:
                best = (dd, s)
        return best


class Pretzel:
    """The surface and the curves of the three-dimensional picture of a diagram."""

    def __init__(self, diagram, layout=None):
        d = self.d = diagram
        if d.L != 1 or d.g > 3:
            raise ValueError("the three-dimensional picture needs L = 1 and genus at most 3")
        layout = layout or d.layouts[0]
        plane = RimLayout(d, layout)
        if not plane.flat or plane.crossings:
            plane = PlaneLayout(d, layout)
        if plane.crossings:
            raise ValueError("the layout has crossings")
        self.plane = plane
        self.flat = plane.flat
        # the feet and the hoops, one per generator
        self.feet = {v: self.sphere(plane.centres[v]) for v in plane.centres}
        self.hoops = []
        for i in range(d.g):
            a, abar = 2 * i, 2 * i + 1
            self.hoops.append(Hoop(self.feet[a], self.feet[abar], plane.radii[a] / plane.scale))
        self._separate_hoops()
        self.alpha = [self._belt(i) for i in range(d.g)]
        # the strands twist as the ports at the two feet dictate; most of the twist is
        # moved into the plane (the arcs curl around the feet instead), which is the
        # same curve up to isotopy and far easier to read on the hoops
        self.strands = [self._strands(i, sample=False) for i in range(d.g)]
        shifts = {}
        for i in range(d.g):
            for e, st in self.strands[i].items():
                excess = st["twist"] - max(-MAX_TWIST, min(MAX_TWIST, st["twist"]))
                if abs(excess) > 1e-9:
                    shifts[e] = 0.5 * excess
                    shifts[d.B(e)] = 0.5 * excess
        if shifts:
            saved = (plane.curves, plane.port_angle)
            _untwist(plane, shifts)
            if plane.crossings:
                # the arcs could not be turned (tiny disks, say): leave the twist on the hoops
                plane.curves, plane.port_angle = saved
                plane.crossings = 0
        self.strands = [self._strands(i) for i in range(d.g)]
        self.arcs = {idx: self._arc(idx, pts) for idx, pts in plane.curves.items()}
        self.beta = [self._beta(r) for r in range(len(d.rels))]

    # -- geometry
    def sphere(self, p):
        """Point of the unit sphere for a point of the plane (stereographic projection
        from the south pole: the circle of radius ``scale`` around ``centre`` is the
        equator, its inside the top)."""
        (cx, cy), R = self.plane.centre, self.plane.scale
        x, y = (p[0] - cx) / R, (p[1] - cy) / R
        n = 1 + x * x + y * y
        return (2 * x / n, 2 * y / n, (1 - x * x - y * y) / n)

    def sphere_direction(self, p, angle, h=1e-5):
        """Unit tangent vector of the sphere in which a planar direction at p points."""
        q = (p[0] + h * math.cos(angle), p[1] + h * math.sin(angle))
        return _unit(_add(self.sphere(q), self.sphere(p), -1.0))

    def _separate_hoops(self):
        """Tilt hoops that would run into each other out of their planes."""
        for i in range(len(self.hoops)):
            for attempt in range(4):
                if all(self._hoops_apart(i, j) for j in range(i)):
                    break
                h = self.hoops[i]
                tilt = (0.55 if attempt % 2 == 0 else -0.55) * (attempt // 2 + 1)
                self.hoops[i] = Hoop(h.P, h.Q, h.r_foot, tilt=tilt)

    def _hoops_apart(self, i, j):
        hi, hj = self.hoops[i], self.hoops[j]
        for k in range(1, 24):
            p = hi.centre(hi.length * k / 24)
            if hj.distance(p)[0] < 1.15 * (hi.r + hj.r):
                return False
        return True

    def field(self, p):
        """Blend of the body and the hoops; the surface is its zero set."""
        f = 1.0 / (1.0 + math.exp(min(50.0, (math.sqrt(_dot(p, p)) - 1.0) / BLEND)))
        for h in self.hoops:
            d, s = h.distance(p)
            f += 1.0 / (1.0 + math.exp(min(50.0, (d - h.radius(s)) / BLEND)))
        return f - 0.5

    def gradient(self, p, h=1e-4):
        f0 = self.field(p)
        return tuple((self.field((p[0] + (k == 0) * h, p[1] + (k == 1) * h, p[2] + (k == 2) * h)) - f0) / h
                     for k in range(3))

    def project(self, p, lift=LIFT):
        """Push p onto the surface (Newton steps along the gradient), then lift it a
        little above it."""
        x, y, z = p
        for _ in range(40):
            f = self.field((x, y, z))
            if abs(f) < 1e-6:
                break
            gx, gy, gz = self.gradient((x, y, z))
            g2 = gx * gx + gy * gy + gz * gz
            if g2 < 1e-30:
                # deep inside or far outside, where the field is flat: head for the body
                n = math.sqrt(x * x + y * y + z * z) + 1e-12
                step = 0.2 if f > 0 else -0.2
                x, y, z = x + step * x / n, y + step * y / n, z + step * z / n
                continue
            # a Newton step, never longer than a short stride (the field saturates)
            t = min(abs(f) / g2, 0.2 / math.sqrt(g2)) * (1 if f > 0 else -1)
            x, y, z = x - t * gx, y - t * gy, z - t * gz
        g = self.gradient((x, y, z))
        if _dot(g, g) < 1e-18:
            return (x, y, z)
        return _add((x, y, z), _unit(g), lift)

    # -- curves
    def _belt(self, i):
        h = self.hoops[i]
        pts = []
        for k in range(SAMPLES + 1):
            phi = 2 * math.pi * k / SAMPLES
            pts.append(self.project(h.point(h.length / 2, phi)))
        return {"generator": self.d.gens[i][0], "points": pts,
                "label_at": self.project(h.point(h.length / 2, math.pi), lift=0.16)}

    def _strands(self, i, sample=True):
        d = self.d
        h = self.hoops[i]
        plane = self.plane
        strands = []
        for e in d.vertex_darts[2 * i]:
            e_bar = d.B(e)
            dir_p = self.sphere_direction(plane.centres[2 * i], plane.port_angle[e])
            dir_q = self.sphere_direction(plane.centres[2 * i + 1], plane.port_angle[e_bar])
            phi_p = h.angle_of(0.0, dir_p) % (2 * math.pi)
            phi_q = h.angle_of(h.length, dir_q) % (2 * math.pi)
            strands.append({"sheet": d.dart_info[e][1], "phi_p": phi_p, "phi_q": phi_q, "dart": e})
        strands.sort(key=lambda s: s["phi_p"])
        if strands:
            base_p, base_q = strands[0]["phi_p"], strands[0]["phi_q"]
            for s in strands:
                s["lift_p"] = base_p + ((s["phi_p"] - base_p) % (2 * math.pi))
                s["lift_q"] = base_q + ((s["phi_q"] - base_q) % (2 * math.pi))
            best = None
            for m in range(-2, 3):
                score = max(abs(s["lift_q"] + 2 * math.pi * m - s["lift_p"]) for s in strands)
                if best is None or score < best[0]:
                    best = (score, m)
            for s in strands:
                s["twist"] = s["lift_q"] + 2 * math.pi * best[1] - s["lift_p"]
                if not sample:
                    continue
                pts = []
                for k in range(SAMPLES + 1):
                    t = k / SAMPLES
                    phi = s["lift_p"] + s["twist"] * (1 - math.cos(math.pi * t)) / 2
                    pts.append(self.project(h.point(h.length * t, phi)))
                s["points"] = pts
        return {s["dart"]: s for s in strands}

    def _strand_end(self, e):
        """The point where the strand of the dart e meets the body."""
        i, sheet, barred = self.d.dart_info[e]
        s = self.strands[i][self.d.B(e) if barred else e]
        return s["points"][-1] if barred else s["points"][0]

    def _arc(self, idx, pts):
        """The arc on the body, its ends bent onto the ends of the strands it continues
        (the tubes need not meet the body exactly at the disks)."""
        d = self.d
        arc = d.arcs[idx]
        out = [self.sphere(p) for p in pts]
        for end, e in ((0, arc["exit"]), (-1, arc["entry"])):
            v = d.vertex_of[e]
            qx, qy = self.plane.centres[v]
            reach = 3.0 * self.plane.radii[v]
            target = self._strand_end(e)
            delta = _add(target, self.project(out[end]), -1.0)
            for k, p in enumerate(pts):
                dist = math.hypot(p[0] - qx, p[1] - qy) - self.plane.radii[v]
                if dist >= reach:
                    continue
                t = 1.0 - dist / reach
                out[k] = _add(out[k], delta, t * t * (3 - 2 * t))
        out = [self.project(p) for p in out]
        out[0] = self._strand_end(arc["exit"])
        out[-1] = self._strand_end(arc["entry"])
        return out

    def _beta(self, r):
        """The curve of relator r: arcs and strands in the order of the letters."""
        d = self.d
        pieces = []
        for idx, arc in enumerate(d.arcs):
            if arc["rel"] != r:
                continue
            pieces.append({"kind": "arc", "index": idx, "points": self.arcs[idx], "closing": arc["closing"]})
            e = arc["entry"]
            i, sheet, barred = d.dart_info[e]
            unbarred = d.B(e) if barred else e
            s = self.strands[i][unbarred]
            pts = s["points"] if not barred else s["points"][::-1]
            pieces.append({"kind": "strand", "generator": i, "sheet": sheet, "points": pts})
        return pieces

    def mesh(self, h=0.15):
        """A triangle mesh of the surface (marching tetrahedra on a grid of spacing h), as
        (vertices, triangles); the viewer builds its own, this one is for checking."""
        reach = self.bounds()
        return _marching_tetrahedra(self.field, ((-reach, reach), (-reach, reach), (-reach, reach)), h)

    def genus(self, h=0.15):
        """The genus of the surface, from the Euler characteristic of ``mesh``."""
        vertices, triangles = self.mesh(h)
        edges = {(min(a, b), max(a, b)) for t in triangles for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0]))}
        used = {v for t in triangles for v in t}
        return (2 - (len(used) - len(edges) + len(triangles))) // 2

    # -- export
    def bounds(self):
        reach = 1.0
        for h in self.hoops:
            for k in range(33):
                p = h.centre(h.length * k / 32)
                reach = max(reach, abs(p[0]), abs(p[1]), abs(p[2]))
        return reach + max(h.r for h in self.hoops) + 0.25

    def scene(self):
        """Everything the viewer needs, as plain data."""
        d = self.d
        label_arcs = {d.arcs[idx]["rel"]: idx for idx in d._label_arcs()}
        hoops = []
        for h in self.hoops:
            hoops.append({"r": h.r, "r_foot": h.r_foot, "neck": h.neck, "length": h.length, "c": h.C,
                          "e1": h.e1, "e2": h.e2, "m": h.m, "rho": h.rho, "psi": h.psi})
        palette = surface.SurfaceDrawing.PALETTE
        betas = []
        for r, pieces in enumerate(self.beta):
            betas.append({"colour": palette[r % len(palette)], "label": "β<sub>%d</sub>" % (r + 1),
                          "pieces": pieces, "label_arc": label_arcs.get(r)})
        alphas = [{"colour": ALPHA, "label": "α<sub>%s</sub>" % a["generator"],
                   "points": a["points"], "label_at": a["label_at"]} for a in self.alpha]
        return {"flat": FLAT, "blend": BLEND, "tube": TUBE, "hoops": hoops, "alpha": alphas, "beta": betas,
                "bounds": self.bounds()}

    def to_html_fragment(self, element_id="pretzel", fallback=""):
        """The viewer: a <div> with the WebGL canvas and a module script that renders the
        scene with three.js (loaded from a CDN).  ``fallback`` is HTML shown inside the
        box until the script has drawn the surface, and for good if it cannot (no
        network, no WebGL)."""
        return (_VIEWER_HTML.replace("__ID__", element_id).replace("__SCENE__", json.dumps(self.scene()))
                .replace("__FALLBACK__", fallback))


_TET_SPLIT = [(0, 1, 3, 7), (0, 1, 5, 7), (0, 4, 5, 7), (0, 4, 6, 7), (0, 2, 6, 7), (0, 2, 3, 7)]


def _marching_tetrahedra(field, bounds, h):
    """Triangle mesh of the zero set of ``field`` on a grid of spacing h: (vertices,
    triangles) with shared vertices."""
    (x0, x1), (y0, y1), (z0, z1) = bounds
    nx, ny, nz = int((x1 - x0) / h) + 1, int((y1 - y0) / h) + 1, int((z1 - z0) / h) + 1
    values = {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                values[(i, j, k)] = field((x0 + i * h, y0 + j * h, z0 + k * h))
    vertices = []
    edge_vertex = {}
    triangles = []

    def vertex_on_edge(a, b):
        key = (a, b) if a < b else (b, a)
        if key not in edge_vertex:
            fa, fb = values[a], values[b]
            t = fa / (fa - fb)
            pa = (x0 + a[0] * h, y0 + a[1] * h, z0 + a[2] * h)
            pb = (x0 + b[0] * h, y0 + b[1] * h, z0 + b[2] * h)
            vertices.append((pa[0] + t * (pb[0] - pa[0]), pa[1] + t * (pb[1] - pa[1]), pa[2] + t * (pb[2] - pa[2])))
            edge_vertex[key] = len(vertices) - 1
        return edge_vertex[key]

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                cube = [(i, j, k), (i + 1, j, k), (i, j + 1, k), (i + 1, j + 1, k),
                        (i, j, k + 1), (i + 1, j, k + 1), (i, j + 1, k + 1), (i + 1, j + 1, k + 1)]
                vals = [values[c] for c in cube]
                if all(v > 0 for v in vals) or all(v <= 0 for v in vals):
                    continue
                for tet in _TET_SPLIT:
                    pts = [cube[t] for t in tet]
                    inside = [p for p in pts if values[p] > 0]
                    outside = [p for p in pts if values[p] <= 0]
                    if not inside or not outside:
                        continue
                    if len(inside) == 1 or len(outside) == 1:
                        apex = inside[0] if len(inside) == 1 else outside[0]
                        others = outside if len(inside) == 1 else inside
                        triangles.append(tuple(vertex_on_edge(apex, o) for o in others))
                    else:
                        a, b = inside
                        c, d = outside
                        v_ac, v_ad, v_bc, v_bd = (vertex_on_edge(a, c), vertex_on_edge(a, d),
                                                  vertex_on_edge(b, c), vertex_on_edge(b, d))
                        triangles.append((v_ac, v_ad, v_bd))
                        triangles.append((v_ac, v_bd, v_bc))
    return vertices, triangles


_VIEWER_HTML = r'''
<div id="__ID__" class="pretzel-viewer" style="position:relative;width:100%;max-width:980px;height:560px;background:#f7f7f9;border:1px solid #ddd;border-radius:6px;overflow:hidden">
<div id="__ID__-fallback" style="position:absolute;inset:0;overflow:auto">__FALLBACK__</div>
</div>
<script type="importmap">
{ "imports": { "three": "https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js",
               "three/addons/": "https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/" } }
</script>
<script type="module">
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { MarchingCubes } from 'three/addons/objects/MarchingCubes.js';
import { CSS2DRenderer, CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';

const data = __SCENE__;
const container = document.getElementById('__ID__');
const fallback = document.getElementById('__ID__-fallback');
if (fallback) fallback.remove();
const hint = document.createElement('div');
hint.textContent = 'drag to turn, scroll to zoom, right-drag to pan';
hint.style.cssText = 'position:absolute;left:10px;bottom:8px;font:12px sans-serif;color:#666;pointer-events:none';
container.appendChild(hint);
const W = container.clientWidth, H = container.clientHeight;
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(window.devicePixelRatio);
renderer.setSize(W, H);
container.appendChild(renderer.domElement);
const labelRenderer = new CSS2DRenderer();
labelRenderer.setSize(W, H);
labelRenderer.domElement.style.position = 'absolute';
labelRenderer.domElement.style.top = '0';
labelRenderer.domElement.style.pointerEvents = 'none';
container.appendChild(labelRenderer.domElement);

const scene = new THREE.Scene();
scene.background = new THREE.Color(0xf7f7f9);
const camera = new THREE.PerspectiveCamera(36, W / H, 0.1, 100);
const reach = data.bounds;
camera.position.set(0, -2.3 * reach, 1.6 * reach);
camera.up.set(0, 0, 1);
const controls = new OrbitControls(camera, renderer.domElement);
controls.target.set(0, 0, 0);
controls.enableDamping = true;
scene.add(new THREE.HemisphereLight(0xffffff, 0x8899aa, 0.9));
const sun = new THREE.DirectionalLight(0xffffff, 1.4);
sun.position.set(2, -3, 5);
scene.add(sun);
const back = new THREE.DirectionalLight(0xffffff, 0.5);
back.position.set(-3, 2, -2);
scene.add(back);

// ---- the surface: marching cubes over the blended field
const flat = data.flat, blend = data.blend;
function hoopDistance(h, p) {
  const dx = p[0] - h.c[0], dy = p[1] - h.c[1], dz = p[2] - h.c[2];
  const x = dx * h.e1[0] + dy * h.e1[1] + dz * h.e1[2];
  const y = dx * h.e2[0] + dy * h.e2[1] + dz * h.e2[2];
  const z = dx * h.m[0] + dy * h.m[1] + dz * h.m[2];
  const a = Math.atan2(y, x);
  if (-h.psi <= a && a <= h.psi) return [Math.hypot(Math.hypot(x, y) - h.rho, z), (h.psi - a) * h.rho];
  let best = Infinity, sbest = 0;
  for (const [sgn, s] of [[1, 0], [-1, h.length]]) {
    const d = Math.hypot(x - h.rho * Math.cos(sgn * h.psi), y - h.rho * Math.sin(sgn * h.psi), z);
    if (d < best) { best = d; sbest = s; }
  }
  return [best, sbest];
}
function hoopRadius(h, s) {
  let t = Math.min(s, h.length - s) / h.neck;
  if (t >= 1) return h.r;
  t = Math.max(t, 0);
  return h.r_foot + (h.r - h.r_foot) * t * t * (3 - 2 * t);
}
function field(p) {
  let f = 1 / (1 + Math.exp(Math.min(50, (Math.hypot(p[0], p[1], p[2]) - 1) / blend)));
  for (const h of data.hoops) {
    const [d, s] = hoopDistance(h, p);
    f += 1 / (1 + Math.exp(Math.min(50, (d - hoopRadius(h, s)) / blend)));
  }
  return f;
}
const RES = 96;
const material = new THREE.MeshStandardMaterial({ color: 0xd9d9de, roughness: 0.55, metalness: 0.05 });
const mc = new MarchingCubes(RES, material, false, false, 400000);
const hx = reach, hy = reach, hz = reach;
for (let k = 0; k < RES; k++) {
  const z = ((k - RES / 2) / (RES / 2)) * hz;
  for (let j = 0; j < RES; j++) {
    const y = ((j - RES / 2) / (RES / 2)) * hy;
    for (let i = 0; i < RES; i++) {
      const x = ((i - RES / 2) / (RES / 2)) * hx;
      mc.field[k * RES * RES + j * RES + i] = field([x, y, z]);
    }
  }
}
mc.isolation = 0.5;
mc.update();
mc.scale.set(hx, hy, hz * flat);
scene.add(mc);

// ---- curves
function tube(points, colour, radius) {
  const pts = points.map(p => new THREE.Vector3(p[0], p[1], p[2] * flat));
  const curve = new THREE.CatmullRomCurve3(pts, false, 'centripetal');
  const geometry = new THREE.TubeGeometry(curve, Math.max(8, pts.length * 2), radius, 8, false);
  scene.add(new THREE.Mesh(geometry, new THREE.MeshStandardMaterial({ color: colour, roughness: 0.5 })));
  return curve;
}
function label(text, p, colour) {
  const div = document.createElement('div');
  div.innerHTML = text;
  div.style.cssText = 'font:italic 17px serif;color:' + colour + ';background:rgba(255,255,255,0.75);padding:0 3px;border-radius:3px';
  const obj = new CSS2DObject(div);
  obj.position.set(p[0], p[1], p[2] * flat);
  scene.add(obj);
}
function cone(curve, t, colour) {
  const p = curve.getPointAt(t), d = curve.getTangentAt(t);
  const mesh = new THREE.Mesh(new THREE.ConeGeometry(2.4 * data.tube, 6 * data.tube, 10),
                              new THREE.MeshStandardMaterial({ color: colour }));
  mesh.position.copy(p);
  mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), d.normalize());
  scene.add(mesh);
}
function dot(p, colour) {
  const mesh = new THREE.Mesh(new THREE.SphereGeometry(2.3 * data.tube, 12, 12), new THREE.MeshStandardMaterial({ color: colour }));
  mesh.position.set(p[0], p[1], p[2] * flat);
  scene.add(mesh);
}
for (const a of data.alpha) {
  tube(a.points, a.colour, data.tube * 0.9);
  label(a.label, a.label_at, a.colour);
}
for (const b of data.beta) {
  for (const piece of b.pieces) {
    const curve = tube(piece.points, b.colour, data.tube);
    if (piece.kind === 'arc') {
      cone(curve, 0.5, b.colour);
      if (piece.closing) dot(piece.points[Math.floor(piece.points.length * 0.75)], b.colour);
      if (piece.index === b.label_arc) {
        const p = piece.points[Math.floor(piece.points.length / 2)];
        label(b.label, [p[0] * 1.12, p[1] * 1.12, p[2] + 0.12], b.colour);
      }
    }
  }
}

function animate() {
  requestAnimationFrame(animate);
  controls.update();
  renderer.render(scene, camera);
  labelRenderer.render(scene, camera);
}
animate();
window.__setView = (v) => {
  camera.position.set(v[0] * reach, v[1] * reach, v[2] * reach);
  if (v.length > 3) controls.target.set(v[3], v[4], v[5]);
  controls.update();
};
window.addEventListener('resize', () => {
  const w = container.clientWidth, h = container.clientHeight;
  camera.aspect = w / h; camera.updateProjectionMatrix();
  renderer.setSize(w, h); labelRenderer.setSize(w, h);
});
</script>
'''
