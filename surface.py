"""The Heegaard diagram drawn on the genus-g surface itself.

The handlebody is drawn in the classical way: a ball with g handles attached to it.
The planar diagram of ``heegaard.py`` is the front of the ball with the 2g attaching
disks of the handles; here every pair of disks ``a_i``, ``a_i-bar`` is joined by a tube
(the i-th 1-handle), drawn in an oblique projection in which a point at height z above
the ball appears displaced by z times the vector Z on the page.  The belt sphere of the handle, the curve
alpha_i, is the circle around the middle of the tube; it appears as an ellipse, solid
on the side facing the viewer and dashed behind.  A beta-curve is made of its arcs on
the ball (the arcs of the planar diagram) and, for every letter of the relator, a
strand that climbs the tube from the point where the arc enters the disk, crosses the
belt circle once and comes down to the disk on the other side.  The two disks a_i and
a_i-bar are the two ends of the tube, so the strand is determined by the diagram; it
twists around the tube as much as the positions of its two ends require.  Parts of the
surface hidden behind a tube are drawn dashed.

When Neuwirth's group <A, C> has L > 1 orbits, the ball is really L balls joined by
tubes in the regions that the algorithm chose (the starred regions of the planar
diagram); those tubes carry no belt circle and no beta-curve.
"""

import html
import math

from heegaard import _chaikin, _polyline_point

Z = (0.34, 0.40)  # oblique projection: a point at height z above the ball appears at + z*Z
SAMPLES = 48      # points per strand / rail along a tube


def _rotate(p, theta):
    c, s = math.cos(theta), math.sin(theta)
    return (c * p[0] - s * p[1], s * p[0] + c * p[1])


def _smooth_without_crossings(polylines, count_crossings):
    """Smoothed versions of the polylines; an arc whose smoothed version would cross
    something keeps its polyline."""
    smoothed = {idx: _chaikin(pts) for idx, pts in polylines.items()}
    current = dict(smoothed)
    for _ in range(len(polylines)):
        if count_crossings(current) == 0:
            return current
        # find an arc involved in a crossing and fall back to its polyline
        changed = False
        for idx in list(current):
            if current[idx] is polylines[idx]:
                continue
            trial = dict(current)
            trial[idx] = polylines[idx]
            if count_crossings(trial) < count_crossings(current):
                current = trial
                changed = True
                break
        if not changed:
            break
    return current if count_crossings(current) == 0 else dict(polylines)


class Tube:
    """A tube rising from the ball at the disk centred at P and landing at Q (an arch),
    of radius r and height h, in the oblique projection.

    Points of the tube are parametrised by the position s in [0, pi] along the arch
    (s = 0 at P, s = pi at Q) and the angle phi around the tube.  At the foot P the
    direction phi = 0 points away from Q and angles increase counter-clockwise on the
    page; at the foot Q the direction phi = 0 points away from P and angles increase
    counter-clockwise as well.  That is what a tube does: the two attaching circles are
    identified by an orientation-reversing map, exactly as Neuwirth's C requires."""

    def __init__(self, P, Q, r, h):
        self.P, self.Q, self.r, self.h = P, Q, r, h
        dx, dy = Q[0] - P[0], Q[1] - P[1]
        self.span = math.hypot(dx, dy) / 2
        self.u = (dx / (2 * self.span), dy / (2 * self.span))
        self.n = (-self.u[1], self.u[0])
        self.M = ((P[0] + Q[0]) / 2, (P[1] + Q[1]) / 2)
        self.theta_u = math.atan2(self.u[1], self.u[0])
        self.zu = Z[0] * self.u[0] + Z[1] * self.u[1]
        self.zn = Z[0] * self.n[0] + Z[1] * self.n[1]

    def axis(self, s):
        R, h, u, M = self.span, self.h, self.u, self.M
        return (M[0] - R * math.cos(s) * u[0] + h * math.sin(s) * Z[0],
                M[1] - R * math.cos(s) * u[1] + h * math.sin(s) * Z[1])

    def radial(self, s):
        """Projection of the outward radial direction of the arch at s."""
        u = self.u
        return (-math.cos(s) * u[0] + math.sin(s) * Z[0], -math.cos(s) * u[1] + math.sin(s) * Z[1])

    def point(self, s, phi):
        a = self.axis(s)
        rad = self.radial(s)
        r, n = self.r, self.n
        return (a[0] + r * math.cos(phi) * rad[0] + r * math.sin(phi) * n[0],
                a[1] + r * math.cos(phi) * rad[1] + r * math.sin(phi) * n[1])

    def facing(self, s, phi):
        """Positive where the tube's surface faces the viewer.  The normal of the tube is
        cos(phi) (-cos s u + sin s z) + sin(phi) n and the viewer sits along (-Z, 1)."""
        return math.cos(phi) * (math.sin(s) + math.cos(s) * self.zu) - math.sin(phi) * self.zn

    def silhouette(self, s):
        """The two angles phi at which the tube's outline runs, at position s."""
        phi = math.atan2(math.sin(s) + math.cos(s) * self.zu, self.zn)
        return phi, phi + math.pi

    def outline(self):
        """Closed polygon (list of points) bounding the projection of the tube: the two
        silhouette rails, closed at each foot by the half of the foot circle that lies on
        the side away from the direction Z in which the tube rises."""
        P, Q, r = self.P, self.Q, self.r
        zang = math.atan2(Z[1], Z[0])

        def ang(C, p):
            return math.atan2(p[1] - C[1], p[0] - C[0])

        def near(a, b):
            return abs((a - b + math.pi) % (2 * math.pi) - math.pi) < 0.5

        rails = ([], [])
        for k in range(SAMPLES + 1):
            s = math.pi * k / SAMPLES
            a, b = self.silhouette(s)
            rails[0].append(self.point(s, a))
            rails[1].append(self.point(s, b))
        rail_a, rail_b = rails
        if not near(ang(P, rail_a[0]), zang + math.pi / 2):
            rail_a, rail_b = rail_b, rail_a

        def cap(C, start, sweep):
            return [(C[0] + r * math.cos(start + sweep * k / 12), C[1] + r * math.sin(start + sweep * k / 12))
                    for k in range(13)]

        path = cap(P, zang - math.pi / 2, -math.pi)       # from rail_b's foot to rail_a's, away from Z
        path += rail_a
        a1 = ang(Q, rail_a[-1])
        direction = 1.0 if near(a1, zang + math.pi / 2) else -1.0
        path += cap(Q, a1, direction * math.pi)
        path += rail_b[::-1]
        return path

    def visible_runs(self, samples):
        """Split a list of (s, phi) into runs of consecutive visible / hidden points."""
        runs = []
        current = None
        for s, phi in samples:
            vis = self.facing(s, phi) > 0
            p = self.point(s, phi)
            if current is None or current[0] != vis:
                if current is not None:
                    current[1].append(p)       # share the point so that the runs join
                current = (vis, [p])
                runs.append(current)
            else:
                current[1].append(p)
        return runs

    def meridian(self):
        """The belt circle at the top of the arch as visible / hidden runs."""
        samples = [(math.pi / 2, 2 * math.pi * k / 72) for k in range(73)]
        return self.visible_runs(samples)


class SurfaceDrawing:
    """Geometry of the Heegaard diagram on the surface, computed from a HeegaardDiagram."""

    def __init__(self, diagram):
        self.d = diagram
        self._build()

    # ------------------------------------------------------------------ geometry

    def _build(self):
        d = self.d
        self.balls = []          # (centre, radius)
        self.feet = {}           # vertex -> centre
        self.dart_angle = {}     # dart -> angle on the page
        self.arcs = {}           # arc index -> smoothed polyline on the ball
        self.face_points = {}    # face index -> a point inside it (inner faces)
        self.outer_faces = {}    # component -> its outer face
        offset = 0.0
        self.gap = 1.1
        self.r = None
        layouts = []
        for ci, layout in enumerate(d.layouts):
            theta = self._best_rotation(layout)
            centres = {v: _rotate(p, theta) for v, p in layout["centres"].items()}
            edges = {idx: [_rotate(p, theta) for p in pts] for idx, pts in layout["edges"].items()}
            face_points = {fi: _rotate(p, theta) for fi, p in layout["face_points"].items()}
            angles = {e: a + theta for e, a in layout["dart_angle"].items()}
            pts = [p for pl in edges.values() for p in pl] + list(centres.values())
            cx = sum(p[0] for p in pts) / len(pts)
            cy = sum(p[1] for p in pts) / len(pts)
            radius = max(math.hypot(p[0] - cx, p[1] - cy) for p in pts) + 2.4 * layout["rho"]
            layouts.append((centres, edges, face_points, angles, (cx, cy), radius, layout))
        total = 0.0
        for k, (centres, edges, face_points, angles, (cx, cy), radius, layout) in enumerate(layouts):
            shift = (offset + radius - cx, -cy)
            self.balls.append(((offset + radius, 0.0), radius))
            for v, p in centres.items():
                self.feet[v] = (p[0] + shift[0], p[1] + shift[1])
            for fi, p in face_points.items():
                self.face_points[fi] = (p[0] + shift[0], p[1] + shift[1])
            self.outer_faces[k] = layout["outer"]
            self.dart_angle.update(angles)
            r = 0.7 * min(layout["radii"].values())
            self.r = r if self.r is None else min(self.r, r)
            for idx, pl in edges.items():
                self.arcs[idx] = [(p[0] + shift[0], p[1] + shift[1]) for p in pl]
            offset += 2 * radius + self.gap
            total = offset - self.gap
        self.width = total
        # arcs trimmed at the (uniform) foot radius and smoothed
        trimmed = {}
        for idx, pl in self.arcs.items():
            arc = d.arcs[idx]
            va, vb = d.vertex_of[arc["exit"]], d.vertex_of[arc["entry"]]
            pts = list(pl)
            for end, v in ((0, va), (-1, vb)):
                c = self.feet[v]
                q = pts[1 if end == 0 else -2]
                dist = math.hypot(q[0] - c[0], q[1] - c[1])
                pts[end] = (c[0] + self.r * (q[0] - c[0]) / dist, c[1] + self.r * (q[1] - c[1]) / dist)
            trimmed[idx] = pts
        self.arcs = _smooth_without_crossings(trimmed, d._count_crossings)
        # the handles
        self.handles = []
        for i in range(d.g):
            self.handles.append(self._handle(i))
        # tubes joining the balls (L > 1): ball ci to ball ci + 1, attached in the
        # regions that the algorithm's permutation F chose
        self.joins = []
        self.join_feet = []
        if d.L > 1:
            region = {}
            for fi in d.tube_faces:
                region[d.component_of[d.faces[fi][0]]] = fi
            for ci in range(d.L - 1):
                P = self._join_foot(ci, region[ci], +1)
                Q = self._join_foot(ci + 1, region[ci + 1], -1)
                self.join_feet += [P, Q]
                self.joins.append(Tube(P, Q, self.r, 3.0 * self.r))

    def _best_rotation(self, layout):
        """Rotate a component so that its handles are as horizontal as possible (a tube
        whose feet are aligned vertically on the page would be foreshortened)."""
        dirs = []
        for i in range(self.d.g):
            a, b = 2 * i, 2 * i + 1
            if a in layout["centres"] and b in layout["centres"]:
                pa, pb = layout["centres"][a], layout["centres"][b]
                dirs.append(math.atan2(pb[1] - pa[1], pb[0] - pa[0]))
        if not dirs:
            return 0.0
        zang = math.atan2(Z[1], Z[0])
        best = (-1.0, 0.0)
        for k in range(180):
            theta = math.pi * k / 180
            # a tube parallel to Z rises along itself (the arch is foreshortened to a
            # straight tube) and one perpendicular to Z has a degenerate belt ellipse:
            # keep every tube close to 45 degrees from Z
            worst = min(abs(math.sin(2 * (a + theta - zang))) for a in dirs)
            if worst > best[0]:
                best = (worst, theta)
        return best[1]

    def _join_foot(self, ci, face, side):
        """Where a joining tube attaches on ball ci, in the region ``face``; ``side`` is
        +1 for the tube to the next ball (on the right) and -1 for the previous one."""
        centre, radius = self.balls[ci]
        if face != self.outer_faces[ci]:
            return self.face_points[face]
        # the outer face is the rest of the ball: attach near the rim, facing the neighbour
        return (centre[0] + side * 0.8 * radius, centre[1] - 0.45 * radius)

    def _handle(self, i):
        d = self.d
        P, Q = self.feet[2 * i], self.feet[2 * i + 1]
        tube = Tube(P, Q, self.r, max(3.5 * self.r, 0.5 * math.hypot(Q[0] - P[0], Q[1] - P[1]) / 2))
        # strands: one per sheet; the unbarred dart at P, its B-partner at Q
        strands = []
        for e in d.vertex_darts[2 * i]:
            e_bar = d.B(e)
            j = d.dart_info[e][1]
            phi_p = (math.pi - (self.dart_angle[e] - tube.theta_u)) % (2 * math.pi)
            phi_q = (self.dart_angle[e_bar] - tube.theta_u) % (2 * math.pi)
            rel = d.arcs[d.edge_of[e]]["rel"]
            # the letter is a_i if the curve enters at the unbarred disk, a_i^-1 otherwise
            forward = d.arcs[d.edge_of[e]]["entry"] == e
            strands.append({"sheet": j, "phi_p": phi_p, "phi_q": phi_q, "rel": rel, "forward": forward})
        # lift the angles so that the strands do not cross on the tube: both ends are in
        # the same cyclic order, so one consistent choice of twists exists
        strands.sort(key=lambda s: s["phi_p"])
        if strands:
            base_p = strands[0]["phi_p"]
            base_q = strands[0]["phi_q"]
            for s in strands:
                s["lift_p"] = base_p + ((s["phi_p"] - base_p) % (2 * math.pi))
                s["lift_q"] = base_q + ((s["phi_q"] - base_q) % (2 * math.pi))
            # choose the overall twist of the bundle with the smallest maximal twist
            best = None
            for m in range(-2, 3):
                twist = [s["lift_q"] + 2 * math.pi * m - s["lift_p"] for s in strands]
                score = max(abs(t) for t in twist)
                if best is None or score < best[0]:
                    best = (score, m)
            m = best[1]
            for s in strands:
                s["twist"] = s["lift_q"] + 2 * math.pi * m - s["lift_p"]
                # the twist happens in the middle of the arch, where the tube is seen from
                # the side; near the feet the strands climb straight up
                samples = []
                for k in range(SAMPLES + 1):
                    t = math.pi * k / SAMPLES
                    samples.append((t, s["lift_p"] + s["twist"] * (1 - math.cos(t)) / 2))
                s["runs"] = tube.visible_runs(samples)
                s["crossing"] = tube.point(math.pi / 2, s["lift_p"] + s["twist"] / 2)
                s["crossing_visible"] = tube.facing(math.pi / 2, s["lift_p"] + s["twist"] / 2) > 0
        return {"tube": tube, "strands": strands, "generator": d.gens[i][0]}

    # ------------------------------------------------------------------ helpers for renderers

    def label_arcs(self):
        return self.d._label_arcs()

    def bounds(self):
        xs, ys = [], []
        for centre, radius in self.balls:
            xs += [centre[0] - radius, centre[0] + radius]
            ys += [centre[1] - radius, centre[1] + radius]
        for h in self.handles:
            for p in h["tube"].outline():
                xs.append(p[0])
                ys.append(p[1])
        for t in self.joins:
            for p in t.outline():
                xs.append(p[0])
                ys.append(p[1])
        return min(xs), max(xs), min(ys), max(ys)

    # ------------------------------------------------------------------ SVG

    PALETTE = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf", "#8c564b", "#e377c2"]
    ALPHA = "#1a4fb0"

    def to_svg(self, width=760, standalone=True):
        d = self.d
        x0, x1, y0, y1 = self.bounds()
        pad = 0.06 * (x1 - x0)
        x0, x1, y0, y1 = x0 - pad, x1 + pad, y0 - pad, y1 + pad
        scale = width / (x1 - x0)
        height = (y1 - y0) * scale

        def P(p):
            return (scale * (p[0] - x0), scale * (y1 - p[1]))

        def path(points, close=False):
            return "M " + " L ".join("%.2f %.2f" % P(p) for p in points) + (" Z" if close else "")

        out = []
        if standalone:
            out.append('<?xml version="1.0" encoding="UTF-8"?>')
        out.append('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d" '
                   'font-family="Helvetica, Arial, sans-serif">' % (width, height, width, height))
        out.append('<defs><radialGradient id="ball" cx="38%" cy="32%" r="75%">'
                   '<stop offset="0%" stop-color="#ffffff"/><stop offset="100%" stop-color="#d9dde3"/>'
                   '</radialGradient>')
        clips = []
        for k, h in enumerate(self.handles):
            clips.append((k, h["tube"]))
        for k, t in enumerate(self.joins):
            clips.append(("j%d" % k, t))
        for key, tube in clips:
            out.append('<clipPath id="band%s"><path d="%s"/></clipPath>' % (key, path(tube.outline(), True)))
        out.append('</defs>')
        out.append('<rect width="100%" height="100%" fill="white"/>')
        # balls
        for centre, radius in self.balls:
            cx, cy = P(centre)
            out.append('<circle cx="%.2f" cy="%.2f" r="%.2f" fill="url(#ball)" stroke="#333" stroke-width="1.5"/>'
                       % (cx, cy, radius * scale))
        # arcs of the beta-curves on the balls
        label_arcs = self.label_arcs()
        arc_paths = {}
        for idx, pts in self.arcs.items():
            colour = self.PALETTE[d.arcs[idx]["rel"] % len(self.PALETTE)]
            arc_paths[idx] = (path(pts), colour)
            out.append('<path d="%s" fill="none" stroke="%s" stroke-width="2" stroke-linejoin="round" '
                       'stroke-linecap="round"/>' % (arc_paths[idx][0], colour))
            screen = [P(p) for p in pts]
            (ax, ay), (tx, ty) = _polyline_point(screen, 0.5)
            out.append(self._arrow(ax, ay, tx, ty, colour))
            if d.arcs[idx]["closing"]:
                (bx, by), _ = _polyline_point(screen, 0.72)
                out.append('<circle cx="%.2f" cy="%.2f" r="4" fill="%s"/>' % (bx, by, colour))
            if idx in label_arcs:
                (lx, ly), (tx, ty) = _polyline_point(screen, 0.3)
                lx, ly = lx - 13 * ty, ly + 13 * tx
                out.append('<text x="%.2f" y="%.2f" font-size="16" font-style="italic" fill="%s" '
                           'text-anchor="middle" dominant-baseline="middle">&#946;<tspan font-size="11" '
                           'baseline-shift="sub">%d</tspan></text>' % (lx, ly, colour, d.arcs[idx]["rel"] + 1))
        # feet
        for v, c in self.feet.items():
            cx, cy = P(c)
            out.append('<circle cx="%.2f" cy="%.2f" r="%.2f" fill="#f0f0f0" stroke="#333" stroke-width="1.3"/>'
                       % (cx, cy, self.r * scale))
        if self.join_feet:
            for c in self.join_feet:
                cx, cy = P(c)
                out.append('<circle cx="%.2f" cy="%.2f" r="%.2f" fill="#f0f0f0" stroke="#333" stroke-width="1.3"/>'
                           % (cx, cy, self.r * scale))
        # tubes joining balls
        for k, tube in enumerate(self.joins):
            out += self._svg_tube(tube, "j%d" % k, arc_paths, P, path, scale, strands=None, label=None)
        # handles, the lower ones first so that a taller arch hides what it passes over
        order = sorted(range(len(self.handles)), key=lambda k: self.handles[k]["tube"].h)
        for k in order:
            h = self.handles[k]
            out += self._svg_tube(h["tube"], k, arc_paths, P, path, scale, strands=h["strands"],
                                  label=h["generator"])
        out.append("</svg>")
        return "\n".join(out)

    def _svg_tube(self, tube, key, arc_paths, P, path, scale, strands, label):
        out = []
        outline = tube.outline()
        out.append('<path d="%s" fill="#ececec" stroke="none"/>' % path(outline, True))
        # what lies under the tube, dashed
        for idx, (pth, colour) in arc_paths.items():
            out.append('<path d="%s" fill="none" stroke="%s" stroke-width="1.5" stroke-dasharray="4 4" '
                       'clip-path="url(#band%s)" opacity="0.8"/>' % (pth, colour, key))
        for v, c in self.feet.items():
            cx, cy = P(c)
            out.append('<circle cx="%.2f" cy="%.2f" r="%.2f" fill="none" stroke="#333" stroke-width="1" '
                       'stroke-dasharray="4 4" clip-path="url(#band%s)"/>' % (cx, cy, self.r * scale, key))
        # hidden strands
        if strands:
            for s in strands:
                colour = self.PALETTE[s["rel"] % len(self.PALETTE)]
                for vis, pts in s["runs"]:
                    if not vis and len(pts) > 1:
                        out.append('<path d="%s" fill="none" stroke="%s" stroke-width="1.5" '
                                   'stroke-dasharray="4 4" opacity="0.8"/>' % (path(pts), colour))
        # outline of the tube
        out.append('<path d="%s" fill="none" stroke="#333" stroke-width="1.5" stroke-linejoin="round"/>'
                   % path(outline, True))
        # hidden part of the belt circle, visible strands, visible part of the belt circle
        if label is not None:
            for vis, pts in tube.meridian():
                if not vis and len(pts) > 1:
                    out.append('<path d="%s" fill="none" stroke="%s" stroke-width="1.8" stroke-dasharray="5 4"/>'
                               % (path(pts), self.ALPHA))
        if strands:
            for s in strands:
                colour = self.PALETTE[s["rel"] % len(self.PALETTE)]
                for vis, pts in s["runs"]:
                    if vis and len(pts) > 1:
                        out.append('<path d="%s" fill="none" stroke="%s" stroke-width="2" stroke-linejoin="round" '
                                   'stroke-linecap="round"/>' % (path(pts), colour))
        if label is not None:
            for vis, pts in tube.meridian():
                if vis and len(pts) > 1:
                    out.append('<path d="%s" fill="none" stroke="%s" stroke-width="2.2"/>' % (path(pts), self.ALPHA))
            apex = tube.axis(math.pi / 2)
            lx, ly = P((apex[0] + 1.6 * tube.r * Z[0] / 0.5, apex[1] + 1.6 * tube.r * Z[1] / 0.5))
            out.append('<text x="%.2f" y="%.2f" font-size="16" font-style="italic" fill="%s" text-anchor="middle" '
                       'dominant-baseline="middle">&#945;<tspan font-size="11" baseline-shift="sub">%s</tspan>'
                       '</text>' % (lx, ly, self.ALPHA, html.escape(label)))
            # sheet numbers at the crossings with the belt circle
            for s in strands:
                cx, cy = P(s["crossing"])
                colour = self.PALETTE[s["rel"] % len(self.PALETTE)]
                out.append('<circle cx="%.2f" cy="%.2f" r="2.4" fill="%s"%s/>'
                           % (cx, cy, colour, "" if s["crossing_visible"] else ' opacity="0.5"'))
                out.append('<text x="%.2f" y="%.2f" font-size="9" fill="#555" text-anchor="middle" '
                           'dominant-baseline="middle">%d</text>' % (cx + 6, cy - 6, s["sheet"]))
        return out

    @staticmethod
    def _arrow(x, y, tx, ty, colour, size=7):
        nx, ny = -ty, tx
        p1 = (x + size * tx, y + size * ty)
        p2 = (x - size * 0.6 * tx + size * 0.55 * nx, y - size * 0.6 * ty + size * 0.55 * ny)
        p3 = (x - size * 0.6 * tx - size * 0.55 * nx, y - size * 0.6 * ty - size * 0.55 * ny)
        return '<polygon points="%.2f,%.2f %.2f,%.2f %.2f,%.2f" fill="%s"/>' % (p1 + p2 + p3 + (colour,))

    # ------------------------------------------------------------------ TikZ

    def to_tikz(self, width_cm=10.0):
        d = self.d
        colours = ["red", "blue", "green!60!black", "violet", "orange", "cyan!70!black", "brown", "magenta"]
        x0, x1, y0, y1 = self.bounds()
        scale = width_cm / (x1 - x0)

        def P(p):
            return (scale * (p[0] - x0), scale * (p[1] - y0))

        def coords(points):
            return " ".join("(%.3f,%.3f)" % P(p) for p in points)

        out = [
            "%% Heegaard diagram for %s on the genus-%d surface" % (d.presentation_text(), d.g),
            "%% generated by check.py; needs \\usetikzlibrary{arrows.meta, decorations.markings}",
            "\\begin{tikzpicture}[line cap=round, line join=round,",
            "    beta/.style={thick},",
            "    hidden/.style={thin, dashed},",
            "    alpha/.style={thick, blue!70!black},",
            "    mid arrow/.style={postaction={decorate, decoration={markings,",
            "        mark=at position #1 with {\\arrow{Stealth}}}}}]",
        ]
        for centre, radius in self.balls:
            cx, cy = P(centre)
            out.append("  \\shade[inner color=white, outer color=black!18] (%.3f,%.3f) circle (%.3f);" % (cx, cy, radius * scale))
            out.append("  \\draw (%.3f,%.3f) circle (%.3f);" % (cx, cy, radius * scale))
        label_arcs = self.label_arcs()
        for idx, pts in self.arcs.items():
            colour = colours[d.arcs[idx]["rel"] % len(colours)]
            out.append("  \\draw[beta, %s, mid arrow=0.5] plot coordinates {%s};"
                       % (colour, coords(pts)))
            if d.arcs[idx]["closing"]:
                (bx, by), _ = _polyline_point(pts, 0.72)
                out.append("  \\fill[%s] (%.3f,%.3f) circle (0.06);" % ((colour,) + P((bx, by))))
            if idx in label_arcs:
                (lx, ly), (tx, ty) = _polyline_point(pts, 0.3)
                lx, ly = P((lx, ly))
                lx, ly = lx + 0.35 * ty, ly - 0.35 * tx
                out.append("  \\node[%s] at (%.3f,%.3f) {$\\beta_{%d}$};" % (colour, lx, ly, d.arcs[idx]["rel"] + 1))
        feet = list(self.feet.values()) + list(self.join_feet)
        for c in feet:
            out.append("  \\draw[fill=black!5] (%.3f,%.3f) circle (%.3f);" % (P(c) + (self.r * scale,)))
        tubes = [(t, None, None) for t in self.joins]
        tubes += sorted([(h["tube"], h["strands"], h["generator"]) for h in self.handles], key=lambda t: t[0].h)
        for tube, strands, label in tubes:
            outline = tube.outline()
            out.append("  \\fill[black!8] plot coordinates {%s} -- cycle;" % coords(outline))
            out.append("  \\begin{scope}")
            out.append("    \\clip plot coordinates {%s} -- cycle;" % coords(outline))
            for idx, pts in self.arcs.items():
                colour = colours[d.arcs[idx]["rel"] % len(colours)]
                out.append("    \\draw[hidden, %s] plot coordinates {%s};" % (colour, coords(pts)))
            for c in feet:
                out.append("    \\draw[hidden] (%.3f,%.3f) circle (%.3f);" % (P(c) + (self.r * scale,)))
            out.append("  \\end{scope}")
            if strands:
                for s in strands:
                    colour = colours[s["rel"] % len(colours)]
                    for vis, pts in s["runs"]:
                        if not vis and len(pts) > 1:
                            out.append("  \\draw[hidden, %s] plot coordinates {%s};" % (colour, coords(pts)))
            out.append("  \\draw plot coordinates {%s} -- cycle;" % coords(outline))
            if label is not None:
                for vis, pts in tube.meridian():
                    if not vis and len(pts) > 1:
                        out.append("  \\draw[alpha, dashed] plot coordinates {%s};" % coords(pts))
            if strands:
                for s in strands:
                    colour = colours[s["rel"] % len(colours)]
                    for vis, pts in s["runs"]:
                        if vis and len(pts) > 1:
                            out.append("  \\draw[beta, %s] plot coordinates {%s};" % (colour, coords(pts)))
            if label is not None:
                for vis, pts in tube.meridian():
                    if vis and len(pts) > 1:
                        out.append("  \\draw[alpha] plot coordinates {%s};" % coords(pts))
                apex = tube.axis(math.pi / 2)
                lx, ly = P((apex[0] + 1.6 * tube.r * Z[0] / 0.5, apex[1] + 1.6 * tube.r * Z[1] / 0.5))
                out.append("  \\node[blue!70!black] at (%.3f,%.3f) {$\\alpha_{%s}$};" % (lx, ly, label))
                for s in strands:
                    cx, cy = P(s["crossing"])
                    colour = colours[s["rel"] % len(colours)]
                    out.append("  \\fill[%s%s] (%.3f,%.3f) circle (0.04);"
                               % (colour, "" if s["crossing_visible"] else ", opacity=0.5", cx, cy))
                    out.append("  \\node[font=\\tiny, gray] at (%.3f,%.3f) {%d};" % (cx + 0.18, cy + 0.18, s["sheet"]))
        out.append("\\end{tikzpicture}")
        return "\n".join(out) + "\n"
