#!/usr/bin/env python3
"""Check whether the presentation complex of a balanced group presentation is the spine
of a closed orientable 3-manifold (Neuwirth's algorithm), and optionally draw the
corresponding Heegaard diagram.

    uv run check.py BBa Ab              # <a, b | b^-2 a, a^-1 b>
    uv run check.py BBa Ab --heegaard   # ... and open the Heegaard diagram in the browser

Relators are words in the generators; a capital letter is the inverse of the lowercase
generator (A = a^-1).  The generators are the letters that occur, in alphabetical order;
a balanced presentation needs as many relators as generators.

Exit status: 0 if the presentation complex is a spine, 1 if it is not, 2 for bad input.
"""

import argparse
import os
import sys
import time
import webbrowser

import neuwirth as nw
import neuwirthcheck


def parse_presentation(words, generators=None):
    """Turn the relator words into the (gen, rel) input of neuwirthcheck."""
    if not words:
        raise ValueError("no relators given")
    for w in words:
        if not w.isalpha() or not w.isascii():
            raise ValueError("relator %r: only the letters a-z and A-Z are allowed" % w)
    letters = sorted({c.lower() for w in words for c in w})
    if generators is None:
        gens = letters
    else:
        gens = list(generators.lower())
        if not generators.isalpha() or not generators.isascii() or len(set(gens)) != len(gens):
            raise ValueError("--generators needs distinct letters, e.g. 'abc'")
        missing = [c for c in letters if c not in gens]
        if missing:
            raise ValueError("the relators use the generator(s) %s, which are not in --generators %r"
                             % (", ".join(missing), generators))
    if len(words) != len(gens):
        raise ValueError("a balanced presentation needs as many relators as generators: "
                         "got %d relator(s) for the generators %s" % (len(words), ", ".join(gens)))
    unused = [g for g in gens if g not in letters]
    if unused:
        raise ValueError("the generator(s) %s do not occur in any relator; Neuwirth's algorithm needs "
                         "every generator to occur" % ", ".join(unused))
    return [[g, g.upper()] for g in gens], list(words)


def presentation_text(gens, rels):
    def word(w):
        out = []
        k = 0
        while k < len(w):
            letter = w[k]
            j = k
            while j < len(w) and w[j] == letter:
                j += 1
            count = j - k
            exponent = count if letter.islower() else -count
            out.append(letter.lower() if exponent == 1 else "%s^%d" % (letter.lower(), exponent))
            k = j
        return " ".join(out)
    return "<%s | %s>" % (", ".join(g[0] for g in gens), ", ".join(word(w) for w in rels))


def output_name(rels, extension):
    return "heegaard-%s.%s" % ("-".join(rels), extension)


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="check.py",
        description="Decide with Neuwirth's algorithm whether the presentation complex of a balanced "
                    "group presentation is the spine of a closed orientable 3-manifold.",
        epilog="Example: `uv run check.py BBa Ab` checks <a, b | b^-2 a, a^-1 b>.  Exit status 0 means "
               "spine, 1 means not a spine, 2 means bad input.",
    )
    parser.add_argument("relators", nargs="+", metavar="RELATOR",
                        help="relator words; capital letters are inverses (A = a^-1)")
    parser.add_argument("-g", "--generators", metavar="LETTERS",
                        help="the generators, e.g. 'abc' (default: the letters that occur)")
    parser.add_argument("--heegaard", action="store_true",
                        help="draw the Heegaard diagram of the spine found on the genus-g surface and open "
                             "it in the browser (written to heegaard-<relators>.html in the current directory; "
                             "for genus at most 3 the surface is shown in space and can be turned around, "
                             "which needs three.js from a CDN)")
    parser.add_argument("--svg", metavar="FILE", help="also write the Heegaard diagram as SVG")
    parser.add_argument("--tikz", metavar="FILE", help="also write the Heegaard diagram as a TikZ picture")
    parser.add_argument("--html", metavar="FILE", help="write the Heegaard diagram page to this file")
    parser.add_argument("--planar", action="store_true",
                        help="with --svg/--tikz: the planar diagram (surface cut open along the alpha-curves) "
                             "instead of the drawing on the surface")
    parser.add_argument("--no-open", action="store_true", help="with --heegaard: do not open the browser")
    parser.add_argument("-q", "--quiet", action="store_true", help="print only the verdict")
    args = parser.parse_args(argv)

    try:
        gens, rels = parse_presentation(args.relators, args.generators)
    except ValueError as exc:
        parser.error(str(exc))

    print("Presentation %s" % presentation_text(gens, rels))
    A, B, C = neuwirthcheck.neuwirth_data(gens, rels)
    prev = nw.prevalences(gens, rels)
    configurations = 1
    for p in prev:
        for j in range(1, p):
            configurations *= j
    if not args.quiet:
        print("Sheets per generator: %s; configurations to try: %d" % (
            ", ".join("%s: %d" % (g[0], p) for g, p in zip(gens, prev)), configurations))
        sys.stdout.flush()

    t0 = time.time()
    found = neuwirthcheck.find_2_spine(A, B, C)
    elapsed = time.time() - t0

    if found is None:
        print("K_P is NOT the spine of a closed orientable 3-manifold (no configuration works; %.1fs)" % elapsed)
        if args.heegaard or args.svg or args.tikz or args.html:
            print("No Heegaard diagram to draw.")
        return 1

    p, Ap = found
    L = neuwirthcheck.count_orbits(Ap, C)
    print("K_P IS the spine of a closed orientable 3-manifold (%.1fs)" % elapsed)
    if not args.quiet:
        identity = list(p) == list(range(1, 2 * sum(prev) + 1))
        print("Witnessing configuration:")
        print("  reordering p = %s" % ("identity" if identity else p.cycle_string()))
        print("  A = %s" % Ap.cycle_string())
        print("  B = %s" % B.cycle_string())
        print("  C = %s" % C.cycle_string(singletons=True))
        print("  L = |<A, C>| = %d%s" % (L, "" if L == 1 else "  (%d spheres joined by tubes)" % L))

    if args.heegaard or args.svg or args.tikz or args.html:
        import heegaard
        t0 = time.time()
        diagram = heegaard.HeegaardDiagram(gens, rels, p, Ap, B, C)
        if not args.quiet:
            print("Heegaard diagram laid out in %.1fs (%d arc%s, %d region%s%s)" % (
                time.time() - t0, len(diagram.arcs), "" if len(diagram.arcs) == 1 else "s",
                len(diagram.faces), "" if len(diagram.faces) == 1 else "s",
                "" if diagram.crossings() == 0 else "; WARNING: the drawing has crossings"))
        written = []
        if args.planar:
            renderer = diagram
        else:
            from surface import SurfaceDrawing
            renderer = SurfaceDrawing(diagram)
        if args.svg:
            with open(args.svg, "w", encoding="utf-8") as fh:
                fh.write(renderer.to_svg())
            written.append(args.svg)
        if args.tikz:
            with open(args.tikz, "w", encoding="utf-8") as fh:
                fh.write(renderer.to_tikz())
            written.append(args.tikz)
        html_path = args.html
        if args.heegaard and html_path is None:
            html_path = output_name(rels, "html")
        if html_path:
            with open(html_path, "w", encoding="utf-8") as fh:
                fh.write(diagram.to_html())
            written.append(html_path)
            if not args.quiet and diagram.pretzel() is not None:
                print("The page shows the surface in space (three.js from cdn.jsdelivr.net; drag to turn it).")
        for path in written:
            print("Wrote %s" % path)
        if args.heegaard and not args.no_open:
            url = "file://" + os.path.abspath(html_path)
            if not webbrowser.open(url):
                print("Could not open a browser; open %s yourself." % html_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
