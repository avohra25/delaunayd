/**
 * Convex hull by Andrew's monotone chain, using the robust orient2d
 * predicate so collinear and nearly collinear inputs are handled exactly.
 */
module delaunayd.hull;

import std.algorithm : sort;
import std.exception : enforce;

import delaunayd.point;
import delaunayd.predicates;

/**
 * Indices of the convex hull vertices of `pts` in counterclockwise order,
 * starting at the lexicographically smallest point (smallest x, then y).
 *
 * Only strict corners are reported: points lying on a hull edge are
 * dropped, and duplicates are reported once (lowest index). Degenerate
 * inputs give one index (all points equal) or two (all collinear).
 */
size_t[] convexHull(const(Point)[] pts)
{
    foreach (p; pts)
        enforce(p.isValid, "convexHull: non-finite coordinate");
    if (pts.length == 0)
        return [];
    auto idx = new size_t[pts.length];
    foreach (i, ref v; idx)
        v = i;
    // Stable lexicographic order; ties broken by index so the first duplicate wins.
    idx.sort!((a, b) => pts[a] < pts[b] || (pts[a] == pts[b] && a < b));
    size_t[] uniq;
    foreach (i; idx)
        if (uniq.length == 0 || pts[uniq[$ - 1]] != pts[i])
            uniq ~= i;
    if (uniq.length <= 2)
        return uniq;

    auto h = new size_t[2 * uniq.length];
    size_t k = 0;
    // lower hull
    foreach (i; uniq)
    {
        while (k >= 2 && orient2d(pts[h[k - 2]], pts[h[k - 1]], pts[i]) <= 0)
            --k;
        h[k++] = i;
    }
    // upper hull
    const lowerLen = k + 1;
    foreach_reverse (i; uniq[0 .. $ - 1])
    {
        while (k >= lowerLen && orient2d(pts[h[k - 2]], pts[h[k - 1]], pts[i]) <= 0)
            --k;
        h[k++] = i;
    }
    // the last point equals the first
    auto result = h[0 .. k - 1].dup;
    if (result.length == 1) // all collinear: lower chain collapsed to [first, last]
        return [uniq[0], uniq[$ - 1]];
    return result;
}

/// Signed area of a simple polygon (positive when counterclockwise).
double polygonArea(const(Point)[] poly) @safe pure nothrow @nogc
{
    double a = 0;
    foreach (i; 0 .. poly.length)
    {
        const p = poly[i], q = poly[(i + 1) % poly.length];
        a += p.x * q.y - q.x * p.y;
    }
    return a / 2;
}

/// True if q is inside or on the boundary of the convex CCW polygon (exact).
bool inConvexPolygon(const(Point)[] poly, Point q)
{
    if (poly.length == 0)
        return false;
    if (poly.length == 1)
        return poly[0] == q;
    foreach (i; 0 .. poly.length)
        if (orient2d(poly[i], poly[(i + 1) % poly.length], q) < 0)
            return false;
    if (poly.length == 2)
    {
        // on the segment: collinear and within the bounding box
        const a = poly[0], b = poly[1];
        return orient2d(a, b, q) == 0 && dot(q - a, b - a) >= 0 && dot(q - b, a - b) >= 0;
    }
    return true;
}

version (unittest)
{
    import std.random : Mt19937, uniform;
    import std.algorithm : canFind;

    // O(n^3) reference: a directed edge (a,b) is on the hull iff no point is
    // strictly to its right and every point collinear with it lies within [a,b].
    bool[size_t] bruteHullVertices(const(Point)[] pts)
    {
        bool[size_t] verts;
        size_t[] uniq;
        foreach (i, p; pts)
        {
            bool seen = false;
            foreach (j; uniq)
                if (pts[j] == p) { seen = true; break; }
            if (!seen) uniq ~= i;
        }
        if (uniq.length <= 2)
        {
            foreach (i; uniq) verts[i] = true;
            return verts;
        }
        foreach (a; uniq)
            foreach (b; uniq)
            {
                if (a == b) continue;
                bool ok = true;
                foreach (c; uniq)
                {
                    const o = orient2d(pts[a], pts[b], pts[c]);
                    if (o < 0) { ok = false; break; }
                    if (o == 0 && c != a && c != b)
                    {
                        const ab = pts[b] - pts[a];
                        const t = dot(pts[c] - pts[a], ab);
                        if (t < 0 || t > dot(ab, ab)) { ok = false; break; }
                    }
                }
                if (ok) { verts[a] = true; verts[b] = true; }
            }
        return verts;
    }

    void checkHull(const(Point)[] pts)
    {
        auto h = convexHull(pts);
        auto brute = bruteHullVertices(pts);
        assert(h.length == brute.length);
        foreach (i; h)
            assert(i in brute);
        if (h.length >= 3)
        {
            foreach (i; 0 .. h.length)
                assert(orient2d(pts[h[i]], pts[h[(i + 1) % $]], pts[h[(i + 2) % $]]) > 0);
            Point[] poly;
            foreach (i; h) poly ~= pts[i];
            foreach (p; pts)
                assert(inConvexPolygon(poly, p));
        }
    }
}

unittest
{
    assert(convexHull([]) == []);
    assert(convexHull([Point(1, 1)]) == [0]);
    assert(convexHull([Point(1, 1), Point(1, 1)]) == [0]);
    assert(convexHull([Point(0, 0), Point(2, 2), Point(1, 1)]) == [0, 1]);
    // square with interior, edge-midpoint and duplicate points
    auto sq = [Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1), Point(0.5, 0.5),
        Point(0.5, 0), Point(1, 1)];
    assert(convexHull(sq) == [0, 1, 2, 3]);
    assert(polygonArea([Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1)]) == 1);
}

unittest
{
    auto rng = Mt19937(99);
    // random points, small integer grids (lots of collinearity), duplicates
    foreach (trial; 0 .. 1500)
    {
        const n = uniform(1, 40, rng);
        Point[] pts;
        const mode = trial % 3;
        foreach (_; 0 .. n)
        {
            if (mode == 0)
                pts ~= Point(uniform(-1.0, 1.0, rng), uniform(-1.0, 1.0, rng));
            else if (mode == 1)
                pts ~= Point(uniform(0, 5, rng), uniform(0, 5, rng));
            else
            {
                const t = uniform(0, 10, rng);
                pts ~= Point(t * 0.1, t * 0.3); // nearly/exactly collinear after rounding
            }
        }
        checkHull(pts);
    }
}

unittest
{
    // Points on a circle: every point is a hull vertex.
    import std.math : cos, sin, PI;

    Point[] pts;
    foreach (i; 0 .. 360)
        pts ~= Point(cos(i * PI / 180), sin(i * PI / 180));
    auto h = convexHull(pts);
    assert(h.length == 360);
    assert(h[0] == 180); // (-1, 0) is lexicographically smallest
}
