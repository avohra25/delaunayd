/**
 * Incremental Delaunay triangulation with a compact half-edge structure.
 *
 * Points are inserted in lexicographic order (a sweep): every new point lies
 * strictly outside the current convex hull, so it is connected to the hull
 * edges it can see and the new edges are made Delaunay by Lawson flips.
 * All decisions use the exact predicates from `delaunayd.predicates`, so
 * duplicates, collinear runs and cocircular sets (e.g. integer grids) are
 * handled without epsilons.
 *
 * Storage follows the "flat half-edge" layout: half-edges `3t`, `3t+1`,
 * `3t+2` form triangle `t` in counterclockwise order, `triangles[e]` is the
 * vertex half-edge `e` starts from and `halfedges[e]` is its twin (or
 * `NONE` on the convex hull).
 */
module delaunayd.triangulation;

import std.algorithm : sort;
import std.exception : enforce;

import delaunayd.point;
import delaunayd.predicates;

/// Marker for "no half-edge / no vertex".
enum int NONE = -1;

/// Next half-edge (counterclockwise) within the same triangle.
int nextHalfedge(int e) @safe pure nothrow @nogc { return e % 3 == 2 ? e - 2 : e + 1; }
/// Previous half-edge within the same triangle.
int prevHalfedge(int e) @safe pure nothrow @nogc { return e % 3 == 0 ? e + 2 : e - 1; }

/// Result of point location.
struct Location
{
    enum Kind
    {
        outside,    /// outside the convex hull (`edge` is a hull edge facing the point, if any)
        inTriangle, /// strictly inside `triangle`
        onEdge,     /// on the interior of half-edge `edge` (of `triangle`)
        onVertex,   /// coincides with `vertex`
    }

    Kind kind;
    int triangle = NONE;
    int edge = NONE;
    int vertex = NONE;
}

/// A Delaunay triangulation of a point set.
final class Triangulation
{
    /// The input points, unchanged. Indices below refer to this array.
    Point[] points;
    /// Start vertex of each half-edge; triangle t is half-edges 3t..3t+2 (CCW).
    int[] triangles;
    /// Twin of each half-edge, or NONE for convex hull edges.
    int[] halfedges;
    /// Convex hull vertices in counterclockwise order, starting at the
    /// lexicographically smallest point. Includes vertices lying on hull
    /// edges. For collinear input it is just the two extreme points.
    int[] hull;
    /// rep[i] is the index of the first point equal to points[i]; duplicates
    /// are not vertices of the triangulation.
    int[] rep;
    /// The distinct vertices in lexicographic order.
    int[] vertices;
    /// One outgoing half-edge per vertex (NONE for duplicates or when there
    /// are no triangles).
    int[] vertexEdge;

    private int[] hullNext, hullPrev, hullTri;
    private bool[ulong] fixedEdges; // constrained edges (unordered vertex pairs)

    /// Builds the Delaunay triangulation of `pts`. Throws on NaN/infinite input.
    this(const(Point)[] pts)
    {
        foreach (p; pts)
            enforce(p.isValid, "Triangulation: non-finite coordinate");
        enforce(pts.length < int.max / 8, "Triangulation: too many points");
        points = pts.dup;
        const n = cast(int) pts.length;
        rep = new int[n];
        vertexEdge = new int[n];
        vertexEdge[] = NONE;
        hullNext = new int[n];
        hullPrev = new int[n];
        hullTri = new int[n];

        auto order = new int[n];
        foreach (i, ref v; order)
            v = cast(int) i;
        order.sort!((a, b) => points[a] < points[b] || (points[a] == points[b] && a < b));
        foreach (i; order)
        {
            if (vertices.length && points[vertices[$ - 1]] == points[i])
                rep[i] = vertices[$ - 1];
            else
            {
                vertices ~= i;
                rep[i] = i;
            }
        }
        build();
    }

    /// Number of triangles.
    size_t triangleCount() const @safe pure nothrow @nogc { return triangles.length / 3; }

    /// Vertex indices of triangle t (counterclockwise).
    int[3] triangle(size_t t) const @safe pure nothrow @nogc
    {
        return [triangles[3 * t], triangles[3 * t + 1], triangles[3 * t + 2]];
    }

    /// True when no triangle exists (fewer than 3 distinct points or all collinear).
    bool isDegenerate() const @safe pure nothrow @nogc { return triangles.length == 0; }

    /// Each undirected edge once, as [from, to] vertex pairs. For collinear
    /// input this is the path through the points in order.
    int[2][] edges() const
    {
        int[2][] r;
        if (isDegenerate)
        {
            foreach (i; 1 .. vertices.length)
                r ~= [vertices[i - 1], vertices[i]];
            return r;
        }
        foreach (e; 0 .. cast(int) triangles.length)
            if (halfedges[e] == NONE || e < halfedges[e])
                r ~= [triangles[e], triangles[nextHalfedge(e)]];
        return r;
    }

    /**
     * Outgoing half-edges of vertex v in counterclockwise order. For hull
     * vertices the list starts with the hull edge leaving v; the edge to the
     * previous hull vertex has no outgoing half-edge (use `neighbors`).
     */
    int[] outgoing(int v) const
    {
        int start = vertexEdge[v];
        if (start == NONE)
            return [];
        // rotate clockwise to the hull boundary (if any)
        int e = start;
        while (true)
        {
            const tw = halfedges[e];
            if (tw == NONE)
                break;
            e = nextHalfedge(tw);
            if (e == start)
                break;
        }
        start = e;
        int[] r;
        do
        {
            r ~= e;
            const tw = halfedges[prevHalfedge(e)];
            if (tw == NONE)
                break;
            e = tw;
        }
        while (e != start);
        return r;
    }

    /// Neighbouring vertices of v in counterclockwise order.
    int[] neighbors(int v) const
    {
        if (rep[v] != v)
            return [];
        if (isDegenerate)
        {
            int[] r;
            foreach (i, u; vertices)
                if (u == v)
                {
                    if (i > 0) r ~= vertices[i - 1];
                    if (i + 1 < vertices.length) r ~= vertices[i + 1];
                }
            return r;
        }
        auto out_ = outgoing(v);
        int[] r;
        foreach (e; out_)
            r ~= triangles[nextHalfedge(e)];
        const last = prevHalfedge(out_[$ - 1]);
        if (halfedges[last] == NONE)
            r ~= triangles[last];
        return r;
    }

    /// True if the edge (u, v) exists in the triangulation.
    bool hasEdge(int u, int v) const
    {
        return findHalfedge(u, v) != NONE || findHalfedge(v, u) != NONE
            || (isDegenerate && isChainEdge(u, v));
    }

    /// The half-edge from u to v, or NONE.
    int findHalfedge(int u, int v) const
    {
        if (isDegenerate || vertexEdge[u] == NONE)
            return NONE;
        foreach (e; outgoing(u))
            if (triangles[nextHalfedge(e)] == v)
                return e;
        return NONE;
    }

    private bool isChainEdge(int u, int v) const
    {
        foreach (i; 1 .. vertices.length)
            if ((vertices[i - 1] == u && vertices[i] == v) || (vertices[i - 1] == v && vertices[i] == u))
                return true;
        return false;
    }

    /// Circumcenter of triangle t (floating point, computed relative to a vertex).
    Point circumcenter(size_t t) const
    {
        const a = points[triangles[3 * t]], b = points[triangles[3 * t + 1]],
            c = points[triangles[3 * t + 2]];
        return circumcenterOf(a, b, c);
    }

    /**
     * Locates q: walks from `hint` (a triangle index) using a stochastic
     * visibility walk, which terminates on any triangulation, including
     * constrained ones. Falls back to a linear scan if the walk is cut off.
     */
    Location locate(Point q, int hint = 0) const
    {
        Location loc;
        if (isDegenerate)
        {
            foreach (v; vertices)
                if (points[v] == q)
                {
                    loc.kind = Location.Kind.onVertex;
                    loc.vertex = v;
                    return loc;
                }
            return loc; // outside
        }
        int t = (hint >= 0 && hint < triangleCount) ? hint : 0;
        int entered = NONE;
        uint seed = 0x9E3779B9u ^ cast(uint) t;
        const maxSteps = 4 * triangleCount + 64;
        foreach (step; 0 .. maxSteps)
        {
            seed = seed * 1664525u + 1013904223u;
            const r = (seed >> 16) % 3;
            bool moved = false;
            foreach (k; 0 .. 3)
            {
                const e = 3 * t + cast(int)((r + k) % 3);
                if (e == entered)
                    continue;
                if (orient2d(points[triangles[e]], points[triangles[nextHalfedge(e)]], q) < 0)
                {
                    const tw = halfedges[e];
                    if (tw == NONE)
                    {
                        loc.kind = Location.Kind.outside;
                        loc.edge = e;
                        loc.triangle = t;
                        return loc;
                    }
                    entered = tw;
                    t = tw / 3;
                    moved = true;
                    break;
                }
            }
            if (!moved)
                return classify(q, t);
        }
        // Linear fallback (not expected to be reached).
        foreach (tt; 0 .. cast(int) triangleCount)
        {
            bool inside = true;
            foreach (k; 0 .. 3)
            {
                const e = 3 * tt + k;
                if (orient2d(points[triangles[e]], points[triangles[nextHalfedge(e)]], q) < 0)
                {
                    inside = false;
                    break;
                }
            }
            if (inside)
                return classify(q, tt);
        }
        return loc;
    }

    private Location classify(Point q, int t) const
    {
        Location loc;
        loc.triangle = t;
        int[3] zero;
        int nz = 0;
        foreach (k; 0 .. 3)
        {
            const e = 3 * t + k;
            if (orient2d(points[triangles[e]], points[triangles[nextHalfedge(e)]], q) == 0)
                zero[nz++] = e;
        }
        if (nz == 0)
            loc.kind = Location.Kind.inTriangle;
        else if (nz == 1)
        {
            loc.kind = Location.Kind.onEdge;
            loc.edge = zero[0];
        }
        else
        {
            // two edges through q: q is their shared vertex
            loc.kind = Location.Kind.onVertex;
            const e0 = zero[0], e1 = zero[1];
            loc.vertex = nextHalfedge(e0) == e1 ? triangles[e1] : triangles[e0];
        }
        return loc;
    }

    // ---------------------------------------------------------------- build

    private void build()
    {
        const m = vertices.length;
        if (m < 3)
        {
            hull = vertices.dup;
            return;
        }
        const u = vertices;
        size_t k = 2;
        while (k < m && orient2d(points[u[0]], points[u[1]], points[u[k]]) == 0)
            ++k;
        if (k == m)
        {
            hull = [u[0], u[m - 1]];
            return;
        }
        triangles.reserve(6 * m);
        halfedges.reserve(6 * m);
        const p = u[k];
        const left = orient2d(points[u[0]], points[u[1]], points[p]) > 0;
        int prevLink = NONE;
        int firstTri = NONE, lastTri = NONE;
        foreach (i; 0 .. k - 1)
        {
            int t;
            if (left)
            {
                // (u_i, u_i+1, p): e0 hull edge, e1 = u_i+1 -> p, e2 = p -> u_i
                t = addTriangle(u[i], u[i + 1], p, NONE, NONE, prevLink);
                prevLink = t + 1;
                hullTri[u[i]] = t;
            }
            else
            {
                // (u_i+1, u_i, p): e0 hull edge, e1 = u_i -> p, e2 = p -> u_i+1
                t = addTriangle(u[i + 1], u[i], p, NONE, prevLink, NONE);
                prevLink = t + 2;
                hullTri[u[i + 1]] = t;
            }
            if (i == 0) firstTri = t;
            lastTri = t;
        }
        if (left)
        {
            foreach (i; 0 .. k - 1)
            {
                hullNext[u[i]] = u[i + 1];
                hullPrev[u[i + 1]] = u[i];
            }
            hullNext[u[k - 1]] = p; hullPrev[p] = u[k - 1];
            hullNext[p] = u[0]; hullPrev[u[0]] = p;
            hullTri[u[k - 1]] = lastTri + 1;
            hullTri[p] = firstTri + 2;
        }
        else
        {
            foreach (i; 0 .. k - 1)
            {
                hullNext[u[i + 1]] = u[i];
                hullPrev[u[i]] = u[i + 1];
            }
            hullNext[u[0]] = p; hullPrev[p] = u[0];
            hullNext[p] = u[k - 1]; hullPrev[u[k - 1]] = p;
            hullTri[u[0]] = firstTri + 1;
            hullTri[p] = lastTri + 2;
        }

        int[] stack;
        foreach (i; k + 1 .. m)
            insertOutside(u[i], u[i - 1], stack);

        // hull list from the lexicographically smallest vertex
        int h = u[0];
        do
        {
            hull ~= h;
            h = hullNext[h];
        }
        while (h != u[0]);

        foreach (e; 0 .. cast(int) triangles.length)
        {
            const v = triangles[e];
            if (vertexEdge[v] == NONE || halfedges[e] == NONE)
                vertexEdge[v] = e;
        }
    }

    private int addTriangle(int a, int b, int c, int ta, int tb, int tc)
    {
        const t = cast(int) triangles.length;
        triangles ~= [a, b, c];
        halfedges ~= [NONE, NONE, NONE];
        link(t, ta);
        link(t + 1, tb);
        link(t + 2, tc);
        return t;
    }

    private void link(int a, int b)
    {
        halfedges[a] = b;
        if (b != NONE)
            halfedges[b] = a;
    }

    /// Inserts p, which lies strictly outside the hull; q is the previously
    /// inserted point (always on the hull and visible from p).
    private void insertOutside(int p, int q, ref int[] stack)
    {
        const P = points[p];
        int s = q;
        while (orient2d(points[hullPrev[s]], points[s], P) < 0)
            s = hullPrev[s];
        int t = q;
        while (orient2d(points[t], points[hullNext[t]], P) < 0)
            t = hullNext[t];
        assert(s != t, "sweep invariant violated: no visible hull edge");

        int prevE2 = NONE, firstE1 = NONE;
        stack.length = 0;
        for (int w = s; w != t; w = hullNext[w])
        {
            const wn = hullNext[w];
            const tri = addTriangle(wn, w, p, hullTri[w], prevE2, NONE);
            if (w == s)
                firstE1 = tri + 1;
            prevE2 = tri + 2;
            stack ~= tri;
        }
        hullTri[s] = firstE1;
        hullTri[p] = prevE2;
        hullNext[s] = p; hullPrev[p] = s;
        hullNext[p] = t; hullPrev[t] = p;

        // Lawson legalization of the edges opposite p.
        while (stack.length)
        {
            const a = stack[$ - 1];
            stack.length -= 1;
            if (!isIllegal(a))
                continue;
            const b = halfedges[a];
            flip(a);
            stack ~= a;
            stack ~= nextHalfedge(b);
        }
    }

    /// True if the edge of half-edge a fails the (constrained) Delaunay test.
    package bool isIllegal(int a) const
    {
        const b = halfedges[a];
        if (b == NONE)
            return false;
        const v0 = triangles[a], v1 = triangles[nextHalfedge(a)];
        if (fixedEdges.length && isConstrained(v0, v1))
            return false;
        const pA = triangles[prevHalfedge(a)], pB = triangles[prevHalfedge(b)];
        return incircle(points[v0], points[v1], points[pA], points[pB]) > 0;
    }

    /**
     * Flips the edge of half-edge a (the quad must be strictly convex).
     * Triangle A = (v0, v1, pA) and B = (v1, v0, pB) become
     * (pB, v1, pA) and (pA, v0, pB); half-edge ids are reused.
     */
    package void flip(int a)
    {
        const b = halfedges[a];
        const an = nextHalfedge(a), ap = prevHalfedge(a);
        const bn = nextHalfedge(b), bp = prevHalfedge(b);
        const v0 = triangles[a], v1 = triangles[an];
        const pA = triangles[ap], pB = triangles[bp];
        const ha = halfedges[ap], hb = halfedges[bp];
        triangles[a] = pB;
        triangles[b] = pA;
        link(a, hb);
        link(b, ha);
        link(ap, bp);
        if (hb == NONE)
            hullTri[pB] = a;
        if (ha == NONE)
            hullTri[pA] = b;
        vertexEdge[v0] = bn;
        vertexEdge[v1] = an;
        vertexEdge[pA] = ap;
        vertexEdge[pB] = bp;
    }

    package static ulong edgeKey(int u, int v) @safe pure nothrow @nogc
    {
        if (u > v) { const t = u; u = v; v = t; }
        return (cast(ulong) cast(uint) u << 32) | cast(uint) v;
    }

    /// True if (u, v) has been inserted as a constrained edge.
    bool isConstrained(int u, int v) const
    {
        return (edgeKey(u, v) in fixedEdges) !is null;
    }
}

/// Circumcenter of triangle abc (floating point; translated to a for accuracy).
Point circumcenterOf(Point a, Point b, Point c) @safe pure nothrow @nogc
{
    const bx = b.x - a.x, by = b.y - a.y;
    const cx = c.x - a.x, cy = c.y - a.y;
    const bl = bx * bx + by * by, cl = cx * cx + cy * cy;
    const d = 2 * (bx * cy - by * cx);
    return Point(a.x + (cy * bl - by * cl) / d, a.y + (bx * cl - cx * bl) / d);
}

/// Convenience wrapper: `new Triangulation(pts)`.
Triangulation delaunay(const(Point)[] pts)
{
    return new Triangulation(pts);
}

// ------------------------------------------------------------------ tests

version (unittest)
{
    import std.random : Mt19937, uniform, uniform01;
    import std.math : cos, sin, PI;
    import delaunayd.hull : convexHull, inConvexPolygon;

    /// Structural and Delaunay validation; `brute` adds the O(n*T) empty
    /// circumcircle check against every vertex.
    void validate(const Triangulation tr, bool brute)
    {
        const pts = tr.points;
        const m = tr.vertices.length;
        if (tr.isDegenerate)
        {
            // fewer than 3 distinct points or all collinear
            if (m >= 3)
                foreach (v; tr.vertices)
                    assert(orient2d(pts[tr.vertices[0]], pts[tr.vertices[$ - 1]], pts[v]) == 0);
            return;
        }
        const T = tr.triangleCount;
        const h = tr.hull.length;
        assert(T == 2 * m - 2 - h, "Euler count mismatch");
        auto used = new bool[pts.length];
        size_t hullEdges = 0;
        foreach (e; 0 .. cast(int) tr.triangles.length)
        {
            used[tr.triangles[e]] = true;
            assert(tr.rep[tr.triangles[e]] == tr.triangles[e]);
            const tw = tr.halfedges[e];
            if (tw == NONE)
            {
                ++hullEdges;
                continue;
            }
            assert(tr.halfedges[tw] == e);
            assert(tr.triangles[tw] == tr.triangles[nextHalfedge(e)]);
            assert(tr.triangles[nextHalfedge(tw)] == tr.triangles[e]);
            if (e < tw && !tr.isConstrained(tr.triangles[e], tr.triangles[tw]))
            {
                const v0 = tr.triangles[e], v1 = tr.triangles[tw];
                const pA = tr.triangles[prevHalfedge(e)], pB = tr.triangles[prevHalfedge(tw)];
                assert(incircle(pts[v0], pts[v1], pts[pA], pts[pB]) <= 0, "edge not locally Delaunay");
            }
        }
        assert(hullEdges == h);
        foreach (v; tr.vertices)
            assert(used[v], "vertex missing from triangulation");
        foreach (t; 0 .. T)
        {
            const tri = tr.triangle(t);
            assert(orient2d(pts[tri[0]], pts[tri[1]], pts[tri[2]]) > 0, "triangle not CCW");
        }
        foreach (i; 0 .. h)
            assert(orient2d(pts[tr.hull[i]], pts[tr.hull[(i + 1) % h]], pts[tr.hull[(i + 2) % h]]) >= 0);
        // strict hull corners match the monotone chain hull
        int[] corners;
        foreach (i; 0 .. h)
            if (orient2d(pts[tr.hull[(i + h - 1) % h]], pts[tr.hull[i]], pts[tr.hull[(i + 1) % h]]) > 0)
                corners ~= tr.hull[i];
        auto ref_ = convexHull(pts);
        assert(corners.length == ref_.length);
        foreach (i, c; corners)
            assert(c == ref_[i]);
        // vertex stars are consistent with the edge list
        size_t degreeSum = 0;
        foreach (v; tr.vertices)
            degreeSum += tr.neighbors(v).length;
        assert(degreeSum == 2 * tr.edges().length);
        if (brute)
            foreach (t; 0 .. T)
            {
                const tri = tr.triangle(t);
                foreach (v; tr.vertices)
                    assert(incircle(pts[tri[0]], pts[tri[1]], pts[tri[2]], pts[v]) <= 0,
                            "point inside a circumcircle");
            }
    }

    Point[] randomPoints(ref Mt19937 rng, size_t n, int mode)
    {
        Point[] pts;
        foreach (i; 0 .. n)
        {
            final switch (mode)
            {
            case 0: // uniform
                pts ~= Point(uniform01(rng), uniform01(rng));
                break;
            case 1: // small integer grid: massively cocircular + duplicates
                pts ~= Point(uniform(0, 6, rng), uniform(0, 6, rng));
                break;
            case 2: // near-cocircular: points on a circle (rounded)
                const a = uniform01(rng) * 2 * PI;
                pts ~= Point(cos(a), sin(a));
                break;
            case 3: // collinear runs plus a few strays
                if (uniform(0, 5, rng) == 0)
                    pts ~= Point(uniform01(rng), uniform01(rng));
                else
                {
                    const t = uniform(0, 20, rng);
                    pts ~= Point(t * 0.1, t * 0.1 * 3 + 0.25);
                }
                break;
            case 4: // clustered with large offset (cancellation-heavy)
                pts ~= Point(1e6 + uniform01(rng) * 1e-3, -1e6 + uniform01(rng) * 1e-3);
                break;
            case 5: // duplicates of a few random points
                pts ~= Point(uniform(0, 3, rng) * 0.37, uniform(0, 3, rng) * 0.73);
                break;
            }
        }
        return pts;
    }
}

unittest
{
    // tiny and degenerate inputs
    assert(delaunay([]).isDegenerate);
    assert(delaunay([Point(0, 0)]).hull == [0]);
    auto two = delaunay([Point(0, 0), Point(1, 1), Point(0, 0)]);
    assert(two.isDegenerate && two.hull == [0, 1] && two.rep == [0, 1, 0]);
    auto line = delaunay([Point(2, 2), Point(0, 0), Point(1, 1), Point(3, 3)]);
    assert(line.isDegenerate && line.hull == [1, 3]);
    assert(line.edges().length == 3 && line.hasEdge(1, 2) && !line.hasEdge(1, 3));
    assert(line.neighbors(2) == [1, 0]);

    auto tri = delaunay([Point(0, 0), Point(1, 0), Point(0, 1)]);
    assert(tri.triangleCount == 1 && tri.hull == [0, 1, 2]);
    validate(tri, true);

    // a square has two triangles; its corners are cocircular so either
    // diagonal is Delaunay
    auto sq = delaunay([Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1)]);
    assert(sq.triangleCount == 2);
    validate(sq, true);

    // collinear points followed by one apex: fan
    auto fan = delaunay([Point(0, 0), Point(1, 0), Point(2, 0), Point(3, 0), Point(1.5, -5)]);
    assert(fan.triangleCount == 3);
    validate(fan, true);
}

unittest
{
    // thousands of seeded random and degenerate inputs, brute-force checked
    auto rng = Mt19937(20261111);
    foreach (trial; 0 .. 3000)
    {
        const mode = trial % 6;
        const n = uniform(1, 70, rng);
        auto pts = randomPoints(rng, n, mode);
        validate(delaunay(pts), true);
    }
}

unittest
{
    // full integer grids: every unit square is cocircular
    foreach (w; 2 .. 12)
        foreach (hgt; 1 .. 12)
        {
            Point[] pts;
            foreach (i; 0 .. w)
                foreach (j; 0 .. hgt)
                    pts ~= Point(i, j);
            auto tr = delaunay(pts);
            validate(tr, true);
            if (hgt > 1)
                assert(tr.triangleCount == 2 * (w - 1) * (hgt - 1));
        }
}

unittest
{
    // regular polygon (exactly symmetric points are not exactly cocircular
    // in floating point; the exact predicates decide consistently)
    foreach (n; [5, 12, 64, 257])
    {
        Point[] pts;
        foreach (i; 0 .. n)
            pts ~= Point(cos(2 * PI * i / n), sin(2 * PI * i / n));
        pts ~= Point(0, 0);
        validate(delaunay(pts), true);
    }
    // extreme scales
    auto rng = Mt19937(5);
    foreach (scale; [1e-30, 1e-8, 1.0, 1e8, 1e30])
    {
        Point[] pts;
        foreach (_; 0 .. 200)
            pts ~= Point(uniform01(rng) * scale, uniform01(rng) * scale);
        validate(delaunay(pts), true);
    }
}

unittest
{
    // large input: structural + local Delaunay check
    auto rng = Mt19937(42);
    auto pts = randomPoints(rng, 20_000, 0);
    auto tr = delaunay(pts);
    validate(tr, false);
    // grid with exact duplicates appended
    Point[] g;
    foreach (i; 0 .. 100)
        foreach (j; 0 .. 100)
            g ~= Point(i * 0.5, j * 0.25);
    g ~= g[0 .. 500];
    validate(delaunay(g), false);
}

unittest
{
    // point location against a brute-force classification
    auto rng = Mt19937(77);
    foreach (trial; 0 .. 300)
    {
        const mode = trial % 2 == 0 ? 0 : 1;
        auto pts = randomPoints(rng, uniform(3, 80, rng), mode);
        auto tr = delaunay(pts);
        if (tr.isDegenerate)
            continue;
        Point[] hullPoly;
        foreach (v; tr.hull)
            hullPoly ~= pts[v];
        Point[] queries;
        foreach (_; 0 .. 50)
            queries ~= mode == 0 ? Point(uniform(-0.2, 1.2, rng), uniform(-0.2, 1.2, rng))
                : Point(uniform(-1, 13, rng) * 0.5, uniform(-1, 13, rng) * 0.5);
        foreach (v; tr.vertices)
            queries ~= pts[v];
        foreach (q; queries)
        {
            const loc = tr.locate(q, uniform(0, cast(int) tr.triangleCount, rng));
            final switch (loc.kind)
            {
            case Location.Kind.outside:
                assert(!inConvexPolygon(hullPoly, q));
                break;
            case Location.Kind.onVertex:
                assert(pts[loc.vertex] == q);
                break;
            case Location.Kind.onEdge:
                const a = pts[tr.triangles[loc.edge]], b = pts[tr.triangles[nextHalfedge(loc.edge)]];
                assert(orient2d(a, b, q) == 0 && q != a && q != b);
                assert(dot(q - a, b - a) > 0 && dot(q - b, a - b) > 0);
                break;
            case Location.Kind.inTriangle:
                const t = tr.triangle(loc.triangle);
                foreach (k; 0 .. 3)
                    assert(orient2d(pts[t[k]], pts[t[(k + 1) % 3]], q) > 0);
                break;
            }
        }
    }
}
