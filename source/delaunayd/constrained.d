/**
 * Constrained Delaunay triangulation: forcing segments into a triangulation.
 *
 * `insertConstraint` traces the segment through the triangulation, removes
 * the crossing edges with Sloan's flip algorithm (flip any crossing edge
 * whose quadrilateral is strictly convex, requeue the others) and then
 * restores the constrained Delaunay property with Lawson flips that never
 * touch constrained edges. Segments passing exactly through other vertices
 * are split there, so collinear inputs need no special handling.
 */
module delaunayd.constrained;

import std.exception : enforce;
import std.format : format;

import delaunayd.point;
import delaunayd.predicates;
import delaunayd.triangulation;

/// Thrown when a constraint would cross an existing constraint.
class ConstraintError : Exception
{
    this(string msg, string file = __FILE__, size_t line = __LINE__) @safe pure nothrow
    {
        super(msg, file, line);
    }
}

/**
 * Forces the segment between points u and v (input indices; duplicates are
 * mapped to their representative) to be a union of triangulation edges, all
 * marked as constrained. Throws ConstraintError if it properly crosses an
 * existing constrained edge; in that case the triangulation is unchanged.
 */
void insertConstraint(Triangulation tr, int u, int v)
{
    enforce(u >= 0 && u < tr.points.length && v >= 0 && v < tr.points.length,
            "insertConstraint: vertex index out of range");
    u = tr.rep[u];
    v = tr.rep[v];
    if (u == v)
        return;
    if (tr.isDegenerate)
    {
        // all vertices lie on one line (or there are only two): mark the
        // chain pieces between u and v
        size_t iu = size_t.max, iv = size_t.max;
        foreach (i, w; tr.vertices)
        {
            if (w == u) iu = i;
            if (w == v) iv = i;
        }
        if (iu > iv) { const t = iu; iu = iv; iv = t; }
        foreach (i; iu .. iv)
            tr.fixedEdges[Triangulation.edgeKey(tr.vertices[i], tr.vertices[i + 1])] = true;
        return;
    }

    // Pass 1 (read-only): split at collinear vertices and make sure no
    // piece crosses an existing constraint, so a failing insertion leaves
    // the triangulation untouched.
    int[2][] pieces;
    int[2][] work = [[u, v]];
    while (work.length)
    {
        const a = work[$ - 1][0], b = work[$ - 1][1];
        work.length -= 1;
        if (a == b)
            continue;
        if (tr.hasEdge(a, b))
        {
            pieces ~= [a, b];
            continue;
        }
        int splitAt;
        trace(tr, a, b, splitAt); // throws ConstraintError on a crossing
        if (splitAt != NONE)
        {
            work ~= [splitAt, b];
            work ~= [a, splitAt];
            continue;
        }
        pieces ~= [a, b];
    }
    // Pass 2: force each piece in.
    foreach (pc; pieces)
    {
        const a = pc[0], b = pc[1];
        if (!tr.hasEdge(a, b))
        {
            int splitAt;
            auto crossed = trace(tr, a, b, splitAt);
            assert(splitAt == NONE);
            removeCrossings(tr, a, b, crossed);
        }
        tr.fixedEdges[Triangulation.edgeKey(a, b)] = true;
    }
}

/// Inserts every segment of `segs` (pairs of point indices).
void insertConstraints(Triangulation tr, const(int[2])[] segs)
{
    foreach (s; segs)
        insertConstraint(tr, s[0], s[1]);
}

private bool onRay(const Triangulation tr, int a, int b, int w)
{
    const P = tr.points;
    return orientSign(P[a], P[b], P[w]) == 0 && dot(P[w] - P[a], P[b] - P[a]) > 0;
}

/**
 * Edges crossed by the open segment a-b, in order from a, each as
 * [rightVertex, leftVertex]. If a vertex lies exactly on the segment,
 * returns nothing and sets splitAt to the first such vertex.
 */
private int[2][] trace(Triangulation tr, int a, int b, out int splitAt)
{
    const P = tr.points;
    splitAt = NONE;
    int h = NONE;
    int x, y;
    foreach (e; tr.outgoing(a))
    {
        const ex = tr.triangles[nextHalfedge(e)], ey = tr.triangles[prevHalfedge(e)];
        if (onRay(tr, a, b, ex)) { splitAt = ex; return null; }
        if (onRay(tr, a, b, ey)) { splitAt = ey; return null; }
        if (orientSign(P[a], P[b], P[ex]) < 0 && orientSign(P[a], P[b], P[ey]) > 0)
        {
            h = nextHalfedge(e);
            x = ex;
            y = ey;
            break;
        }
    }
    enforce(h != NONE, "insertConstraint: could not start tracing the segment");
    int[2][] crossed;
    while (true)
    {
        if (tr.isConstrained(x, y))
            throw new ConstraintError(format("constraint %s-%s crosses constrained edge %s-%s", a, b, x, y));
        crossed ~= [x, y];
        const tw = tr.halfedges[h];
        enforce(tw != NONE, "insertConstraint: segment left the triangulation");
        const w = tr.triangles[prevHalfedge(tw)];
        if (w == b)
            break;
        const ow = orientSign(P[a], P[b], P[w]);
        if (ow == 0)
        {
            splitAt = w;
            return null;
        }
        if (ow < 0)
        {
            h = prevHalfedge(tw); // w -> y
            x = w;
        }
        else
        {
            h = nextHalfedge(tw); // x -> w
            y = w;
        }
    }
    return crossed;
}

/// True if the open segments p-q and a-b cross at a single interior point.
private bool crossesProperly(const Triangulation tr, int p, int q, int a, int b)
{
    const P = tr.points;
    const o1 = orientSign(P[a], P[b], P[p]), o2 = orientSign(P[a], P[b], P[q]);
    if (o1 == 0 || o2 == 0 || o1 == o2)
        return false;
    const o3 = orientSign(P[p], P[q], P[a]), o4 = orientSign(P[p], P[q], P[b]);
    return o3 != 0 && o4 != 0 && o3 != o4;
}

private void removeCrossings(Triangulation tr, int a, int b, int[2][] queue)
{
    const P = tr.points;
    int[] touched;
    size_t head = 0;
    size_t sinceProgress = 0;
    while (head < queue.length)
    {
        const x = queue[head][0], y = queue[head][1];
        ++head;
        int e = tr.findHalfedge(x, y);
        if (e == NONE)
            e = tr.findHalfedge(y, x);
        enforce(e != NONE, "insertConstraint: lost a crossing edge");
        const f = tr.halfedges[e];
        const v0 = tr.triangles[e], v1 = tr.triangles[f];
        const pA = tr.triangles[prevHalfedge(e)], pB = tr.triangles[prevHalfedge(f)];
        const s0 = orientSign(P[pA], P[pB], P[v0]), s1 = orientSign(P[pA], P[pB], P[v1]);
        if (s0 == 0 || s1 == 0 || s0 == s1)
        {
            // not strictly convex: try again later
            queue ~= [x, y];
            enforce(++sinceProgress <= queue.length - head + 1,
                    "insertConstraint: no flippable crossing edge (invalid triangulation?)");
            continue;
        }
        sinceProgress = 0;
        tr.flip(e);
        const e0 = e - e % 3, f0 = f - f % 3;
        touched ~= [e0, e0 + 1, e0 + 2, f0, f0 + 1, f0 + 2];
        if (crossesProperly(tr, pA, pB, a, b))
            queue ~= [pA, pB];
    }
    tr.fixedEdges[Triangulation.edgeKey(a, b)] = true;
    legalize(tr, touched);
}

/// Lawson flips until every unconstrained edge reachable from `stack` is
/// locally Delaunay.
package void legalize(Triangulation tr, int[] stack)
{
    while (stack.length)
    {
        const e = stack[$ - 1];
        stack.length -= 1;
        if (!tr.isIllegal(e))
            continue;
        const f = tr.halfedges[e];
        tr.flip(e);
        stack ~= [e, nextHalfedge(e), f, nextHalfedge(f)];
    }
}

/// Builds the constrained Delaunay triangulation of points and segments.
Triangulation constrainedDelaunay(const(Point)[] pts, const(int[2])[] segs)
{
    auto tr = new Triangulation(pts);
    insertConstraints(tr, segs);
    return tr;
}

version (unittest)
{
    import std.random : Mt19937, uniform, uniform01;
    import std.algorithm : sort;

    // Every vertex strictly inside a triangle's circumcircle must be hidden
    // from the triangle's interior by some constrained segment.
    void bruteCdtCheck(const Triangulation tr, const(int[2])[] segs)
    {
        const P = tr.points;
        foreach (t; 0 .. tr.triangleCount)
        {
            const tri = tr.triangle(t);
            const c = Point((P[tri[0]].x + P[tri[1]].x + P[tri[2]].x) / 3,
                    (P[tri[0]].y + P[tri[1]].y + P[tri[2]].y) / 3);
            bool centroidInside = true;
            foreach (k; 0 .. 3)
                if (orient2d(P[tri[k]], P[tri[(k + 1) % 3]], c) <= 0)
                    centroidInside = false;
            if (!centroidInside)
                continue; // sliver: rounding moved the centroid out
            foreach (v; tr.vertices)
            {
                if (incircle(P[tri[0]], P[tri[1]], P[tri[2]], P[v]) <= 0)
                    continue;
                bool blocked = false;
                foreach (s; segs)
                {
                    const a = P[s[0]], b = P[s[1]];
                    const o1 = orientSign(a, b, c), o2 = orientSign(a, b, P[v]);
                    const o3 = orientSign(c, P[v], a), o4 = orientSign(c, P[v], b);
                    if (o1 * o2 < 0 && o3 * o4 <= 0)
                    {
                        blocked = true;
                        break;
                    }
                }
                assert(blocked, "visible vertex inside circumcircle");
            }
        }
    }

    // The vertices lying on segment a-b (inclusive), sorted along it.
    int[] piecesOf(const Triangulation tr, int a, int b)
    {
        const P = tr.points;
        int[] on;
        foreach (v; tr.vertices)
            if (orientSign(P[a], P[b], P[v]) == 0 && dot(P[v] - P[a], P[b] - P[a]) >= 0
                    && dot(P[v] - P[b], P[a] - P[b]) >= 0)
                on ~= v;
        on.sort!((p, q) => dot(P[p] - P[a], P[b] - P[a]) < dot(P[q] - P[a], P[b] - P[a]));
        return on;
    }

    bool onSegment(const(Point)[] P, int[2] s, Point w)
    {
        return orientSign(P[s[0]], P[s[1]], w) == 0 && dot(w - P[s[0]], P[s[1]] - P[s[0]]) >= 0
            && dot(w - P[s[1]], P[s[0]] - P[s[1]]) >= 0;
    }

    // True if s and t cross at a single point that is not an input vertex
    // (crossings through a vertex are legal: both get split there).
    bool segmentsCross(const(Point)[] P, int[2] s, int[2] t)
    {
        const o1 = orientSign(P[s[0]], P[s[1]], P[t[0]]), o2 = orientSign(P[s[0]], P[s[1]], P[t[1]]);
        const o3 = orientSign(P[t[0]], P[t[1]], P[s[0]]), o4 = orientSign(P[t[0]], P[t[1]], P[s[1]]);
        if (!(o1 * o2 < 0 && o3 * o4 < 0))
            return false;
        foreach (w; P)
            if (onSegment(P, s, w) && onSegment(P, t, w))
                return false;
        return true;
    }
}

unittest
{
    // A square with a diagonal forced the "wrong" way round.
    auto pts = [Point(0, 0), Point(4, 0), Point(4, 1), Point(0, 1), Point(2, 0.5)];
    auto tr = delaunay(pts);
    assert(!tr.hasEdge(0, 2));
    insertConstraint(tr, 0, 2); // passes through vertex 4 at (2, 0.5)
    assert(tr.hasEdge(0, 4) && tr.hasEdge(4, 2));
    assert(tr.isConstrained(0, 4) && tr.isConstrained(4, 2));
    assert(tr.constraintCount == 2);
    validate(tr, false);

    // A long thin constraint across a dense row of points.
    Point[] row;
    foreach (i; 0 .. 20)
    {
        row ~= Point(i, 1);
        row ~= Point(i + 0.5, -1);
    }
    row ~= Point(-1, 0);
    row ~= Point(21, 0);
    auto t2 = delaunay(row);
    insertConstraint(t2, cast(int) row.length - 2, cast(int) row.length - 1);
    assert(t2.hasEdge(cast(int) row.length - 2, cast(int) row.length - 1));
    validate(t2, false);
    bruteCdtCheck(t2, t2.constrainedEdges());

    // crossing an existing constraint is rejected
    const before = t2.constrainedEdges();
    bool threw = false;
    try
        insertConstraint(t2, 0, 1);
    catch (ConstraintError)
        threw = true;
    assert(threw);
    assert(t2.constrainedEdges() == before); // failed insertions change nothing

    // collinear-only input
    auto line = delaunay([Point(0, 0), Point(1, 1), Point(3, 3), Point(2, 2)]);
    insertConstraint(line, 0, 2);
    assert(line.constraintCount == 3);
}

unittest
{
    // Random point sets with random non-crossing segments.
    auto rng = Mt19937(1111);
    foreach (trial; 0 .. 600)
    {
        const grid = trial % 2 == 1;
        const n = uniform(4, 60, rng);
        Point[] pts;
        foreach (_; 0 .. n)
            pts ~= grid ? Point(uniform(0, 7, rng), uniform(0, 7, rng))
                : Point(uniform01(rng), uniform01(rng));
        auto tr = delaunay(pts);
        int[2][] accepted;
        foreach (_; 0 .. uniform(1, 12, rng))
        {
            int[2] s = [uniform(0, cast(int) n, rng), uniform(0, cast(int) n, rng)];
            if (pts[s[0]] == pts[s[1]])
                continue;
            bool bad = false;
            foreach (t; accepted)
                if (segmentsCross(pts, s, t))
                {
                    bad = true;
                    break;
                }
            if (bad)
            {
                bool threw = false;
                try
                    insertConstraint(tr, s[0], s[1]);
                catch (ConstraintError)
                    threw = true;
                if (threw)
                    continue;
                assert(false, "crossing constraint accepted");
            }
            insertConstraint(tr, s[0], s[1]);
            accepted ~= s;
        }
        validate(tr, false);
        foreach (s; accepted)
        {
            auto on = piecesOf(tr, tr.rep[s[0]], tr.rep[s[1]]);
            foreach (i; 1 .. on.length)
            {
                assert(tr.hasEdge(on[i - 1], on[i]));
                assert(tr.isConstrained(on[i - 1], on[i]));
            }
        }
        if (!tr.isDegenerate)
            bruteCdtCheck(tr, tr.constrainedEdges());
    }
}

unittest
{
    // A polygon outline (star) forced into a random cloud.
    import std.math : cos, sin, PI;

    auto rng = Mt19937(31);
    Point[] pts;
    const k = 24;
    foreach (i; 0 .. k)
    {
        const r = i % 2 ? 0.2 : 0.45;
        pts ~= Point(0.5 + r * cos(2 * PI * i / k), 0.5 + r * sin(2 * PI * i / k));
    }
    foreach (_; 0 .. 400)
        pts ~= Point(uniform01(rng), uniform01(rng));
    int[2][] segs;
    foreach (i; 0 .. k)
        segs ~= [i, (i + 1) % k];
    auto tr = constrainedDelaunay(pts, segs);
    validate(tr, false);
    foreach (s; segs)
        assert(tr.hasEdge(s[0], s[1]) || piecesOf(tr, s[0], s[1]).length > 2);
    bruteCdtCheck(tr, tr.constrainedEdges());
}
