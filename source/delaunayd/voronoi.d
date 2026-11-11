/**
 * Voronoi diagram as the dual of a Delaunay triangulation.
 *
 * Voronoi vertices are triangle circumcenters; each Delaunay edge yields a
 * Voronoi edge (a segment between the circumcenters of its two triangles,
 * or a ray for convex hull edges). Cells clipped to a box are computed by
 * intersecting the box with the bisector half-planes of the site's Delaunay
 * neighbours, which also covers collinear inputs (where the diagram is a set
 * of parallel lines and there are no triangles).
 */
module delaunayd.voronoi;

import std.exception : enforce;

import delaunayd.point;
import delaunayd.predicates;
import delaunayd.triangulation;

/// Axis-aligned rectangle.
struct Box
{
    double minX, minY, maxX, maxY;

    /// Bounding box of points, grown by `margin` times its larger side (at least 1 unit if empty).
    static Box around(const(Point)[] pts, double margin = 0.1) @safe pure nothrow
    {
        if (pts.length == 0)
            return Box(-1, -1, 1, 1);
        Box b = Box(pts[0].x, pts[0].y, pts[0].x, pts[0].y);
        foreach (p; pts)
        {
            if (p.x < b.minX) b.minX = p.x;
            if (p.y < b.minY) b.minY = p.y;
            if (p.x > b.maxX) b.maxX = p.x;
            if (p.y > b.maxY) b.maxY = p.y;
        }
        double s = b.maxX - b.minX;
        if (b.maxY - b.minY > s) s = b.maxY - b.minY;
        if (s == 0) s = 1;
        const m = s * margin;
        return Box(b.minX - m, b.minY - m, b.maxX + m, b.maxY + m);
    }

    double area() const @safe pure nothrow @nogc { return (maxX - minX) * (maxY - minY); }

    bool contains(Point p) const @safe pure nothrow @nogc
    {
        return p.x >= minX && p.x <= maxX && p.y >= minY && p.y <= maxY;
    }
}

/// One Voronoi edge, dual to the Delaunay edge (siteA, siteB).
struct VoronoiEdge
{
    enum Kind { segment, ray, line }

    Kind kind;
    int siteA, siteB;
    /// segment: from -> to; ray: starts at `from` going along `dir`;
    /// line: passes through `from` with direction `dir`.
    Point from, to, dir;
}

/// The Voronoi diagram of a triangulation's (distinct) points.
final class Voronoi
{
    const Triangulation tri;
    /// Circumcenter of every Delaunay triangle (the Voronoi vertices).
    Point[] centers;

    this(const Triangulation t)
    {
        enforce(t.constraintCount == 0,
                "Voronoi: requires an unconstrained Delaunay triangulation");
        tri = t;
        centers = new Point[t.triangleCount];
        foreach (i; 0 .. t.triangleCount)
            centers[i] = t.circumcenter(i);
    }

    /// Voronoi edges, one per Delaunay edge.
    VoronoiEdge[] edges() const
    {
        VoronoiEdge[] r;
        const P = tri.points;
        if (tri.isDegenerate)
        {
            foreach (e; tri.edges())
            {
                const a = P[e[0]], b = P[e[1]];
                VoronoiEdge ve;
                ve.kind = VoronoiEdge.Kind.line;
                ve.siteA = e[0];
                ve.siteB = e[1];
                ve.from = Point((a.x + b.x) / 2, (a.y + b.y) / 2);
                ve.dir = Point(-(b.y - a.y), b.x - a.x);
                r ~= ve;
            }
            return r;
        }
        foreach (e; 0 .. cast(int) tri.triangles.length)
        {
            const tw = tri.halfedges[e];
            if (tw != NONE && tw < e)
                continue;
            VoronoiEdge ve;
            ve.siteA = tri.triangles[e];
            ve.siteB = tri.triangles[nextHalfedge(e)];
            ve.from = centers[e / 3];
            if (tw != NONE)
            {
                ve.kind = VoronoiEdge.Kind.segment;
                ve.to = centers[tw / 3];
            }
            else
            {
                // hull edge a->b has the interior on its left: go right
                const a = P[ve.siteA], b = P[ve.siteB];
                ve.kind = VoronoiEdge.Kind.ray;
                ve.dir = Point(b.y - a.y, -(b.x - a.x));
            }
            r ~= ve;
        }
        return r;
    }

    /**
     * Indices into `centers` of the Voronoi vertices of site v's cell in
     * counterclockwise order. `unbounded` is set for hull sites, whose cell
     * continues to infinity after the last vertex and returns before the
     * first one.
     */
    int[] cellVertices(int v, out bool unbounded) const
    {
        unbounded = true;
        if (tri.isDegenerate || tri.rep[v] != v)
            return [];
        auto out_ = tri.outgoing(v);
        int[] r;
        foreach (e; out_)
            r ~= e / 3;
        unbounded = tri.halfedges[out_[0]] == NONE;
        return r;
    }

    /// Site v's cell intersected with `box` as a CCW polygon (may be empty).
    Point[] clippedCell(int v, Box box) const
    {
        if (tri.rep[v] != v)
            return [];
        Point[] poly = [Point(box.minX, box.minY), Point(box.maxX, box.minY),
            Point(box.maxX, box.maxY), Point(box.minX, box.maxY)];
        const s = tri.points[v];
        foreach (u; tri.neighbors(v))
        {
            const q = tri.points[u];
            poly = clipHalfPlane(poly, Point((s.x + q.x) / 2, (s.y + q.y) / 2), q - s);
            if (poly.length == 0)
                break;
        }
        return poly;
    }

    /// All clipped cells, indexed by point (empty for duplicates).
    Point[][] clippedCells(Box box) const
    {
        auto r = new Point[][tri.points.length];
        foreach (v; tri.vertices)
            r[v] = clippedCell(v, box);
        return r;
    }
}

/// Convenience: Voronoi diagram of a point set.
Voronoi voronoi(const(Point)[] pts)
{
    return new Voronoi(new Triangulation(pts));
}

/// Keeps the part of convex polygon `poly` where dot(p - m, n) <= 0.
Point[] clipHalfPlane(const(Point)[] poly, Point m, Point n) @safe pure nothrow
{
    Point[] r;
    if (poly.length == 0)
        return r;
    foreach (i; 0 .. poly.length)
    {
        const p = poly[i], q = poly[(i + 1) % poly.length];
        const dp = dot(p - m, n), dq = dot(q - m, n);
        if (dp <= 0)
            r ~= p;
        if ((dp < 0 && dq > 0) || (dp > 0 && dq < 0))
        {
            const t = dp / (dp - dq);
            r ~= Point(p.x + t * (q.x - p.x), p.y + t * (q.y - p.y));
        }
    }
    return r;
}

version (unittest)
{
    import std.math : abs, sqrt, cos, sin, PI;
    import std.random : Mt19937, uniform, uniform01;
    import delaunayd.hull : polygonArea, inConvexPolygon;

    void checkVoronoi(const(Point)[] pts, ref Mt19937 rng)
    {
        auto tr = delaunay(pts);
        auto vd = new Voronoi(tr);
        // circumcenters are equidistant from their triangle's corners
        foreach (t; 0 .. tr.triangleCount)
        {
            const tt = tr.triangle(t);
            const c = vd.centers[t];
            const r0 = dist2(c, pts[tt[0]]);
            if (r0 > 1e12)
                continue; // near-degenerate sliver: center far away, skip
            foreach (k; 1 .. 3)
                assert(abs(dist2(c, pts[tt[k]]) - r0) <= 1e-6 * (r0 + 1e-300));
        }
        assert(vd.edges().length == tr.edges().length);
        const box = Box.around(pts, 0.25);
        auto cells = vd.clippedCells(box);
        double total = 0;
        foreach (v; tr.vertices)
        {
            const a = polygonArea(cells[v]);
            assert(a >= 0);
            total += a;
            assert(inConvexPolygon(cells[v], pts[v]) || cells[v].length < 3);
        }
        assert(abs(total - box.area) <= 1e-9 * box.area);
        // random queries: the cell that contains q belongs to a nearest site
        foreach (_; 0 .. 100)
        {
            const q = Point(uniform(box.minX, box.maxX, rng), uniform(box.minY, box.maxY, rng));
            double best = double.infinity, second = double.infinity;
            foreach (v; tr.vertices)
            {
                const d = dist2(q, pts[v]);
                if (d < best) { second = best; best = d; }
                else if (d < second) second = d;
            }
            if (second - best <= 1e-9 * second)
                continue; // near a cell boundary
            foreach (v; tr.vertices)
                if (cells[v].length >= 3 && inConvexPolygon(cells[v], q))
                    assert(dist2(q, pts[v]) == best);
        }
        // bounded cells built from circumcenters match the clipped cells
        const big = Box(box.minX - 1e3 * (box.maxX - box.minX), box.minY - 1e3 * (box.maxY - box.minY),
                box.maxX + 1e3 * (box.maxX - box.minX), box.maxY + 1e3 * (box.maxY - box.minY));
        foreach (v; tr.vertices)
        {
            bool unbounded;
            auto cv = vd.cellVertices(v, unbounded);
            if (unbounded || cv.length == 0)
                continue;
            Point[] poly;
            foreach (c; cv)
                poly ~= vd.centers[c];
            const a1 = polygonArea(poly), a2 = polygonArea(vd.clippedCell(v, big));
            assert(abs(a1 - a2) <= 1e-7 * (a2 + 1e-12));
        }
    }
}

unittest
{
    // Two points: two half-planes split by the bisector.
    auto vd = voronoi([Point(0, 0), Point(2, 0)]);
    assert(vd.centers.length == 0);
    auto e = vd.edges();
    assert(e.length == 1 && e[0].kind == VoronoiEdge.Kind.line && e[0].from == Point(1, 0));
    const box = Box(-1, -1, 3, 1);
    assert(abs(polygonArea(vd.clippedCell(0, box)) - 4) < 1e-12);
    assert(abs(polygonArea(vd.clippedCell(1, box)) - 4) < 1e-12);

    // Square plus center: the center's cell is a diamond.
    auto sq = voronoi([Point(-1, -1), Point(1, -1), Point(1, 1), Point(-1, 1), Point(0, 0)]);
    bool unb;
    auto cv = sq.cellVertices(4, unb);
    assert(!unb && cv.length == 4);
    Point[] diamond;
    foreach (c; cv)
        diamond ~= sq.centers[c];
    assert(abs(polygonArea(diamond) - 2) < 1e-12);
    auto corner = sq.cellVertices(0, unb);
    assert(unb && corner.length == 2);
    size_t rays = 0;
    foreach (ed; sq.edges())
        if (ed.kind == VoronoiEdge.Kind.ray)
            ++rays;
    assert(rays == 4);
}

unittest
{
    auto rng = Mt19937(808);
    foreach (trial; 0 .. 200)
    {
        Point[] pts;
        const n = uniform(2, 60, rng);
        const mode = trial % 4;
        foreach (_; 0 .. n)
        {
            if (mode == 0)
                pts ~= Point(uniform01(rng), uniform01(rng));
            else if (mode == 1)
                pts ~= Point(uniform(0, 6, rng), uniform(0, 6, rng)); // cocircular grid
            else if (mode == 2)
            {
                const t = uniform(0, 30, rng);
                pts ~= Point(t * 0.5, t * 0.25); // collinear
            }
            else
            {
                const a = uniform01(rng) * 2 * PI;
                pts ~= Point(cos(a), sin(a)); // near-cocircular
            }
        }
        checkVoronoi(pts, rng);
    }
}

unittest
{
    // Voronoi needs an unconstrained triangulation.
    import delaunayd.constrained : insertConstraint;

    auto tr = delaunay([Point(0, 0), Point(1, 0), Point(0, 1), Point(1, 1)]);
    insertConstraint(tr, 0, 3);
    bool threw = false;
    try
        new Voronoi(tr);
    catch (Exception)
        threw = true;
    assert(threw);
}
