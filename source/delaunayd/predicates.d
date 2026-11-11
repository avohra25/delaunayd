/**
 * Robust geometric predicates: orient2d and incircle.
 *
 * Both follow Shewchuk's adaptive approach: a cheap floating-point
 * evaluation guarded by a forward error bound, and only when the result is
 * too close to zero to trust, progressively more exact evaluations, ending
 * in an exact computation with expansion arithmetic. The returned value has
 * the exact sign of the true determinant; its magnitude is approximate.
 */
module delaunayd.predicates;

import delaunayd.expansion;
import delaunayd.point;

@safe pure nothrow
{

private enum double resultErrBound = (3.0 + 8.0 * epsilon) * epsilon;
private enum double ccwErrBoundA = (3.0 + 16.0 * epsilon) * epsilon;
private enum double ccwErrBoundB = (2.0 + 12.0 * epsilon) * epsilon;
private enum double ccwErrBoundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon;
private enum double iccErrBoundA = (10.0 + 96.0 * epsilon) * epsilon;

/**
 * Positive if a, b, c are in counterclockwise order, negative if clockwise,
 * zero if collinear. Approximates twice the signed triangle area.
 */
double orient2d(Point a, Point b, Point c)
{
    const double detLeft = (a.x - c.x) * (b.y - c.y);
    const double detRight = (a.y - c.y) * (b.x - c.x);
    const double det = detLeft - detRight;
    double detSum;
    if (detLeft > 0.0)
    {
        if (detRight <= 0.0)
            return det;
        detSum = detLeft + detRight;
    }
    else if (detLeft < 0.0)
    {
        if (detRight >= 0.0)
            return det;
        detSum = -detLeft - detRight;
    }
    else
        return det;
    const double errBound = ccwErrBoundA * detSum;
    if (det >= errBound || -det >= errBound)
        return det;
    return orient2dAdapt(a, b, c, detSum);
}

private double orient2dAdapt(Point a, Point b, Point c, double detSum)
{
    const double acx = a.x - c.x, bcx = b.x - c.x;
    const double acy = a.y - c.y, bcy = b.y - c.y;

    // Stage B: exact determinant of the rounded differences.
    auto bExp = difference(productExpansion(acx, bcy), productExpansion(acy, bcx));
    double det = estimate(bExp);
    double errBound = ccwErrBoundB * detSum;
    if (det >= errBound || -det >= errBound)
        return det;

    double t;
    double acxTail, bcxTail, acyTail, bcyTail;
    twoDiff(a.x, c.x, t, acxTail);
    twoDiff(b.x, c.x, t, bcxTail);
    twoDiff(a.y, c.y, t, acyTail);
    twoDiff(b.y, c.y, t, bcyTail);
    if (acxTail == 0.0 && acyTail == 0.0 && bcxTail == 0.0 && bcyTail == 0.0)
        return det; // the differences were exact, so stage B was exact

    // Stage C: first-order correction from the difference tails.
    errBound = ccwErrBoundC * detSum + resultErrBound * (det >= 0 ? det : -det);
    det += (acx * bcyTail + bcy * acxTail) - (acy * bcxTail + bcx * acyTail);
    if (det >= errBound || -det >= errBound)
        return det;

    // Stage D: fully exact.
    return orient2dExactValue(a, b, c);
}

/// The determinant computed exactly; returns its rounded value (exact sign).
double orient2dExactValue(Point a, Point b, Point c)
{
    auto acx = diffExpansion(a.x, c.x), bcx = diffExpansion(b.x, c.x);
    auto acy = diffExpansion(a.y, c.y), bcy = diffExpansion(b.y, c.y);
    auto d = difference(product(acx, bcy), product(acy, bcx));
    return signedEstimate(d);
}

/**
 * Positive if d lies inside the circle through a, b, c (which must be in
 * counterclockwise order), negative if outside, zero if cocircular.
 * The sign is reversed if a, b, c are clockwise.
 */
double incircle(Point a, Point b, Point c, Point d)
{
    const double adx = a.x - d.x, bdx = b.x - d.x, cdx = c.x - d.x;
    const double ady = a.y - d.y, bdy = b.y - d.y, cdy = c.y - d.y;

    const double bdxcdy = bdx * cdy, cdxbdy = cdx * bdy;
    const double alift = adx * adx + ady * ady;
    const double cdxady = cdx * ady, adxcdy = adx * cdy;
    const double blift = bdx * bdx + bdy * bdy;
    const double adxbdy = adx * bdy, bdxady = bdx * ady;
    const double clift = cdx * cdx + cdy * cdy;

    const double det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy)
        + clift * (adxbdy - bdxady);
    const double permanent = (abs(bdxcdy) + abs(cdxbdy)) * alift
        + (abs(cdxady) + abs(adxcdy)) * blift + (abs(adxbdy) + abs(bdxady)) * clift;
    const double errBound = iccErrBoundA * permanent;
    if (det > errBound || -det > errBound)
        return det;
    return incircleExactValue(a, b, c, d);
}

/// The incircle determinant computed exactly (exact sign).
double incircleExactValue(Point a, Point b, Point c, Point d)
{
    auto adx = diffExpansion(a.x, d.x), bdx = diffExpansion(b.x, d.x);
    auto cdx = diffExpansion(c.x, d.x);
    auto ady = diffExpansion(a.y, d.y), bdy = diffExpansion(b.y, d.y);
    auto cdy = diffExpansion(c.y, d.y);
    auto alift = sum(product(adx, adx), product(ady, ady));
    auto blift = sum(product(bdx, bdx), product(bdy, bdy));
    auto clift = sum(product(cdx, cdx), product(cdy, cdy));
    auto bc = difference(product(bdx, cdy), product(cdx, bdy));
    auto ca = difference(product(cdx, ady), product(adx, cdy));
    auto ab = difference(product(adx, bdy), product(bdx, ady));
    auto det = sum(sum(product(alift, bc), product(blift, ca)), product(clift, ab));
    return signedEstimate(det);
}

/// Exact sign of orient2d as -1, 0 or +1.
int orientSign(Point a, Point b, Point c)
{
    const v = orient2d(a, b, c);
    return v > 0 ? 1 : v < 0 ? -1 : 0;
}

/// Exact sign of incircle as -1, 0 or +1.
int incircleSign(Point a, Point b, Point c, Point d)
{
    const v = incircle(a, b, c, d);
    return v > 0 ? 1 : v < 0 ? -1 : 0;
}

/// The plain floating-point orient2d (non-robust; for comparison/testing).
double orient2dFast(Point a, Point b, Point c) @nogc
{
    return (a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x);
}

private double abs(double v) @nogc { return v < 0 ? -v : v; }

/// Estimate whose sign is guaranteed to match the expansion's exact sign.
private double signedEstimate(const(double)[] e)
{
    const s = signOf(e);
    if (s == 0)
        return 0.0;
    const v = estimate(e);
    // The largest component dominates, so the estimate already has the right
    // sign; guard anyway against pathological underflow.
    if ((v > 0) == (s > 0) && v != 0)
        return v;
    return s > 0 ? double.min_normal : -double.min_normal;
}

} // attributes

version (unittest)
{
    import std.bigint : BigInt;

    // Points in [0.25, 1) or so are integer multiples of 2^-54.
    private BigInt big(double v)
    {
        const scaled = v * 0x1p54;
        assert(scaled == cast(double) cast(long) scaled);
        return BigInt(cast(long) scaled);
    }

    private int bigSign(BigInt v)
    {
        return v > 0 ? 1 : v < 0 ? -1 : 0;
    }

    private int orientRef(Point a, Point b, Point c)
    {
        auto d = (big(a.x) - big(c.x)) * (big(b.y) - big(c.y))
            - (big(a.y) - big(c.y)) * (big(b.x) - big(c.x));
        return bigSign(d);
    }

    private int incircleRef(Point a, Point b, Point c, Point d)
    {
        auto adx = big(a.x) - big(d.x), ady = big(a.y) - big(d.y);
        auto bdx = big(b.x) - big(d.x), bdy = big(b.y) - big(d.y);
        auto cdx = big(c.x) - big(d.x), cdy = big(c.y) - big(d.y);
        auto det = (adx * adx + ady * ady) * (bdx * cdy - cdx * bdy)
            + (bdx * bdx + bdy * bdy) * (cdx * ady - adx * cdy)
            + (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady);
        return bigSign(det);
    }
}

unittest
{
    assert(orient2d(Point(0, 0), Point(1, 0), Point(0, 1)) > 0);
    assert(orient2d(Point(0, 0), Point(0, 1), Point(1, 0)) < 0);
    assert(orient2d(Point(0, 0), Point(1, 1), Point(2, 2)) == 0);
    assert(incircle(Point(0, 0), Point(1, 0), Point(0, 1), Point(0.5, 0.5)) > 0);
    assert(incircle(Point(0, 0), Point(1, 0), Point(0, 1), Point(1, 1)) == 0);
    assert(incircle(Point(0, 0), Point(1, 0), Point(0, 1), Point(2, 2)) < 0);
}

unittest
{
    // Kettner et al. "Classroom examples of robustness problems": the naive
    // predicate gets the orientation of points near the line y = x wrong.
    import std.math : nextUp;

    int naiveWrong = 0;
    const p = Point(0.5, 0.5);
    const q = Point(12, 12), r = Point(24, 24);
    foreach (i; 0 .. 64)
        foreach (j; 0 .. 64)
        {
            Point a = p;
            foreach (_; 0 .. i) a.x = nextUp(a.x);
            foreach (_; 0 .. j) a.y = nextUp(a.y);
            // Reference with scaled integers (q, r have larger exponents but
            // are exact multiples of 2^-54 as well).
            const ref_ = orientRef(a, q, r);
            assert(orientSign(a, q, r) == ref_);
            assert(orientSign(q, r, a) == ref_);
            assert(orientSign(r, a, q) == ref_);
            const fast = orient2dFast(a, q, r);
            if ((fast > 0 ? 1 : fast < 0 ? -1 : 0) != ref_)
                ++naiveWrong;
        }
    assert(naiveWrong > 0, "expected the naive predicate to fail somewhere");
}

unittest
{
    // Randomized near-degenerate inputs checked against exact BigInt values.
    import std.random : Mt19937, uniform;

    auto rng = Mt19937(2026);
    double coord(ref Mt19937 g)
    {
        // a multiple of 2^-20 in [0.25, 0.75), perturbed by a few ulps
        double v = 0.25 + uniform(0, 1 << 19, g) * 0x1p-20;
        const k = uniform(-3, 4, g);
        return v + k * 0x1p-54;
    }

    foreach (_; 0 .. 20_000)
    {
        // Nearly collinear: c = a + t (b - a) rounded, then perturbed.
        const a = Point(coord(rng), coord(rng));
        const b = Point(coord(rng), coord(rng));
        const t = uniform(0, 1 << 10, rng) * 0x1p-10;
        auto c = Point(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y));
        // snap c onto the 2^-54 grid so the reference stays exact
        c.x = cast(double) cast(long)(c.x * 0x1p54) * 0x1p-54;
        c.y = cast(double) cast(long)(c.y * 0x1p54) * 0x1p-54;
        assert(orientSign(a, b, c) == orientRef(a, b, c));
    }

    foreach (_; 0 .. 5000)
    {
        // Nearly cocircular: perturbations of points on an axis-aligned square
        // (whose corners are exactly cocircular).
        const cx = 0.5, cy = 0.5;
        const h = uniform(1, 1 << 18, rng) * 0x1p-20;
        Point[4] sq = [Point(cx - h, cy - h), Point(cx + h, cy - h),
            Point(cx + h, cy + h), Point(cx - h, cy + h)];
        foreach (ref pt; sq)
        {
            pt.x += uniform(-2, 3, rng) * 0x1p-54;
            pt.y += uniform(-2, 3, rng) * 0x1p-54;
        }
        assert(incircleSign(sq[0], sq[1], sq[2], sq[3]) == incircleRef(sq[0], sq[1], sq[2], sq[3]));
        assert(incircleSign(sq[1], sq[2], sq[3], sq[0]) == incircleRef(sq[1], sq[2], sq[3], sq[0]));
    }
}

unittest
{
    // Exact paths agree with the filtered paths on clear-cut inputs.
    import std.random : Mt19937, uniform01;

    auto rng = Mt19937(7);
    foreach (_; 0 .. 2000)
    {
        Point[4] p;
        foreach (ref q; p)
            q = Point(uniform01(rng), uniform01(rng));
        const o = orient2d(p[0], p[1], p[2]);
        const oe = orient2dExactValue(p[0], p[1], p[2]);
        assert((o > 0) == (oe > 0) && (o < 0) == (oe < 0));
        const ic = incircle(p[0], p[1], p[2], p[3]);
        const ice = incircleExactValue(p[0], p[1], p[2], p[3]);
        assert((ic > 0) == (ice > 0) && (ic < 0) == (ice < 0));
    }
}
