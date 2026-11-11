/**
 * Exact floating-point expansion arithmetic (Shewchuk, "Adaptive Precision
 * Floating-Point Arithmetic and Fast Robust Geometric Predicates", 1997).
 *
 * An expansion is a sum of non-overlapping doubles stored in increasing order
 * of magnitude; the exact value is the (unrounded) sum of its components.
 * All routines here assume IEEE 754 binary64 with round-to-nearest-even and
 * no extended-precision intermediates or fused multiply-add contraction
 * (the dub build disables FP contraction for LDC).
 */
module delaunayd.expansion;

@safe pure nothrow
{

/// Half an ulp of 1.0: 2^-53.
enum double epsilon = 0x1p-53;
/// 2^ceil(53/2) + 1, used to split a double into two 26-bit halves.
enum double splitter = 0x1p27 + 1.0;

/// x + y == a + b exactly, x = fl(a + b).
void twoSum(double a, double b, out double x, out double y) @nogc
{
    x = a + b;
    const double bv = x - a;
    const double av = x - bv;
    const double br = b - bv;
    const double ar = a - av;
    y = ar + br;
}

/// Like twoSum but requires |a| >= |b| (or a == 0).
void fastTwoSum(double a, double b, out double x, out double y) @nogc
{
    x = a + b;
    const double bv = x - a;
    y = b - bv;
}

/// x + y == a - b exactly, x = fl(a - b).
void twoDiff(double a, double b, out double x, out double y) @nogc
{
    x = a - b;
    const double bv = a - x;
    const double av = x + bv;
    const double br = bv - b;
    const double ar = a - av;
    y = ar + br;
}

/// Splits a into hi + lo where each half has at most 26 significant bits.
void split(double a, out double hi, out double lo) @nogc
{
    const double c = splitter * a;
    const double abig = c - a;
    hi = c - abig;
    lo = a - hi;
}

/// x + y == a * b exactly, x = fl(a * b).
void twoProduct(double a, double b, out double x, out double y) @nogc
{
    x = a * b;
    double ahi, alo, bhi, blo;
    split(a, ahi, alo);
    split(b, bhi, blo);
    const double err1 = x - ahi * bhi;
    const double err2 = err1 - alo * bhi;
    const double err3 = err2 - ahi * blo;
    y = alo * blo - err3;
}

private void twoProductPresplit(double a, double b, double bhi, double blo,
        out double x, out double y) @nogc
{
    x = a * b;
    double ahi, alo;
    split(a, ahi, alo);
    const double err1 = x - ahi * bhi;
    const double err2 = err1 - alo * bhi;
    const double err3 = err2 - ahi * blo;
    y = alo * blo - err3;
}

/// The exact difference a - b as an expansion (zero components removed).
double[] diffExpansion(double a, double b)
{
    double x, y;
    twoDiff(a, b, x, y);
    if (y == 0.0)
        return [x];
    return [y, x];
}

/// The exact product a * b as an expansion.
double[] productExpansion(double a, double b)
{
    double x, y;
    twoProduct(a, b, x, y);
    if (y == 0.0)
        return [x];
    return [y, x];
}

/**
 * h = e + f, exactly, with zero components eliminated
 * (fast_expansion_sum_zeroelim). Never returns an empty expansion.
 */
double[] sum(const(double)[] e, const(double)[] f)
{
    if (e.length == 0)
        return f.length ? f.dup : [0.0];
    if (f.length == 0)
        return e.dup;
    auto h = new double[e.length + f.length];
    size_t hi = 0, ei = 0, fi = 0;
    double enow = e[0], fnow = f[0];
    double q, qnew, hh;
    double nextE() { ++ei; return ei < e.length ? e[ei] : 0.0; }
    double nextF() { ++fi; return fi < f.length ? f[fi] : 0.0; }

    if ((fnow > enow) == (fnow > -enow)) { q = enow; enow = nextE(); }
    else { q = fnow; fnow = nextF(); }
    if (ei < e.length && fi < f.length)
    {
        if ((fnow > enow) == (fnow > -enow)) { fastTwoSum(enow, q, qnew, hh); enow = nextE(); }
        else { fastTwoSum(fnow, q, qnew, hh); fnow = nextF(); }
        q = qnew;
        if (hh != 0.0) h[hi++] = hh;
        while (ei < e.length && fi < f.length)
        {
            if ((fnow > enow) == (fnow > -enow)) { twoSum(q, enow, qnew, hh); enow = nextE(); }
            else { twoSum(q, fnow, qnew, hh); fnow = nextF(); }
            q = qnew;
            if (hh != 0.0) h[hi++] = hh;
        }
    }
    while (ei < e.length)
    {
        twoSum(q, enow, qnew, hh); enow = nextE();
        q = qnew;
        if (hh != 0.0) h[hi++] = hh;
    }
    while (fi < f.length)
    {
        twoSum(q, fnow, qnew, hh); fnow = nextF();
        q = qnew;
        if (hh != 0.0) h[hi++] = hh;
    }
    if (q != 0.0 || hi == 0)
        h[hi++] = q;
    return h[0 .. hi];
}

/// -e.
double[] negate(const(double)[] e)
{
    auto r = new double[e.length];
    foreach (i, c; e)
        r[i] = -c;
    return r;
}

/// e - f, exactly.
double[] difference(const(double)[] e, const(double)[] f)
{
    return sum(e, negate(f));
}

/// h = b * e, exactly (scale_expansion_zeroelim).
double[] scale(const(double)[] e, double b)
{
    if (e.length == 0)
        return [0.0];
    auto h = new double[2 * e.length];
    size_t hi = 0;
    double bhi, blo, q, hh, p1, p0, s;
    split(b, bhi, blo);
    twoProductPresplit(e[0], b, bhi, blo, q, hh);
    if (hh != 0.0) h[hi++] = hh;
    foreach (i; 1 .. e.length)
    {
        twoProductPresplit(e[i], b, bhi, blo, p1, p0);
        twoSum(q, p0, s, hh);
        if (hh != 0.0) h[hi++] = hh;
        fastTwoSum(p1, s, q, hh);
        if (hh != 0.0) h[hi++] = hh;
    }
    if (q != 0.0 || hi == 0)
        h[hi++] = q;
    return h[0 .. hi];
}

/// e * f, exactly (sum of scaled copies).
double[] product(const(double)[] e, const(double)[] f)
{
    double[] acc = [0.0];
    foreach (c; f)
        if (c != 0.0)
            acc = sum(acc, scale(e, c));
    return acc;
}

/// Approximate value (sum of components, smallest first).
double estimate(const(double)[] e) @nogc
{
    double s = 0.0;
    foreach (c; e)
        s += c;
    return s;
}

/// Exact sign (-1, 0, +1): the largest component determines it.
int signOf(const(double)[] e) @nogc
{
    foreach_reverse (c; e)
    {
        if (c > 0.0) return 1;
        if (c < 0.0) return -1;
    }
    return 0;
}

} // attributes

unittest
{
    double x, y;
    twoSum(1.0, 0x1p-60, x, y);
    assert(x == 1.0 && y == 0x1p-60);
    twoDiff(1.0, 0x1p-60, x, y);
    assert(x == 1.0 && y == -0x1p-60);
    twoProduct(1.0 + 0x1p-30, 1.0 + 0x1p-30, x, y);
    assert(x == 1.0 + 0x1p-29 && y == 0x1p-60);
}

unittest
{
    // (1 + 2^-60) - 1 is not representable as a single double sum but the
    // expansion keeps it exactly.
    auto e = sum([0x1p-60, 1.0], [-1.0]);
    assert(e == [0x1p-60]);
    assert(signOf(difference([1.0], [1.0])) == 0);
    assert(signOf(difference([0x1p-80, 1.0], [1.0])) == 1);
    auto p = product(diffExpansion(1.0, 0x1p-60), diffExpansion(1.0, 0x1p-60));
    // (1 - 2^-60)^2 = 1 - 2^-59 + 2^-120
    assert(p[0] == 0x1p-120);
    assert(signOf(difference(p, sum(sum([0x1p-120], [-0x1p-59]), [1.0]))) == 0);
}

unittest
{
    import std.bigint : BigInt;
    import std.random : Mt19937, uniform;

    // Compare expansion products/sums of integer-valued doubles with BigInt.
    // Integer-valued double -> BigInt (handles magnitudes beyond long).
    static BigInt bigOf(double v)
    {
        int shift = 0;
        while (v >= 0x1p62 || v <= -0x1p62) { v /= 2; ++shift; }
        BigInt r = BigInt(cast(long) v);
        foreach (_; 0 .. shift) r *= 2;
        return r;
    }

    static BigInt toBig(const(double)[] e)
    {
        BigInt r = 0;
        foreach (c; e)
            r += bigOf(c);
        return r;
    }

    auto rng = Mt19937(12345);
    foreach (_; 0 .. 2000)
    {
        const a = cast(double) uniform(-(1L << 40), 1L << 40, rng);
        const b = cast(double) uniform(-(1L << 40), 1L << 40, rng);
        const c = cast(double) uniform(-(1L << 40), 1L << 40, rng);
        auto ab = productExpansion(a, b);
        auto expected = BigInt(cast(long) a) * BigInt(cast(long) b);
        assert(toBig(ab) == expected);
        auto abc = product(ab, [c]);
        // components may exceed long range; check sign and value through BigInt sums of halves
        BigInt got = toBig(abc);
        assert(got == expected * BigInt(cast(long) c));
        const s = sum(ab, productExpansion(c, a));
        assert(signOf(s) == ((expected + BigInt(cast(long) c) * BigInt(cast(long) a)) > 0 ? 1
                : (expected + BigInt(cast(long) c) * BigInt(cast(long) a)) < 0 ? -1 : 0));
    }
}
