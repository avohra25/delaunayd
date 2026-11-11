/// Basic 2D point type.
module delaunayd.point;

import std.math : isFinite;

/// A point (or vector) in the plane.
struct Point
{
    double x = 0.0;
    double y = 0.0;

    Point opBinary(string op)(Point o) const @safe pure nothrow @nogc
            if (op == "+" || op == "-")
    {
        return Point(mixin("x " ~ op ~ " o.x"), mixin("y " ~ op ~ " o.y"));
    }

    Point opBinary(string op : "*")(double s) const @safe pure nothrow @nogc
    {
        return Point(x * s, y * s);
    }

    /// Lexicographic order: by x, then by y.
    int opCmp(Point o) const @safe pure nothrow @nogc
    {
        if (x < o.x) return -1;
        if (x > o.x) return 1;
        if (y < o.y) return -1;
        if (y > o.y) return 1;
        return 0;
    }

    bool isValid() const @safe pure nothrow @nogc
    {
        return isFinite(x) && isFinite(y);
    }
}

double dot(Point a, Point b) @safe pure nothrow @nogc
{
    return a.x * b.x + a.y * b.y;
}

double cross(Point a, Point b) @safe pure nothrow @nogc
{
    return a.x * b.y - a.y * b.x;
}

double dist2(Point a, Point b) @safe pure nothrow @nogc
{
    const dx = a.x - b.x, dy = a.y - b.y;
    return dx * dx + dy * dy;
}

unittest
{
    assert(Point(1, 2) < Point(1, 3));
    assert(Point(0, 9) < Point(1, -9));
    assert(Point(1, 2) + Point(3, 4) == Point(4, 6));
    assert(Point(1, 2) * 2 == Point(2, 4));
    assert(cross(Point(1, 0), Point(0, 1)) == 1);
    assert(!Point(double.nan, 0).isValid);
}
