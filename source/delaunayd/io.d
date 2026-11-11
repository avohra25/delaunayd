/**
 * Plain-text point/segment input.
 *
 * One point per line as `x y` (whitespace and/or comma separated);
 * `seg i j` (or `s i j`) lines add a constraint between points i and j
 * (0-based, in file order). `#` starts a comment; blank lines are ignored.
 */
module delaunayd.io;

import std.algorithm : splitter, filter;
import std.array : array;
import std.ascii : isWhite;
import std.conv : to, ConvException;
import std.format : format;
import std.math : isFinite;
import std.string : indexOf, strip;

import delaunayd.point;

/// Parsed input: points and optional constraint segments.
struct Input
{
    Point[] points;
    int[2][] segments;
}

/// Parse error with a 1-based line number.
class ParseError : Exception
{
    size_t line;

    this(size_t line, string msg, string file = __FILE__, size_t l = __LINE__) @safe pure
    {
        this.line = line;
        super(format("line %s: %s", line, msg), file, l);
    }
}

/// Parses the text format described in the module documentation.
Input parseInput(string text)
{
    Input r;
    size_t lineNo = 0;
    foreach (raw; text.splitter('\n'))
    {
        ++lineNo;
        auto line = raw;
        const hash = line.indexOf('#');
        if (hash >= 0)
            line = line[0 .. hash];
        auto fields = line.splitter!(c => isWhite(c) || c == ',')
            .filter!(f => f.length > 0).array;
        if (fields.length == 0)
            continue;
        if (fields[0] == "seg" || fields[0] == "s")
        {
            if (fields.length != 3)
                throw new ParseError(lineNo, "expected `seg i j`");
            int i, j;
            try
            {
                i = fields[1].to!int;
                j = fields[2].to!int;
            }
            catch (ConvException)
                throw new ParseError(lineNo, "segment endpoints must be integers");
            r.segments ~= [i, j];
            continue;
        }
        if (fields.length != 2)
            throw new ParseError(lineNo, format("expected 2 coordinates, got %s field(s)", fields.length));
        double x, y;
        try
        {
            x = fields[0].strip.to!double;
            y = fields[1].strip.to!double;
        }
        catch (ConvException)
            throw new ParseError(lineNo, "invalid number");
        if (!isFinite(x) || !isFinite(y))
            throw new ParseError(lineNo, "coordinates must be finite");
        r.points ~= Point(x, y);
    }
    foreach (k, s; r.segments)
        foreach (v; s)
            if (v < 0 || v >= r.points.length)
                throw new ParseError(0, format("segment %s refers to point %s, but there are %s points",
                        k, v, r.points.length));
    return r;
}

/// Formats a double so that it parses back to the same value.
string fmt(double v)
{
    return format("%.17g", v);
}

unittest
{
    auto inp = parseInput("# demo\n0 0\n1,0\n  0.5   2.5e-1 # inline\n\nseg 0 1\ns 1 2\n");
    assert(inp.points == [Point(0, 0), Point(1, 0), Point(0.5, 0.25)]);
    assert(inp.segments == [[0, 1], [1, 2]]);
    assert(parseInput("").points.length == 0);
    assert(parseInput("1 2\r\n3 4\r\n").points == [Point(1, 2), Point(3, 4)]);
}

unittest
{
    import std.exception : assertThrown, collectException;

    assertThrown!ParseError(parseInput("1 2 3\n"));
    assertThrown!ParseError(parseInput("1 x\n"));
    assertThrown!ParseError(parseInput("nan 1\n"));
    assertThrown!ParseError(parseInput("0 0\nseg 0\n"));
    assertThrown!ParseError(parseInput("0 0\nseg 0 5\n"));
    auto e = collectException!ParseError(parseInput("0 0\n\n1 2 3\n"));
    assert(e !is null && e.line == 3);
    // round-trip formatting
    foreach (v; [0.1, 1.0 / 3, 1e-300, -2.5e17])
        assert(fmt(v).to!double == v);
}
