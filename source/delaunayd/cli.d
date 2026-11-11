/**
 * The `delaunayd` command-line tool.
 *
 * ---
 * delaunayd triangulate [--svg] [-o FILE] [FILE]
 * delaunayd hull        [--svg] [-o FILE] [FILE]
 * delaunayd voronoi     [--svg] [--margin M] [--delaunay] [-o FILE] [FILE]
 * delaunayd locate X Y  [FILE]
 * delaunayd random N    [--seed S] [--dist uniform|grid|circle|normal] [-o FILE]
 * ---
 */
module delaunayd.cli;

import std.array : appender, Appender;
import std.conv : to, ConvException;
import std.format : format, formattedWrite;
import std.stdio : File, stdin, stdout, stderr;

import delaunayd.point;
import delaunayd.predicates;
import delaunayd.hull;
import delaunayd.triangulation;
import delaunayd.constrained;
import delaunayd.voronoi;
import delaunayd.io;
import delaunayd.svg;

enum string versionString = "delaunayd 0.1.0";

enum string usage = `usage: delaunayd <command> [options] [FILE]

Robust 2D Delaunay triangulation, convex hulls and Voronoi diagrams.

commands:
  triangulate      Delaunay triangulation (constrained if the input has
                   'seg i j' lines)
  hull             convex hull (strict corners, counterclockwise)
  voronoi          Voronoi diagram with cells clipped to a box
  locate X Y       locate a point in the triangulation of the input
  random N         print N random points
  help, version

options:
  -o, --output FILE   write to FILE instead of standard output
  --svg               write SVG instead of text
  --margin M          voronoi: clip box margin as a fraction of the input
                      extent (default 0.1)
  --delaunay          voronoi --svg: overlay the Delaunay edges
  --seed S            random: RNG seed (default 1)
  --dist D            random: uniform | grid | circle | normal

input: one point per line ("x y" or "x,y"), 'seg i j' for a constraint
between points i and j (0-based), '#' comments. FILE '-' or none = stdin.
`;

private class UsageError : Exception
{
    this(string msg, string file = __FILE__, size_t line = __LINE__) @safe pure nothrow
    {
        super(msg, file, line);
    }
}

private struct Options
{
    string command;
    string[] positional;
    string output;
    bool svg;
    bool overlay;
    double margin = 0.1;
    uint seed = 1;
    string dist = "uniform";
}

private Options parseArgs(string[] args)
{
    Options o;
    if (args.length < 2)
        throw new UsageError("missing command");
    o.command = args[1];
    for (size_t i = 2; i < args.length; ++i)
    {
        const a = args[i];
        string value()
        {
            if (i + 1 >= args.length)
                throw new UsageError(a ~ " needs a value");
            return args[++i];
        }

        try
        {
            switch (a)
            {
            case "-o", "--output": o.output = value(); break;
            case "--svg": o.svg = true; break;
            case "--delaunay": o.overlay = true; break;
            case "--margin": o.margin = value().to!double; break;
            case "--seed": o.seed = value().to!uint; break;
            case "--dist": o.dist = value(); break;
            default:
                if (a.length > 1 && a[0] == '-' && !isNumber(a))
                    throw new UsageError("unknown option " ~ a);
                o.positional ~= a;
            }
        }
        catch (ConvException)
            throw new UsageError("invalid value for " ~ a);
    }
    return o;
}

private bool isNumber(string s)
{
    try
    {
        s.to!double;
        return true;
    }
    catch (ConvException)
        return false;
}

/**
 * Runs the tool. `readFile` loads an input path ("-" = stdin); output goes
 * to `outp` (or the -o file), diagnostics to `errp`. Returns the exit code.
 */
int run(string[] args, ref Appender!string outp, ref Appender!string errp,
        string delegate(string path) readFile)
{
    Options o;
    try
        o = parseArgs(args);
    catch (UsageError e)
    {
        errp.put("delaunayd: " ~ e.msg ~ "\n\n" ~ usage);
        return 2;
    }
    try
    {
        string result;
        Input loadInput(size_t positionalIndex)
        {
            const path = o.positional.length > positionalIndex ? o.positional[positionalIndex] : "-";
            return parseInput(readFile(path));
        }

        switch (o.command)
        {
        case "help", "--help", "-h":
            outp.put(usage);
            return 0;
        case "version", "--version":
            outp.put(versionString ~ "\n");
            return 0;
        case "triangulate":
            result = cmdTriangulate(loadInput(0), o);
            break;
        case "hull":
            result = cmdHull(loadInput(0), o);
            break;
        case "voronoi":
            result = cmdVoronoi(loadInput(0), o);
            break;
        case "locate":
            if (o.positional.length < 2)
                throw new UsageError("locate needs X and Y");
            const q = Point(o.positional[0].to!double, o.positional[1].to!double);
            result = cmdLocate(loadInput(2), q);
            break;
        case "random":
            if (o.positional.length < 1)
                throw new UsageError("random needs a point count");
            result = cmdRandom(o.positional[0].to!size_t, o);
            break;
        default:
            throw new UsageError("unknown command " ~ o.command);
        }
        if (o.output.length)
            File(o.output, "w").write(result);
        else
            outp.put(result);
        return 0;
    }
    catch (UsageError e)
    {
        errp.put("delaunayd: " ~ e.msg ~ "\n\n" ~ usage);
        return 2;
    }
    catch (ConvException e)
    {
        errp.put("delaunayd: invalid number: " ~ e.msg ~ "\n");
        return 2;
    }
    catch (Exception e)
    {
        errp.put("delaunayd: " ~ e.msg ~ "\n");
        return 1;
    }
}

/// Entry point used by `main`.
int runCli(string[] args)
{
    auto outp = appender!string(), errp = appender!string();
    const code = run(args, outp, errp, (string path) {
        import std.file : readText;

        if (path == "-")
        {
            auto buf = appender!string();
            foreach (chunk; stdin.byChunk(1 << 16))
                buf.put(cast(const(char)[]) chunk);
            return buf[];
        }
        return readText(path);
    });
    stdout.write(outp[]);
    stderr.write(errp[]);
    return code;
}

// ------------------------------------------------------------- commands

private Triangulation buildTriangulation(const Input inp)
{
    auto tr = new Triangulation(inp.points);
    insertConstraints(tr, inp.segments);
    return tr;
}

private string cmdTriangulate(const Input inp, const Options o)
{
    auto tr = buildTriangulation(inp);
    const P = tr.points;
    if (o.svg)
    {
        auto c = SvgCanvas(Box.around(P, 0.05));
        foreach (t; 0 .. tr.triangleCount)
        {
            const v = tr.triangle(t);
            c.polygon([P[v[0]], P[v[1]], P[v[2]]], "#eef4ff", "#5577aa", 0.8);
        }
        if (tr.isDegenerate)
            foreach (e; tr.edges())
                c.line(P[e[0]], P[e[1]], "#5577aa", 0.8);
        foreach (e; tr.constrainedEdges())
            c.line(P[e[0]], P[e[1]], "#d62728", 2.5);
        foreach (v; tr.vertices)
            c.dot(P[v]);
        return c.finish();
    }
    auto w = appender!string();
    w.formattedWrite("# delaunayd triangulate: %s points (%s distinct), %s triangles, %s edges, %s constrained\n",
            P.length, tr.vertices.length, tr.triangleCount, tr.edges().length, tr.constraintCount);
    w.formattedWrite("points %s\n", P.length);
    foreach (p; P)
        w.formattedWrite("%s %s\n", fmt(p.x), fmt(p.y));
    w.formattedWrite("triangles %s\n", tr.triangleCount);
    foreach (t; 0 .. tr.triangleCount)
    {
        const v = tr.triangle(t);
        w.formattedWrite("%s %s %s\n", v[0], v[1], v[2]);
    }
    const ce = tr.constrainedEdges();
    w.formattedWrite("constraints %s\n", ce.length);
    foreach (e; ce)
        w.formattedWrite("%s %s\n", e[0], e[1]);
    return w[];
}

private string cmdHull(const Input inp, const Options o)
{
    const P = inp.points;
    auto h = convexHull(P);
    Point[] poly;
    foreach (i; h)
        poly ~= P[i];
    if (o.svg)
    {
        auto c = SvgCanvas(Box.around(P, 0.05));
        c.polygon(poly, "#e6f5e6", "#2a7a2a", 2);
        foreach (p; P)
            c.dot(p, 2, "#555");
        foreach (p; poly)
            c.dot(p, 3.5, "#2a7a2a");
        return c.finish();
    }
    auto w = appender!string();
    w.formattedWrite("# delaunayd hull: %s vertices of %s points, area %s\n", h.length, P.length,
            fmt(polygonArea(poly)));
    foreach (i; h)
        w.formattedWrite("%s %s %s\n", i, fmt(P[i].x), fmt(P[i].y));
    return w[];
}

private string cmdVoronoi(const Input inp, const Options o)
{
    if (inp.segments.length)
        throw new Exception("voronoi does not support constraint segments");
    auto tr = new Triangulation(inp.points);
    auto vd = new Voronoi(tr);
    const P = tr.points;
    const box = Box.around(P, o.margin);
    auto cells = vd.clippedCells(box);
    if (o.svg)
    {
        auto c = SvgCanvas(box);
        foreach (k, v; tr.vertices)
            c.polygon(cells[v], pastel(k), "#333", 1);
        if (o.overlay)
            foreach (e; tr.edges())
                c.line(P[e[0]], P[e[1]], "#888", 0.6, "3 2");
        foreach (v; tr.vertices)
            c.dot(P[v]);
        return c.finish();
    }
    auto w = appender!string();
    auto edges = vd.edges();
    w.formattedWrite("# delaunayd voronoi: %s sites, %s vertices, %s edges; box %s %s %s %s\n",
            tr.vertices.length, vd.centers.length, edges.length, fmt(box.minX), fmt(box.minY),
            fmt(box.maxX), fmt(box.maxY));
    w.formattedWrite("vertices %s\n", vd.centers.length);
    foreach (p; vd.centers)
        w.formattedWrite("%s %s\n", fmt(p.x), fmt(p.y));
    w.formattedWrite("edges %s\n", edges.length);
    foreach (e; edges)
    {
        final switch (e.kind)
        {
        case VoronoiEdge.Kind.segment:
            w.formattedWrite("segment %s %s %s %s %s %s\n", e.siteA, e.siteB, fmt(e.from.x),
                    fmt(e.from.y), fmt(e.to.x), fmt(e.to.y));
            break;
        case VoronoiEdge.Kind.ray:
            w.formattedWrite("ray %s %s %s %s %s %s\n", e.siteA, e.siteB, fmt(e.from.x),
                    fmt(e.from.y), fmt(e.dir.x), fmt(e.dir.y));
            break;
        case VoronoiEdge.Kind.line:
            w.formattedWrite("line %s %s %s %s %s %s\n", e.siteA, e.siteB, fmt(e.from.x),
                    fmt(e.from.y), fmt(e.dir.x), fmt(e.dir.y));
            break;
        }
    }
    w.formattedWrite("cells %s\n", tr.vertices.length);
    foreach (v; tr.vertices)
    {
        w.formattedWrite("site %s %s", v, cells[v].length);
        foreach (p; cells[v])
            w.formattedWrite(" %s %s", fmt(p.x), fmt(p.y));
        w.put("\n");
    }
    return w[];
}

private string cmdLocate(const Input inp, Point q)
{
    auto tr = buildTriangulation(inp);
    const loc = tr.locate(q);
    final switch (loc.kind)
    {
    case Location.Kind.outside:
        return "outside\n";
    case Location.Kind.onVertex:
        return format("vertex %s\n", loc.vertex);
    case Location.Kind.onEdge:
        return format("edge %s %s\n", tr.triangles[loc.edge], tr.triangles[nextHalfedge(loc.edge)]);
    case Location.Kind.inTriangle:
        const v = tr.triangle(loc.triangle);
        return format("triangle %s %s %s %s\n", loc.triangle, v[0], v[1], v[2]);
    }
}

private string cmdRandom(size_t n, const Options o)
{
    import std.random : Mt19937, uniform, uniform01;
    import std.math : cos, sin, sqrt, log, PI;

    auto rng = Mt19937(o.seed);
    auto w = appender!string();
    w.formattedWrite("# delaunayd random %s --dist %s --seed %s\n", n, o.dist, o.seed);
    foreach (_; 0 .. n)
    {
        Point p;
        switch (o.dist)
        {
        case "uniform":
            p = Point(uniform01(rng), uniform01(rng));
            break;
        case "grid":
            p = Point(uniform(0, 20, rng), uniform(0, 20, rng));
            break;
        case "circle":
            const a = uniform01(rng) * 2 * PI;
            p = Point(cos(a), sin(a));
            break;
        case "normal":
            // Box-Muller
            const u1 = uniform!"(]"(0.0, 1.0, rng), u2 = uniform01(rng);
            const r = sqrt(-2 * log(u1));
            p = Point(r * cos(2 * PI * u2), r * sin(2 * PI * u2));
            break;
        default:
            throw new UsageError("unknown distribution " ~ o.dist);
        }
        w.formattedWrite("%s %s\n", fmt(p.x), fmt(p.y));
    }
    return w[];
}

// ------------------------------------------------------------------ tests

version (unittest)
{
    import std.algorithm : canFind, startsWith;
    import std.string : splitLines;

    struct Result
    {
        int code;
        string output, error;
    }

    Result invoke(string[] args, string input = "")
    {
        auto outp = appender!string(), errp = appender!string();
        const code = run(["delaunayd"] ~ args, outp, errp, (string path) {
            if (path == "-")
                return input;
            import std.file : readText;

            return readText(path);
        });
        return Result(code, outp[], errp[]);
    }
}

unittest
{
    const square = "0 0\n1 0\n1 1\n0 1\n0.5 0.5\n";
    auto r = invoke(["triangulate"], square);
    assert(r.code == 0);
    assert(r.output.startsWith("# delaunayd triangulate: 5 points (5 distinct), 4 triangles, 8 edges, 0 constrained"));
    auto lines = r.output.splitLines;
    assert(lines[1] == "points 5" && lines[7] == "triangles 4" && lines[$ - 1] == "constraints 0");

    r = invoke(["hull"], square);
    assert(r.code == 0 && r.output.canFind("4 vertices of 5 points, area 1"));
    assert(r.output.splitLines[1] == "0 0 0");

    r = invoke(["locate", "0.25", "0.5"], square);
    assert(r.output.startsWith("triangle "));
    assert(invoke(["locate", "0.5", "0.5"], square).output == "vertex 4\n");
    assert(invoke(["locate", "5", "5"], square).output == "outside\n");
    assert(invoke(["locate", "0.5", "0"], square).output == "edge 0 1\n"
            || invoke(["locate", "0.5", "0"], square).output == "edge 1 0\n");
}

unittest
{
    // constrained input and SVG output
    const pts = "0 0\n4 0\n4 1\n0 1\n2 0.5\n1 0.9\n3 0.1\nseg 0 2\n";
    auto r = invoke(["triangulate"], pts);
    assert(r.code == 0, r.error);
    assert(r.output.canFind("2 constrained")); // split at the collinear vertex (2, 0.5)
    r = invoke(["triangulate", "--svg"], pts);
    assert(r.output.startsWith("<svg") && r.output.canFind(`stroke="#d62728"`));
    r = invoke(["voronoi"], "0 0\n2 0\n0 2\n2 2\n1 1\n");
    assert(r.code == 0);
    assert(r.output.canFind("5 sites, 4 vertices, 8 edges"));
    assert(r.output.canFind("\ncells 5\n"));
    r = invoke(["voronoi", "--svg", "--delaunay"], "0 0\n2 0\n0 2\n2 2\n1 1\n");
    assert(r.output.canFind("<polygon") && r.output.canFind("stroke-dasharray"));
    r = invoke(["hull", "--svg"], "0 0\n1 0\n0 1\n");
    assert(r.output.startsWith("<svg"));
}

unittest
{
    // random generation is deterministic and parseable
    auto a = invoke(["random", "50", "--seed", "7", "--dist", "normal"]);
    auto b = invoke(["random", "50", "--seed", "7", "--dist", "normal"]);
    assert(a.code == 0 && a.output == b.output);
    assert(parseInput(a.output).points.length == 50);
    foreach (d; ["uniform", "grid", "circle"])
        assert(parseInput(invoke(["random", "10", "--dist", d]).output).points.length == 10);
}

unittest
{
    // errors
    assert(invoke([]).code == 2);
    assert(invoke(["frobnicate"]).code == 2);
    assert(invoke(["triangulate", "--bogus"]).code == 2);
    assert(invoke(["random"]).code == 2);
    assert(invoke(["random", "5", "--dist", "zipf"]).code == 2);
    auto r = invoke(["triangulate"], "0 0\n1 2 3\n");
    assert(r.code == 1 && r.error.canFind("line 2"));
    // crossing constraints
    r = invoke(["triangulate"], "0 0\n1 1\n1 0\n0 1\nseg 0 1\nseg 2 3\n");
    assert(r.code == 1 && r.error.canFind("crosses"));
    assert(invoke(["voronoi"], "0 0\n1 0\n0 1\nseg 0 1\n").code == 1);
    assert(invoke(["help"]).output.canFind("usage:"));
    assert(invoke(["version"]).output == versionString ~ "\n");
    // negative numbers are positional values, not options
    assert(invoke(["locate", "-1", "-1"], "0 0\n1 0\n0 1\n").output == "outside\n");
}
