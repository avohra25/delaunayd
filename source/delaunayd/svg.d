/// Minimal SVG rendering of points, triangulations, hulls and Voronoi cells.
module delaunayd.svg;

import std.array : appender, Appender;
import std.format : formattedWrite;

import delaunayd.point;
import delaunayd.voronoi : Box;

/// Accumulates SVG elements in world coordinates, mapped into a fixed-width
/// canvas with y pointing up.
struct SvgCanvas
{
    private Box box;
    private double scale;
    private double width, height;
    private Appender!string body_;

    this(Box world, double widthPx = 800)
    {
        box = world;
        const w = world.maxX - world.minX, h = world.maxY - world.minY;
        width = widthPx;
        scale = w > 0 ? widthPx / w : 1;
        height = h * scale;
        if (height < 1)
            height = 1;
        body_ = appender!string();
    }

    private double sx(double x) const { return (x - box.minX) * scale; }
    private double sy(double y) const { return (box.maxY - y) * scale; }

    void polygon(const(Point)[] pts, string fill, string stroke, double strokeWidth = 1)
    {
        if (pts.length == 0)
            return;
        body_.put(`<polygon points="`);
        foreach (i, p; pts)
            body_.formattedWrite("%s%.3f,%.3f", i ? " " : "", sx(p.x), sy(p.y));
        body_.formattedWrite(`" fill="%s" stroke="%s" stroke-width="%s"/>`, fill, stroke, strokeWidth);
        body_.put("\n");
    }

    void line(Point a, Point b, string stroke, double strokeWidth = 1, string dash = null)
    {
        body_.formattedWrite(`<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="%s" stroke-width="%s"`,
                sx(a.x), sy(a.y), sx(b.x), sy(b.y), stroke, strokeWidth);
        if (dash.length)
            body_.formattedWrite(` stroke-dasharray="%s"`, dash);
        body_.put("/>\n");
    }

    void dot(Point p, double r = 2.5, string fill = "#111")
    {
        body_.formattedWrite(`<circle cx="%.3f" cy="%.3f" r="%s" fill="%s"/>`, sx(p.x), sy(p.y), r, fill);
        body_.put("\n");
    }

    /// The complete SVG document.
    string finish() const
    {
        auto o = appender!string();
        o.formattedWrite(`<svg xmlns="http://www.w3.org/2000/svg" width="%.0f" height="%.0f" viewBox="0 0 %.3f %.3f">`,
                width, height, width, height);
        o.put("\n");
        o.formattedWrite(`<rect width="100%%" height="100%%" fill="white"/>`);
        o.put("\n");
        o.put(body_[]);
        o.put("</svg>\n");
        return o[];
    }
}

/// A pleasant pastel fill for index i (deterministic).
string pastel(size_t i)
{
    import std.format : format;

    // golden-angle hue walk
    const hue = (i * 137.508) % 360.0;
    return format("hsl(%.0f,70%%,85%%)", hue);
}

unittest
{
    import std.algorithm : canFind;

    auto c = SvgCanvas(Box(0, 0, 10, 5), 200);
    c.polygon([Point(0, 0), Point(10, 0), Point(10, 5)], "none", "black");
    c.line(Point(0, 0), Point(10, 5), "red", 2, "4 2");
    c.dot(Point(5, 2.5));
    const s = c.finish();
    assert(s.canFind(`width="200" height="100"`));
    assert(s.canFind(`<polygon points="0.000,100.000 200.000,100.000 200.000,0.000"`));
    assert(s.canFind(`stroke-dasharray="4 2"`));
    assert(s.canFind(`<circle cx="100.000" cy="50.000"`));
    assert(pastel(0) != pastel(1));
}
