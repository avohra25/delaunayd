# delaunayd

Robust 2D computational geometry in D, built on exact geometric predicates.

Work in progress. Currently included:

- **Exact expansion arithmetic** (`delaunayd.expansion`): Shewchuk's
  error-free transformations (`twoSum`, `twoProduct`, ...) and
  zero-eliminating expansion sum/scale/product.
- **Robust predicates** (`delaunayd.predicates`): adaptive `orient2d`
  (floating-point filter, then stages B/C, then exact) and `incircle`
  (filter, then exact expansion evaluation). Signs are always exact.
- **Convex hull** (`delaunayd.hull`): Andrew's monotone chain on top of the
  exact predicates, handling duplicates and collinear points; cross-checked
  against an O(n^3) brute force on random and degenerate inputs.
- **Delaunay triangulation** (`delaunayd.triangulation`): lexicographic
  sweep insertion with Lawson flips over a flat half-edge structure
  (`triangles`/`halfedges` arrays), duplicate detection, convex hull with
  boundary vertices, vertex stars, and stochastic-walk point location.
  Verified with an empty-circumcircle brute force on thousands of seeded
  random, grid (cocircular), circular, collinear and duplicate-heavy inputs.
- **Constrained edges** (`delaunayd.constrained`): segment tracing,
  Sloan's flip-based crossing removal, splitting at collinear vertices,
  atomic rejection of crossing constraints, and constrained Lawson flips;
  checked with a brute-force visibility/empty-circle test.
- **Voronoi diagram** (`delaunayd.voronoi`): circumcenters, edges/rays per
  Delaunay edge, CCW cell vertex lists, and box-clipped cells by bisector
  half-plane intersection (also for collinear input); checked against
  brute-force nearest-site queries and total-area conservation.

```sh
dub test
```

Phobos only, no dependencies. MIT licensed.
