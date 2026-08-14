### Graph library

A straightforward graph library written in modern Fortran.

![The flow and potential visualization in a simple network](http://github.com/grofz/ds/assets/a.png)
![The flow in a simple network](http://github.com/grofz/ds/assets/b.png)
![Betweenness centrality of a simple network](http://github.com/grofz/ds/assets/c.png)

### Features

- directed and undirected graphs
- weighted and unweighted graphs
- connected and strongly connected components
- topological ordering and levels
- shortest-path search
- betweenness centrality
- network conductance
- maximum flow and minimum s-t cut
- following methods are used:
  - Brandes algorithm (betweenness)
  - Tarjan and Kahn algorithms
  - Stoer-Wagner, Edmond-Karp and Dinics algorithms (maximum flow)
  - Conjugate-Gradient method

### Quick start

The size of arrays storing vertice and edge data are defined at compile time. Edit the
constants in `src/graph_user.f90` if necessary.

Build the library with:

```console
$ make
```

Optionally build and run the tests and example:

```console
$ make test
$ bin/example
```

To use the library in your own program:

```console
$ gfortran mycode.f90 -I./include -L. -lgraph -o mycode
```

The library is provided as `libgraph.a`; the required module files are in
`include/`.

### Example & Documentation

A short example demonstrating the basic API is provided in
`testsrc/example.f90`.

Graph algorithms and other public procedures are available as type-bound
procedures of `graph_t`. Their interfaces and documentation are provided in
`src/graph.f90`.

The storage sizes can be changed in `src/graph_user.f90` and the library
recompiled.

A simple VTU-format read/write utility is provided in `src/vtuio.f90`.

Other files in `src/` are intended for internal use.

### License

Graph is free software distributed under the GNU Lesser General Public
License, version 3 or later (LGPL-3.0-or-later).
See the LICENSE file for the full licence text.