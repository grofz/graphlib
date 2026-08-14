### Graph library
A straightforward graph library in Fortran

### Features

- supports directed or undirected graphs
- supports weighted or unweighted graphs
- betweenness centrality index calculation
- network conductance calculation
- network maximum flow, minimum ST-cut 
- find connected components
- find strongly connected components, topological ordering
- following algorithms are used:
  - aaa
  - bbb

### Quick start

1. Edit constants in "src/graph_user.f90" to set the size of storage space for
  edges and vertices. The storage space must be known at compile time.
  This design choice is to avoid parametrized defined types and pointers.
  Also, fixed size storage may be faster.

2. Compile library
```console
$ make
```

3. Compile and run tests and an example (optional)
```console
$ make test
$ bin/test_betweenness
$ bin/test_concom
$ bin/test_consts
$ bin/test_maxflow
$ bin/test_scc
$ bin/example
```

4. Use in your code
```console
$ gfortran mycode.f90 -I./include -L. -lgraph -o mycode
```
Link with "./libgraph.a" and use module files in "./include" directory.

### Example & Documentation

The quick use example is in "testsrc/example.f90".

All graph procedures are exported as type bound procedures of "graph_t" type.
Please consult "src/graph.f90" for documentation of the procedures.

The storage size can be set by redefining constants in "src/graph_user.f90" and
recompiling the library.

A simple read/write procedure (VTU format) is in "src/vtuio.f90"

The remaining files in "src/" are intended for internal use.


### License

GNU license.
