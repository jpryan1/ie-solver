# ie-solver
## Dependencies:

* bazel, (industrial open source fast compiler and linker)

* blas, lapack, lapacke, (lin alg libs)

* boost (getting integration weights for parametric boundaries of certain types requires elliptic integral evaluation, but maybe you can avoid the need for a heavy lib by just commenting out that code in the BUILD files and in the includes)

