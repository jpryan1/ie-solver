package(default_visibility = ["//visibility:public"])

cc_binary(
  name = "experiment_one",
  srcs = ["experiment_one.cpp"],
  copts= ["-fopenmp"],
  deps = [
      "//ie_solver:ie_mat",
      "//ie_solver:initialization",
      "//ie_solver/skel_factorization",
      "//ie_solver/quadtree",
      "//ie_solver/boundaries",
      "//ie_solver:linear_solve",
      ],
  linkopts = ["-llapacke",
        "-llapack",
        "-lblas",
        "-lgomp"]

)


cc_binary(
  name = "experiment_two",
  srcs = ["experiment_two.cpp"],
  copts= ["-fopenmp"],
  deps = [
      "//ie_solver:ie_mat",
      "//ie_solver:initialization",
      "//ie_solver/skel_factorization",
      "//ie_solver/quadtree",
      "//ie_solver/boundaries",
      "//ie_solver:linear_solve",
      ],
  linkopts = ["-llapacke",
        "-llapack",
        "-lblas",
        "-lgomp"]

)


cc_binary(
  name = "experiment_three",
  srcs = ["experiment_three.cpp"],
  copts= ["-fopenmp"],
  deps = [
      "//ie_solver:ie_mat",
      "//ie_solver:initialization",
      "//ie_solver/skel_factorization",
      "//ie_solver/quadtree",
      "//ie_solver/boundaries",
      "//ie_solver:linear_solve",
      ],
  linkopts = ["-llapacke",
        "-llapack",
        "-lblas",
        "-lgomp"]

)

cc_binary(
  name = "landscapes",
  srcs = ["landscapes.cpp"],
  copts= ["-fopenmp"],
  deps = [
      "//ie_solver:ie_mat",
      "//ie_solver:initialization",
      "//ie_solver/skel_factorization",
      "//ie_solver/quadtree",
      "//ie_solver/boundaries",
      "//ie_solver:linear_solve",
      ],
  linkopts = ["-llapacke",
        "-llapack",
        "-lblas",
        "-lgomp"]

)
