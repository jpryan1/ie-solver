package(default_visibility = ["//visibility:public"])

cc_library(
	name = "skel_factorization",
	hdrs = ["skel_factorization.h"],
	srcs = ["skel_factorization.cpp"],
  copts = ["-fopenmp"],
	deps = ["//ie_solver:ie_mat",
			"//ie_solver/quadtree",
			"//ie_solver/kernel"]
)