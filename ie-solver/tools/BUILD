package(default_visibility = ["//visibility:public"])

cc_library(
	name = "ie_solver_tools",
	hdrs = ["ie_solver_tools.h"],
	srcs = ["ie_solver_tools.cpp",
			"schur_update.cpp",
			"compression.cpp",
			"boundary_to_domain.cpp",
			"proxy.cpp"],
	deps = ["//ie-solver:ie_mat",
			"//ie-solver:quadtree",
			"//ie-solver:kernel"]
)