cc_binary(
	name = "ie-solver",
	srcs = ["ie-solver.cpp"],
	deps = [":common",
			":initialization",
			":ie_solver_tools",
			":quadtree",
			":circle",
			":rounded_square",
			],
	linkopts = ["-llapacke",
				"-llapack",
				"-lblas",
				"-lgomp"]

)

cc_library(
	name = "common",
	hdrs = ["common.h"],
	srcs = ["vec2.h", 
			"log.h", 
			"interaction_lists.h", 
			"clock.h",
			"matmul.cpp",
			"ie_mat.cpp", 
		],
)


cc_library(
	name = "initialization",
	hdrs = ["initialization.h"],
	srcs = ["initialization.cpp",
			"stokes_init.cpp",
			],
	deps = [":common",
			":quadtree"]
)


cc_library(
	name = "ie_solver_tools",
	hdrs = ["ie_solver_tools.h"],
	srcs = ["ie_solver_tools.cpp"],
	deps = [":common",
			":quadtree"]
)


cc_library(
	name = "quadtree",
	hdrs = ["quadtree.h"],
	srcs = ["quadtree.cpp"],
	deps = [":common"]
)

cc_library(
	name = "circle",
	hdrs = ["circle.h"],
	srcs = ["circle.cpp"],
	deps = [":common"]
)

cc_library(
	name = "rounded_square",
	hdrs = ["rounded_square.h"],
	srcs = ["rounded_square.cpp"],
	deps = [":common"]
)
