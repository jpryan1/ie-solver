package(default_visibility = ["//visibility:public"])

cc_binary(
	name = "ie-solver",
	srcs = ["ie-solver.cpp"],
	deps = [":ie_mat",
			":initialization",
			"//ie-solver/tools:ie_solver_tools",
			":quadtree",
			"//ie-solver/boundaries",
			":helpers",
			],
	linkopts = ["-llapacke",
				"-llapack",
				"-lblas",
				"-lgomp"]

)


cc_binary(
	name = "updating-ie-solver",
	srcs = ["updating-ie-solver.cpp"],
	deps = [":ie_mat",
			":initialization",
			"//ie-solver/tools:ie_solver_tools",
			":quadtree",
			"//ie-solver/boundaries",
			":helpers",
			],
	linkopts = ["-llapacke",
				"-llapack",
				"-lblas",
				"-lgomp"]

)


cc_library(
	name = "helpers",
	hdrs = ["helpers.h"],
	srcs = ["helpers.cpp"],
	deps = [":ie_mat",
			":log",
			"//ie-solver/boundaries"],
)


cc_library(
	name = "ie_mat",
	hdrs = ["ie_mat.h"],
	srcs = ["ie_mat.cpp",
		],
	deps = [":vec2",
			":log"],
)


cc_library(
	name = "vec2",
	hdrs = ["vec2.h"],
	srcs = ["vec2.cpp"],
	deps = [":log"],
)


cc_library(
	name = "log",
	hdrs = ["log.h"],
	srcs = ["log.cpp"],
)


cc_library(
	name = "kernel",
	hdrs = ["kernel.h"],
	srcs = ["kernel.cpp"],
	deps = [":ie_mat",
			"//ie-solver/boundaries"]
)


cc_library(
	name = "initialization",
	hdrs = ["initialization.h"],
	srcs = ["initialization.cpp",
			"stokes_init.cpp",
			],
	deps = [":ie_mat",
			":quadtree"]
)


cc_library(
	name = "quadtree",
	hdrs = ["quadtree.h"],
	srcs = ["quadtree.cpp"],
	deps = [":ie_mat",
			"//ie-solver/boundaries"]
)