cc_binary(
	name = "laplace-ie-solver",
	srcs = ["ex1.cpp"],
	deps = [":common",
			":initialization",
			":skelfac",
			":quadtree",
			":circle"
			],
	linkopts = ["-llapacke",
				"-llapack",
				"-lblas",
				"-lgomp"]

)

cc_binary(
	name = "stokes-ie-solver",
	srcs = ["ex2.cpp"],
	deps = [":common",
			":initialization",
			":skelfac",
			":quadtree",
			":circle"
			],
	linkopts = ["-llapacke",
				"-llapack",
				"-lblas",
				"-lgomp"]

)

cc_library(
	name = "common",
	hdrs = ["common.h"],
	srcs = ["common.cpp"],

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
	name = "skelfac",
	hdrs = ["skelfac.h"],
	srcs = ["skelfac.cpp"],
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