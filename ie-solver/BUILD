cc_binary(
	name = "ie-solver",
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
				"-lgfortran",
				"-lgomp",
				"-lpthread"]

)

cc_library(
	name = "common",
	hdrs = ["common.h"],
	srcs = ["common.cpp"],

)

cc_library(
	name = "initialization",
	hdrs = ["initialization.h"],
	srcs = ["initialization.cpp"],
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