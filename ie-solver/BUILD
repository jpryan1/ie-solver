cc_binary(
	name = "ie-solver",
	srcs = ["ie-solver.cpp"],
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
	srcs = ["vec2.h", 
			"log.h", 
			"box.h", 
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