package(default_visibility = ["//visibility:public"])

cc_library(
	name = "boundaries",
	hdrs = ["circle.h",
			"annulus.h",
			"rounded_square.h",
			"rounded_square_with_bump.h",
			"squiggly.h",
			"cubic_spline.h",
			"ex1boundary.h",
			"ex3boundary.h",
			"boundary.h"],
	srcs = ["circle.cpp",
			"annulus.cpp",
			"rounded_square.cpp",
			"rounded_square_with_bump.cpp",
			"squiggly.cpp",
			"cubic_spline.cpp",
			"cubic_boundary.cpp",
			"ex1boundary.cpp",
			"ex3boundary.cpp",
			],
	deps = ["//ie_solver:ie_mat",
			"//ie_solver:vec2"]
)
