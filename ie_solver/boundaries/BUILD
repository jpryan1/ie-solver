package(default_visibility = ["//visibility:public"])

cc_library(
	name = "boundaries",
	hdrs = [
			"circle.h",
			"annulus.h",
			"donut.h",
			"rounded_square.h",
#			"rounded_square_with_bump.h",
#			"squiggly.h",
			"cubic_spline.h",
			"ex1boundary.h",
			"ex2boundary.h",
			"ex3boundary.h",
			"ex4boundary.h",
			"boundary.h"],
	srcs = [
			"circle.cpp",
			"annulus.cpp",
			"donut.cpp",
			"rounded_square.cpp",
#			"rounded_square_with_bump.cpp",
#			"squiggly.cpp",
			"cubic_spline.cpp",
			"boundary.cpp",
			"cubic_boundary.cpp",
			"ex1boundary.cpp",
			"ex2boundary.cpp",
			"ex3boundary.cpp",
			"ex4boundary.cpp",
			],
	deps = ["//ie_solver:ie_mat",
			"//ie_solver:vec2"]
)
