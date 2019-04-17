package(default_visibility = ["//visibility:public"])

cc_library(
	name = "boundaries",
	hdrs = ["circle.h",
			"annulus.h",
			"rounded_square.h",
			"rounded_square_with_bump.h",
			"squiggly.h",
			"ellipses.h",
			"cubic_spline.h",
			"boundary.h"],
	srcs = ["circle.cpp",
			"annulus.cpp",
			"rounded_square.cpp",
			"rounded_square_with_bump.cpp",
			"squiggly.cpp",
			"ellipses.cpp",
			"cubic_spline.cpp",
			],
	deps = ["//ie-solver:ie_mat",
			"//ie-solver:vec2"]
)
