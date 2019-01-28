package(default_visibility = ["//visibility:public"])

cc_library(
	name = "boundaries",
	hdrs = ["circle.h",
			"rounded_square.h",
			"rounded_square_with_bump.h",
			"squiggly.h",
			"boundary.h"],
	srcs = ["circle.cpp", 
			"rounded_square.cpp", 
			"rounded_square_with_bump.cpp",
			"squiggly.cpp"],
	deps = ["//ie-solver:ie_mat", 
			"//ie-solver:vec2"]
)
