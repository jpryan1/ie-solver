package(default_visibility = ["//visibility:public"])

cc_binary(
	name = "ie_solver_dev",
	srcs = ["ie_solver_dev.cpp"],
	deps = [
			":ie_mat",
			":initialization",
			"//ie_solver/skel_factorization",
			"//ie_solver/quadtree",
			"//ie_solver/boundaries",
			":linear_solve",
			],
	linkopts = ["-llapacke",
				"-llapack",
				"-lblas",
				"-lgomp"]

)


cc_library(
	name = "linear_solve",
	hdrs = ["linear_solve.h"],
	srcs = ["linear_solve.cpp"],
	deps = [":ie_mat",
			":initialization",
			"//ie_solver/skel_factorization",
			"//ie_solver/quadtree",
			":log",
			"//ie_solver/io",
			"//ie_solver/boundaries"],
  linkopts = ["-pthread"],

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
	name = "initialization",
	hdrs = ["initialization.h"],
	srcs = ["initialization.cpp"],
	copts= ["-fopenmp"],
	deps = [":ie_mat",
			"//ie_solver/quadtree",
			"//ie_solver/kernel"]
)


cc_test(
    name = "ie_solver_test",
    srcs = ["ie_solver_test.cpp"],
	deps = [
			":ie_mat",
			"//ie_solver/skel_factorization",
			"//ie_solver/quadtree",
			"//ie_solver/boundaries",
			":linear_solve",
			"@gtest//:gtest",
      "@gtest//:gtest_main"
    ],
  	linkopts = ["-llapacke",
			"-llapack",
			"-lblas",
			"-lgomp"]
)