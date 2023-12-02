cc_library(
  name = "arithmetic",
  hdrs = ["arithmetic.h"]
)

cc_library(
  name = "polynomial",
  hdrs = ["polynomial.h"],
  deps = [
    ":arithmetic",
    "@com_google_absl//absl/log:check",
  ]
)

cc_binary(
  name = "sandbox",
  srcs = ["sandbox.cc"],
  deps = [
    "@com_google_absl//absl/flags:flag",
    "@com_google_absl//absl/flags:parse",
    ":polynomial",
  ]
)