#include <iostream>
#include <string>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

ABSL_FLAG(std::string, echo, "", "Text to echo");

int main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);

  std::cout << absl::GetFlag(FLAGS_echo) << std::endl;

  return 0;
}
