// examples/embed.cpp — minimal C++ embedding of librulemonkey.
//
// Reads a BNGL XML file (produced by BNG2.pl writeXML), runs a
// network-free simulation, and prints the final-time observable values.
//
// Build (from the repository root):
//   cmake --preset release -DRULEMONKEY_BUILD_EXAMPLES=ON
//   cmake --build --preset release --target rm_embed_example
//
// Run:
//   build/release/examples/rm_embed_example path/to/model.xml
//
// Hosting from another CMake project:
//   add_subdirectory(path/to/RuleMonkey)
//   add_executable(my_app my_app.cpp)
//   target_link_libraries(my_app PRIVATE rulemonkey)
//
// or, if RuleMonkey has been installed:
//   find_package(RuleMonkey CONFIG REQUIRED)
//   target_link_libraries(my_app PRIVATE RuleMonkey::rulemonkey)

#include <rulemonkey/simulator.hpp>

#include <cstdio>
#include <exception>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <model.xml> [t_end] [n_points] [seed]\n", argv[0]);
    return 1;
  }

  const std::string xml_path = argv[1];
  const double t_end = (argc > 2) ? std::stod(argv[2]) : 10.0;
  const int n_points = (argc > 3) ? std::stoi(argv[3]) : 11;
  const std::uint64_t seed = (argc > 4) ? std::stoull(argv[4]) : 42;

  try {
    rulemonkey::RuleMonkeySimulator sim(xml_path);

    // Strict BNGL semantics: bimolecular rules only fire across complexes
    // (matches NFsim's -bscb flag).  Set false to allow intramolecular binding.
    sim.set_block_same_complex_binding(true);

    // Surface any features the loaded model uses but RM does not implement.
    for (const auto& w : sim.unsupported_features()) {
      const char* level = (w.severity == rulemonkey::Severity::Error) ? "ERROR" : "WARN";
      std::fprintf(stderr, "%s: %s (XML element: %s)\n", level, w.feature.c_str(),
                   w.element.c_str());
    }

    // One-shot run: stateless, repeatable for a given seed.
    rulemonkey::TimeSpec ts{/*t_start=*/0.0, t_end, n_points};
    auto result = sim.run(ts, seed);

    // Print final time-point as a sanity-check trajectory summary.
    // Global-function columns (if the model has a `begin functions`
    // block) follow the observables, same column-major layout.
    std::cout << "# time";
    for (auto& name : result.observable_names)
      std::cout << "\t" << name;
    for (auto& name : result.function_names)
      std::cout << "\t" << name;
    std::cout << "\n";

    if (result.n_times() > 0) {
      const std::size_t last = result.n_times() - 1;
      std::cout << result.time[last];
      for (std::size_t obs = 0; obs < result.n_observables(); ++obs)
        std::cout << "\t" << result.observable_data[obs][last];
      for (std::size_t fn = 0; fn < result.n_functions(); ++fn)
        std::cout << "\t" << result.function_data[fn][last];
      std::cout << "\n";
    }
    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 2;
  }
}
