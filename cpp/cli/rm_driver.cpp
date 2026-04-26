// RuleMonkey batch driver: run a model from XML and output .gdat format.
// Usage: rm_driver <xml_path> <t_end> <n_steps> [seed] [-no-bscb]
//        [--ignore-unsupported] [--save-state <path>] [--load-state <path>]
//        [--t-start <time>]
// Output: tab-separated .gdat to stdout (# header, then time + observables)
// Stderr: CPU time in seconds
//
// RuleMonkey defaults to strict BNGL semantics (block_same_complex_binding=true,
// matching NFsim -bscb).  Pass -no-bscb to disable the checks for testing
// against NFsim runs that were done without -bscb/-cb.
//
// Pass --ignore-unsupported to suppress errors for unsupported BNGL features
// in the XML (warnings are always printed).
//
// State continuation:
//   --save-state <path>  Save full simulation state at t_end to file.
//   --load-state <path>  Load state from file instead of seed species.
//   --t-start <time>     Override t_start (default 0; use with --load-state).

#include "rulemonkey/simulator.hpp"

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::cerr << "Usage: rm_driver <xml_path> <t_end> <n_steps> [seed] [-no-bscb]"
                 " [--ignore-unsupported] [--save-state <path>]"
                 " [--load-state <path>] [--t-start <time>]\n";
    return 1;
  }

  std::string xml_path = argv[1];
  double t_end = std::stod(argv[2]);
  int n_steps = std::stoi(argv[3]);
  uint64_t seed = 42;
  bool bscb = true; // strict BNGL semantics by default
  bool ignore_unsupported = false;
  std::string save_state_path;
  std::string load_state_path;
  double t_start = 0.0;
  bool t_start_set = false;

  // Parse remaining args: seed (positional) and flags
  for (int i = 4; i < argc; ++i) {
    if (std::strcmp(argv[i], "-no-bscb") == 0) {
      bscb = false;
    } else if (std::strcmp(argv[i], "-bscb") == 0) {
      bscb = true;
    } else if (std::strcmp(argv[i], "--ignore-unsupported") == 0) {
      ignore_unsupported = true;
    } else if (std::strcmp(argv[i], "--save-state") == 0) {
      if (++i >= argc) {
        std::cerr << "--save-state requires a path\n";
        return 1;
      }
      save_state_path = argv[i];
    } else if (std::strcmp(argv[i], "--load-state") == 0) {
      if (++i >= argc) {
        std::cerr << "--load-state requires a path\n";
        return 1;
      }
      load_state_path = argv[i];
    } else if (std::strcmp(argv[i], "--t-start") == 0) {
      if (++i >= argc) {
        std::cerr << "--t-start requires a value\n";
        return 1;
      }
      t_start = std::stod(argv[i]);
      t_start_set = true;
    } else {
      seed = std::stoull(argv[i]);
    }
  }

  try {
    rulemonkey::RuleMonkeySimulator sim(xml_path);
    sim.set_block_same_complex_binding(bscb);

    // Check for unsupported BNGL features
    bool has_errors = false;
    for (auto& w : sim.unsupported_warnings()) {
      const char* level = (w.severity == rulemonkey::Severity::Error) ? "ERROR" : "WARN";
      std::cerr << level << ": Unsupported BNGL feature: " << w.feature
                << " (XML element: " << w.element << ")\n";
      if (w.severity == rulemonkey::Severity::Error)
        has_errors = true;
    }
    if (has_errors && !ignore_unsupported) {
      std::cerr << "ERROR: Model uses unsupported features that produce incorrect results.\n"
                << "       Pass --ignore-unsupported to run anyway.\n";
      return 2;
    }

    auto t0 = std::chrono::high_resolution_clock::now();

    rulemonkey::Result result;
    bool use_session = !save_state_path.empty() || !load_state_path.empty();

    if (use_session) {
      // Session mode: needed for save/load state
      if (!load_state_path.empty()) {
        sim.load_state(load_state_path);
        if (!t_start_set)
          std::cerr << "WARN: --load-state without --t-start; using t_start=0\n";
      } else {
        sim.initialize(seed);
      }

      result = sim.simulate(t_start, t_end, n_steps);

      if (!save_state_path.empty())
        sim.save_state(save_state_path);

      sim.destroy_session();
    } else {
      // One-shot mode (fastest path)
      rulemonkey::TimeSpec ts;
      ts.t_start = t_start;
      ts.t_end = t_end;
      ts.n_points = n_steps;
      result = sim.run(ts, seed);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double cpu_s = std::chrono::duration<double>(t1 - t0).count();

    // Output header
    std::cout << "#";
    std::cout << " time";
    for (auto& name : result.observable_names)
      std::cout << "\t" << name;
    std::cout << "\n";

    // Output data
    for (size_t t = 0; t < result.n_times(); ++t) {
      std::cout << result.time[t];
      for (size_t o = 0; o < result.n_observables(); ++o)
        std::cout << "\t" << result.observable_data[o][t];
      std::cout << "\n";
    }

    // CPU time to stderr
    std::cerr << cpu_s << "\n";

    return 0;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
}
