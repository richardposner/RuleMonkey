// Quick smoke test: parse an NFsim XML file and run a short simulation.
#include "rulemonkey/simulator.hpp"
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: smoke_test <xml_path>\n";
    return 1;
  }

  try {
    std::string xml_path = argv[1];
    std::cout << "Loading model from: " << xml_path << "\n";

    rulemonkey::RuleMonkeySimulator sim(xml_path);

    std::cout << "Observable names:\n";
    for (auto& name : sim.observable_names())
      std::cout << "  " << name << "\n";

    std::cout << "Parameter names:\n";
    for (auto& name : sim.parameter_names())
      std::cout << "  " << name << "\n";

    rulemonkey::TimeSpec ts;
    ts.t_start = 0.0;
    ts.t_end = 10.0;
    ts.n_points = 10;

    std::cout << "\nRunning simulation (t=0..10, 10 points)...\n";
    auto result = sim.run(ts, 12345);

    std::cout << "Events fired: " << result.event_count << "\n";
    std::cout << "Time points: " << result.n_times() << "\n";

    if (result.n_times() > 0 && result.n_observables() > 0) {
      std::cout << "\nResults (time";
      for (auto& name : result.observable_names)
        std::cout << "\t" << name;
      std::cout << "):\n";

      for (size_t t = 0; t < result.n_times(); ++t) {
        std::cout << result.time[t];
        for (size_t o = 0; o < result.n_observables(); ++o)
          std::cout << "\t" << result.observable_data[o][t];
        std::cout << "\n";
      }
    }

    std::cout << "\nDone.\n";
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
}
