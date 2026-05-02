// Quick smoke test: parse an NFsim XML file and run a short simulation.
#include "rulemonkey/simulator.hpp"
#include <cstdlib>
#include <iostream>
#include <string>

namespace {
[[noreturn]] void fail(const std::string& msg) {
  std::cerr << "FAIL: " << msg << "\n";
  std::exit(1);
}
} // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: smoke_test <xml_path>\n";
    return 1;
  }

  try {
    std::string xml_path = argv[1];
    std::cout << "Loading model from: " << xml_path << "\n";

    rulemonkey::RuleMonkeySimulator sim(xml_path);

    auto obs_names = sim.observable_names();
    if (obs_names.empty())
      fail("model exposes zero observables — smoke test cannot verify trajectory shape");

    std::cout << "Observable names:\n";
    for (auto& name : obs_names)
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

    // n_points samples + the t_start row → n_points+1 expected.  An
    // off-by-one in the recorder or a silent abort mid-run would
    // change this; without the assert the test would still print
    // results and exit 0.
    constexpr std::size_t expected_n_times = 11;
    if (result.n_times() != expected_n_times)
      fail("expected " + std::to_string(expected_n_times) +
           " sample rows (n_points + t_start), got " + std::to_string(result.n_times()));
    if (result.n_observables() != obs_names.size())
      fail("observable count mismatch between Simulator metadata and Result rows");

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

    std::cout << "\nDone.\n";
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
}
