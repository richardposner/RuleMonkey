// RuleMonkey parameter-sweep driver: scan a parameter across a value range
// and report the per-point endpoint response in `.scan` format.
//
// Usage:
//   rm_scan <xml_path> --param <name>
//           (--values v1,v2,... | --min <m> --max <x> --n-points <n> [--log])
//           --t-end <t> [--t-start <t>] [--n-steps <k>] [--seed <s>]
//           [--bifurcate] [--reset-conc 0|1]
//           [--print-functions] [-no-bscb] [--ignore-unsupported]
//
// Output: tab-separated `.scan` to stdout (see docs/scan_format.md).  CPU
// time in seconds to stderr.
//
//   parameter_scan (default): one row per swept value —
//       # <param>  <obs...>  [<func...>]
//   bifurcate (--bifurcate): a forward sweep (min->max) and a backward
//   sweep (max->min) run as one continuous trajectory; each observable /
//   function gets a `_fwd` and `_bwd` column —
//       # <param>  <obs>_fwd  <obs>_bwd ...  [<func>_fwd <func>_bwd ...]
//
// Global-function columns appear only with --print-functions, mirroring
// BNGL's `print_functions=>1` (issue #7).  The endpoint recorded at each
// point is the observable / function value at --t-end, matching BNG's
// parameter_scan, which extracts the last `.gdat` row of each run.

#include "rulemonkey/simulator.hpp"

#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

double parse_double(const char* arg, const char* name) {
  try {
    return std::stod(arg);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("rm_scan: ") + name + " requires a numeric value, got '" +
                             arg + "' (" + e.what() + ")");
  }
}

int parse_int(const char* arg, const char* name) {
  try {
    return std::stoi(arg);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("rm_scan: ") + name + " requires an integer, got '" + arg +
                             "' (" + e.what() + ")");
  }
}

std::uint64_t parse_uint64(const char* arg, const char* name) {
  try {
    return std::stoull(arg);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("rm_scan: ") + name +
                             " requires a non-negative integer, got '" + arg + "' (" + e.what() +
                             ")");
  }
}

// Splits a comma-separated list of doubles for --values.
std::vector<double> parse_value_list(const char* arg) {
  std::vector<double> out;
  const std::string s = arg;
  std::size_t pos = 0;
  while (pos <= s.size()) {
    const std::size_t comma = s.find(',', pos);
    const std::string tok = s.substr(pos, comma - pos);
    if (tok.empty())
      throw std::runtime_error("rm_scan: --values has an empty entry in '" + s + "'");
    out.push_back(parse_double(tok.c_str(), "--values"));
    if (comma == std::string::npos)
      break;
    pos = comma + 1;
  }
  return out;
}

// `argv[i]` needs a following value; advances `i` and returns it.
const char* take_value(int& i, int argc, char* argv[]) {
  if (++i >= argc)
    throw std::runtime_error(std::string("rm_scan: ") + argv[i - 1] + " requires a value");
  return argv[i];
}

void usage() {
  std::cerr << "Usage: rm_scan <xml_path> --param <name>\n"
               "         (--values v1,v2,... | --min <m> --max <x> --n-points <n> [--log])\n"
               "         --t-end <t> [--t-start <t>] [--n-steps <k>] [--seed <s>]\n"
               "         [--bifurcate] [--reset-conc 0|1]\n"
               "         [--print-functions] [-no-bscb] [--ignore-unsupported]\n";
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    usage();
    return 1;
  }

  try {
    const std::string xml_path = argv[1];
    if (!xml_path.empty() && xml_path[0] == '-') {
      usage();
      return 1;
    }

    std::string param;
    std::vector<double> values; // explicit --values list
    bool have_values = false;
    double par_min = 0.0, par_max = 0.0;
    bool have_min = false, have_max = false;
    int n_points = 0;
    bool have_n_points = false;
    bool log_scale = false;
    double t_start = 0.0, t_end = 0.0;
    bool have_t_end = false;
    int n_steps = 1;
    std::uint64_t seed = 42;
    bool bifurcate = false;
    bool reset_conc = true;
    bool bscb = true;
    bool ignore_unsupported = false;
    bool print_functions = false;

    for (int i = 2; i < argc; ++i) {
      if (std::strcmp(argv[i], "--param") == 0) {
        param = take_value(i, argc, argv);
      } else if (std::strcmp(argv[i], "--values") == 0) {
        values = parse_value_list(take_value(i, argc, argv));
        have_values = true;
      } else if (std::strcmp(argv[i], "--min") == 0) {
        par_min = parse_double(take_value(i, argc, argv), "--min");
        have_min = true;
      } else if (std::strcmp(argv[i], "--max") == 0) {
        par_max = parse_double(take_value(i, argc, argv), "--max");
        have_max = true;
      } else if (std::strcmp(argv[i], "--n-points") == 0) {
        n_points = parse_int(take_value(i, argc, argv), "--n-points");
        have_n_points = true;
      } else if (std::strcmp(argv[i], "--log") == 0) {
        log_scale = true;
      } else if (std::strcmp(argv[i], "--t-start") == 0) {
        t_start = parse_double(take_value(i, argc, argv), "--t-start");
      } else if (std::strcmp(argv[i], "--t-end") == 0) {
        t_end = parse_double(take_value(i, argc, argv), "--t-end");
        have_t_end = true;
      } else if (std::strcmp(argv[i], "--n-steps") == 0) {
        n_steps = parse_int(take_value(i, argc, argv), "--n-steps");
      } else if (std::strcmp(argv[i], "--seed") == 0) {
        seed = parse_uint64(take_value(i, argc, argv), "--seed");
      } else if (std::strcmp(argv[i], "--bifurcate") == 0) {
        bifurcate = true;
      } else if (std::strcmp(argv[i], "--reset-conc") == 0) {
        reset_conc = parse_int(take_value(i, argc, argv), "--reset-conc") != 0;
      } else if (std::strcmp(argv[i], "--print-functions") == 0) {
        print_functions = true;
      } else if (std::strcmp(argv[i], "-no-bscb") == 0) {
        bscb = false;
      } else if (std::strcmp(argv[i], "-bscb") == 0) {
        bscb = true;
      } else if (std::strcmp(argv[i], "--ignore-unsupported") == 0) {
        ignore_unsupported = true;
      } else {
        std::cerr << "rm_scan: unknown argument '" << argv[i] << "'\n";
        usage();
        return 1;
      }
    }

    if (param.empty()) {
      std::cerr << "rm_scan: --param is required\n";
      usage();
      return 1;
    }
    if (!have_t_end) {
      std::cerr << "rm_scan: --t-end is required\n";
      usage();
      return 1;
    }
    if (n_steps < 1) {
      std::cerr << "rm_scan: --n-steps must be >= 1 (got " << n_steps << ")\n";
      return 1;
    }
    if (have_values) {
      if (have_min || have_max || have_n_points)
        std::cerr << "WARN: --values given; --min/--max/--n-points ignored\n";
    } else if (!(have_min && have_max && have_n_points)) {
      std::cerr << "rm_scan: provide either --values or all of --min, --max, --n-points\n";
      usage();
      return 1;
    }

    rulemonkey::ScanSpec spec;
    spec.parameter = param;
    spec.values = values;
    spec.par_min = par_min;
    spec.par_max = par_max;
    spec.n_points = n_points;
    spec.log_scale = log_scale;
    spec.reset_conc = reset_conc;
    spec.per_point = rulemonkey::TimeSpec{t_start, t_end, n_steps};

    rulemonkey::RuleMonkeySimulator sim(xml_path);
    sim.set_block_same_complex_binding(bscb);

    bool has_errors = false;
    for (const auto& w : sim.unsupported_features()) {
      const char* level = (w.severity == rulemonkey::Severity::Error) ? "ERROR" : "WARN";
      std::cerr << level << ": Unsupported BNGL feature: " << w.feature
                << " (XML element: " << w.element << ")\n";
      if (w.severity == rulemonkey::Severity::Error)
        has_errors = true;
    }
    if (has_errors && !ignore_unsupported) {
      std::cerr << "ERROR: Model uses unsupported features that produce incorrect results.\n"
                   "       Pass --ignore-unsupported to run anyway.\n";
      return 2;
    }

    const auto t0 = std::chrono::high_resolution_clock::now();

    // Round-trip-precise doubles, matching rm_driver's `.gdat` output.
    std::cout << std::setprecision(17);

    if (bifurcate) {
      const rulemonkey::BifurcateResult r = sim.bifurcate(spec, seed);
      const auto& fwd = r.forward;
      const auto& bwd = r.backward;

      std::cout << "# " << fwd.parameter;
      for (const auto& name : fwd.observable_names)
        std::cout << "\t" << name << "_fwd\t" << name << "_bwd";
      if (print_functions)
        for (const auto& name : fwd.function_names)
          std::cout << "\t" << name << "_fwd\t" << name << "_bwd";
      std::cout << "\n";

      for (std::size_t k = 0; k < fwd.n_points(); ++k) {
        std::cout << fwd.param_values[k];
        for (std::size_t o = 0; o < fwd.observable_names.size(); ++o)
          std::cout << "\t" << fwd.observable_endpoints[k][o] << "\t"
                    << bwd.observable_endpoints[k][o];
        if (print_functions)
          for (std::size_t f = 0; f < fwd.function_names.size(); ++f)
            std::cout << "\t" << fwd.function_endpoints[k][f] << "\t"
                      << bwd.function_endpoints[k][f];
        std::cout << "\n";
      }
    } else {
      const rulemonkey::ScanResult r = sim.parameter_scan(spec, seed);

      std::cout << "# " << r.parameter;
      for (const auto& name : r.observable_names)
        std::cout << "\t" << name;
      if (print_functions)
        for (const auto& name : r.function_names)
          std::cout << "\t" << name;
      std::cout << "\n";

      for (std::size_t k = 0; k < r.n_points(); ++k) {
        std::cout << r.param_values[k];
        for (const double v : r.observable_endpoints[k])
          std::cout << "\t" << v;
        if (print_functions)
          for (const double v : r.function_endpoints[k])
            std::cout << "\t" << v;
        std::cout << "\n";
      }
    }

    const auto t1 = std::chrono::high_resolution_clock::now();
    std::cerr << std::chrono::duration<double>(t1 - t0).count() << "\n";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
}
