// RuleMonkey batch driver: run a model from XML and output .gdat format.
// Usage: rm_driver <xml_path> <t_end> <n_steps> [seed] [-no-bscb]
//        [--ignore-unsupported] [--save-state <path>] [--load-state <path>]
//        [--t-start <time>] [--print-functions] [--species <path>]
//        [--reps <n> [--threads <n>] [--sd-out <path>]]
// Output: tab-separated .gdat to stdout (# header, then time +
//         observables; global-function columns appended only with
//         --print-functions)
// Stderr: CPU time in seconds
//
// RuleMonkey defaults to strict BNGL semantics (block_same_complex_binding=true,
// matching NFsim -bscb).  Pass -no-bscb to disable the checks for testing
// against NFsim runs that were done without -bscb/-cb.
//
// Pass --ignore-unsupported to suppress errors for unsupported BNGL features
// in the XML (warnings are always printed).
//
// Pass --print-functions to append the model's global-function values as
// trailing .gdat columns (after the observables).  Off by default, mirroring
// BNGL's `print_functions=>1` opt-in; the in-process Result API exposes the
// values unconditionally regardless of this flag.
//
// State continuation:
//   --save-state <path>  Save full simulation state at t_end to file.
//   --load-state <path>  Load state from file instead of seed species.
//   --t-start <time>     Start time (default 0).  From seed species, the
//                        seed state sits at <time> on the model's clock and
//                        the run spans t_end - <time> (issue #98).  With
//                        --load-state it must match the saved time, which
//                        it defaults to.
//
// Species output:
//   --species <path>     After the run, write the final-state species
//                        census to <path> as a BNG-format `.species`
//                        file (BNG2.pl `readNFspecies`-compatible).
//
// Replicate ensembles (issue #113):
//   --reps <n>           Run n replicates seeded seed, seed+1, ..., and write
//                        their per-time-point MEAN as the .gdat on stdout.
//                        Replicate i is the run `rm_driver ... <seed+i>`
//                        would print.
//   --threads <n>        Worker threads (default 0 = one per hardware
//                        thread).  Changes speed only, never the output.
//   --sd-out <path>      Also write the sample standard deviation (n - 1)
//                        to <path>, in the same .gdat layout as stdout.
//   Not combinable with --save-state / --load-state / --species, which
//   describe one trajectory's final state.

#include "rulemonkey/simulator.hpp"

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

double parse_double(const char* arg, const char* name) {
  try {
    return std::stod(arg);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("rm_driver: ") + name +
                             " requires a numeric value, got '" + arg + "' (" + e.what() + ")");
  }
}

int parse_int(const char* arg, const char* name) {
  try {
    return std::stoi(arg);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("rm_driver: ") + name + " requires an integer, got '" +
                             arg + "' (" + e.what() + ")");
  }
}

std::uint64_t parse_uint64(const char* arg, const char* name) {
  try {
    return std::stoull(arg);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("rm_driver: ") + name +
                             " requires a non-negative integer, got '" + arg + "' (" + e.what() +
                             ")");
  }
}

// Writes one .gdat table: a `# time <columns>` header, then one row per
// sample.  `obs` / `fns` are column-major ([column][time]), as in Result.
// Global-function columns (the model's `begin functions` entries with no
// local arguments) follow the observables, but only when
// `print_functions` is set.  This mirrors BNGL's `print_functions=>1`
// opt-in: the default `.gdat` is observables-only.
void write_gdat(std::ostream& out, const std::vector<double>& time,
                const std::vector<std::string>& obs_names,
                const std::vector<std::vector<double>>& obs,
                const std::vector<std::string>& fn_names,
                const std::vector<std::vector<double>>& fns, bool print_functions) {
  out << "#";
  out << " time";
  for (auto& name : obs_names)
    out << "\t" << name;
  if (print_functions)
    for (auto& name : fn_names)
      out << "\t" << name;
  out << "\n";

  // Round-trip-precise doubles in the .gdat output: 17 decimal digits
  // is the IEEE-754 binary64 round-trip threshold.  Default ostream
  // precision is 6 (silently lossy) — downstream parity checks that
  // diff RM .gdat against NFsim .gdat would otherwise read fuzzed
  // values and could hide drift below the 1e-6 mark.
  out << std::setprecision(17);
  for (size_t t = 0; t < time.size(); ++t) {
    out << time[t];
    for (const auto& col : obs)
      out << "\t" << col[t];
    if (print_functions)
      for (const auto& col : fns)
        out << "\t" << col[t];
    out << "\n";
  }
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::cerr << "Usage: rm_driver <xml_path> <t_end> <n_steps> [seed] [-no-bscb]"
                 " [--ignore-unsupported] [--save-state <path>]"
                 " [--load-state <path>] [--t-start <time>] [--print-functions]"
                 " [--species <path>] [--reps <n> [--threads <n>] [--sd-out <path>]]\n";
    return 1;
  }

  try {
    std::string const xml_path = argv[1];
    double const t_end = parse_double(argv[2], "t_end");
    int const n_steps = parse_int(argv[3], "n_steps");
    if (n_steps < 1) {
      std::cerr << "ERROR: n_steps must be >= 1 (got " << n_steps
                << ").  Output is `n_steps + 1` sample rows from t_start to t_end;\n"
                   "       use a positive integer.  See docs/gdat_format.md.\n";
      return 1;
    }
    uint64_t seed = 42;
    bool bscb = true; // strict BNGL semantics by default
    bool ignore_unsupported = false;
    bool print_functions = false; // opt-in; matches BNGL `print_functions=>1`
    bool seed_set = false;
    std::string save_state_path;
    std::string load_state_path;
    std::string species_path;
    double t_start = 0.0;
    bool t_start_set = false;
    int reps = 0; // 0 = single trajectory (no --reps)
    int threads = 0;
    bool threads_set = false;
    std::string sd_out_path;

    // Parse remaining args: seed (positional, may appear at most once) and
    // flags (named).  Reject unknown `--*` / `-*` flags up front so a typo
    // (e.g. `--seed 7`) is not silently parsed as the seed positional via
    // `std::stoull` and crash with an opaque "stoull: no conversion".
    for (int i = 4; i < argc; ++i) {
      if (std::strcmp(argv[i], "-no-bscb") == 0) {
        bscb = false;
      } else if (std::strcmp(argv[i], "-bscb") == 0) {
        bscb = true;
      } else if (std::strcmp(argv[i], "--ignore-unsupported") == 0) {
        ignore_unsupported = true;
      } else if (std::strcmp(argv[i], "--print-functions") == 0) {
        print_functions = true;
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
      } else if (std::strcmp(argv[i], "--species") == 0) {
        if (++i >= argc) {
          std::cerr << "--species requires a path\n";
          return 1;
        }
        species_path = argv[i];
      } else if (std::strcmp(argv[i], "--reps") == 0) {
        if (++i >= argc) {
          std::cerr << "--reps requires a value\n";
          return 1;
        }
        reps = parse_int(argv[i], "--reps");
        if (reps < 1) {
          std::cerr << "rm_driver: --reps must be >= 1 (got " << reps << ")\n";
          return 1;
        }
      } else if (std::strcmp(argv[i], "--threads") == 0) {
        if (++i >= argc) {
          std::cerr << "--threads requires a value\n";
          return 1;
        }
        threads = parse_int(argv[i], "--threads");
        threads_set = true;
        if (threads < 0) {
          std::cerr << "rm_driver: --threads must be >= 0 (got " << threads << ")\n";
          return 1;
        }
      } else if (std::strcmp(argv[i], "--sd-out") == 0) {
        if (++i >= argc) {
          std::cerr << "--sd-out requires a path\n";
          return 1;
        }
        sd_out_path = argv[i];
      } else if (std::strcmp(argv[i], "--t-start") == 0) {
        if (++i >= argc) {
          std::cerr << "--t-start requires a value\n";
          return 1;
        }
        t_start = parse_double(argv[i], "--t-start");
        t_start_set = true;
      } else if (argv[i][0] == '-') {
        std::cerr << "rm_driver: unknown flag '" << argv[i]
                  << "'.  See usage above for the supported set.\n";
        return 1;
      } else if (!seed_set) {
        seed = parse_uint64(argv[i], "seed");
        seed_set = true;
      } else {
        std::cerr << "rm_driver: unexpected positional argument '" << argv[i]
                  << "' (seed already provided).\n";
        return 1;
      }
    }

    if (reps == 0 && (threads_set || !sd_out_path.empty())) {
      std::cerr << "rm_driver: --threads and --sd-out require --reps\n";
      return 1;
    }
    if (reps > 0 &&
        (!save_state_path.empty() || !load_state_path.empty() || !species_path.empty())) {
      std::cerr << "rm_driver: --reps cannot be combined with --save-state, --load-state or"
                   " --species, which describe a single trajectory\n";
      return 1;
    }

    rulemonkey::RuleMonkeySimulator sim(xml_path);
    sim.set_block_same_complex_binding(bscb);

    // Check for unsupported BNGL features
    bool has_errors = false;
    for (auto& w : sim.unsupported_features()) {
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

    // Opened before the ensemble runs, so a bad path fails at once instead
    // of after every replicate has been simulated.
    std::ofstream sd_out;
    if (!sd_out_path.empty()) {
      sd_out.open(sd_out_path);
      if (!sd_out)
        throw std::runtime_error("cannot open --sd-out path '" + sd_out_path + "'");
    }

    auto t0 = std::chrono::high_resolution_clock::now();

    if (reps > 0) {
      rulemonkey::TimeSpec ts;
      ts.t_start = t_start;
      ts.t_end = t_end;
      ts.n_points = n_steps;
      rulemonkey::EnsembleSpec spec;
      spec.n_reps = reps;
      spec.n_threads = threads;
      auto ens = sim.run_ensemble(ts, spec, seed);
      double const cpu_s =
          std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();

      write_gdat(std::cout, ens.time, ens.observable_names, ens.observable_mean, ens.function_names,
                 ens.function_mean, print_functions);
      if (sd_out.is_open()) {
        write_gdat(sd_out, ens.time, ens.observable_names, ens.observable_sd, ens.function_names,
                   ens.function_sd, print_functions);
        if (!sd_out.flush())
          throw std::runtime_error("failed writing --sd-out path '" + sd_out_path + "'");
      }
      std::cerr << cpu_s << "\n";
      return 0;
    }

    rulemonkey::Result result;
    // --species needs the pool to survive past the run, so it forces
    // session mode (the one-shot run() path keeps no inspectable pool).
    bool const use_session =
        !save_state_path.empty() || !load_state_path.empty() || !species_path.empty();

    if (use_session) {
      // Session mode: needed for save/load state
      if (!load_state_path.empty()) {
        sim.load_state(load_state_path);
        if (!t_start_set) {
          // Default t_start to the loaded session time so simulate()'s
          // t_start ≈ current_time invariant holds.  Without this, a
          // user who forgot --t-start would hit the contract throw.
          t_start = sim.current_time();
          std::cerr << "WARN: --load-state without --t-start; using t_start=" << t_start
                    << " from loaded state\n";
        }
      } else {
        sim.initialize(seed, t_start);
      }

      result = sim.simulate(t_start, t_end, n_steps);

      if (!save_state_path.empty())
        sim.save_state(save_state_path);

      if (!species_path.empty())
        sim.write_species_file(species_path);

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
    double const cpu_s = std::chrono::duration<double>(t1 - t0).count();

    write_gdat(std::cout, result.time, result.observable_names, result.observable_data,
               result.function_names, result.function_data, print_functions);

    // CPU time to stderr
    std::cerr << cpu_s << "\n";

    return 0;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
}
