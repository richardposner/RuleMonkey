// The XML loader reads numbers the way BNG2 writes them, whatever locale the
// host process has set (issue #122).
//
// strtod and std::stod follow the decimal point of the C locale, so under a
// decimal-comma LC_NUMERIC (de_DE, fr_FR, ...) they stopped at the `.` of
// "0.5", returned 0, and reported nothing.  A host that embeds RuleMonkey can
// set one: a GUI, or a Python host that calls setlocale(LC_ALL, "").  The
// rule's symmetry_factor and a TFUN's inline xData/yData went through stod;
// a TFUN read from a file went through an iostream, which follows the C++
// global locale instead.
//
// The oracle is RuleMonkey itself in the "C" locale: under a decimal-comma
// locale, set either way, every model must give the identical trajectory.
//
//   locale_sym_model.xml   ring closure, symmetry_factor="0.5"
//   tfun_time_rate_model   inline TFUN, rewritten here with fractional data
//   ft_tfun.xml            TFUN from rates.tfun, copied with a fractional x
//
// It also pins the other half of the fix: a value that is only partly a
// number is refused at load, where stod read its leading number and dropped
// the rest.
//
// The locale half needs a decimal-comma locale installed.  Without one it is
// skipped (exit 77); the refusal half runs regardless.

#include "rulemonkey/simulator.hpp"

#include <clocale>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

int g_failures = 0;

void check(bool ok, const std::string& msg) {
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s\n", msg.c_str());
    ++g_failures;
  }
}

std::string read_file(const fs::path& p) {
  std::ifstream const in(p, std::ios::binary);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

void write_file(const fs::path& p, const std::string& text) {
  std::ofstream out(p, std::ios::binary);
  out << text;
}

// `text` with its one occurrence of `from` replaced by `to`.
std::string replace_once(std::string text, const std::string& from, const std::string& to) {
  auto const at = text.find(from);
  if (at == std::string::npos || text.find(from, at + 1) != std::string::npos)
    throw std::runtime_error("fixture does not contain exactly one '" + from + "'");
  return text.replace(at, from.size(), to);
}

std::vector<std::vector<double>> trajectory(const fs::path& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml.string());
  return sim.run(rulemonkey::TimeSpec{0.0, 10.0, 3, {}}, 1).observable_data;
}

// Whether loading `xml` throws a runtime_error that names `token`.
bool refuses(const fs::path& xml, const std::string& token) {
  try {
    rulemonkey::RuleMonkeySimulator const sim(xml.string());
  } catch (const std::runtime_error& e) {
    return std::string(e.what()).find(token) != std::string::npos;
  }
  return false;
}

// Set the process's C locale for `category`; false if `name` is not
// installed.  setlocale and localeconv are not thread safe, and this test is
// single-threaded: changing the process locale is the thing under test.
bool set_c_locale(int category, const char* name) {
  return std::setlocale(category, name) != nullptr; // NOLINT(concurrency-mt-unsafe)
}
char c_decimal_point() {
  return std::localeconv()->decimal_point[0]; // NOLINT(concurrency-mt-unsafe)
}

// A locale whose decimal point is a comma, by the name setlocale knows it by,
// or "" if none is installed.
std::string decimal_comma_locale() {
  for (const char* name : {"de_DE.UTF-8", "de_DE.utf8", "fr_FR.UTF-8", "fr_FR.utf8", "de_DE",
                           "fr_FR", "de-DE", "fr-FR", "German_Germany.1252"}) {
    if (!set_c_locale(LC_NUMERIC, name))
      continue;
    bool const comma = c_decimal_point() == ',';
    set_c_locale(LC_NUMERIC, "C");
    if (comma)
      return name;
  }
  return "";
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::fprintf(stderr, "Usage: locale_number_test <locale_sym_model.xml> "
                         "<tfun_time_rate_model.xml> <ft_tfun.xml>\n");
    return 2;
  }
  fs::path const sym_xml = argv[1];
  fs::path const tfun_src = argv[2];
  fs::path const tfun_file_xml = argv[3];

  fs::path const dir = fs::temp_directory_path() / "rulemonkey_locale_number_test";
  fs::remove_all(dir);
  fs::create_directories(dir);

  bool skipped = false;
  try {
    // Inline TFUN data with a fractional point on each axis.
    std::string const tfun_text = read_file(tfun_src);
    fs::path const tfun_inline_xml = dir / "tfun_inline.xml";
    write_file(tfun_inline_xml,
               replace_once(replace_once(tfun_text, "xData=\"0,2\"", "xData=\"0,2.5\""),
                            "yData=\"0,200\"", "yData=\"0.5,200.25\""));

    // --- A value that is only partly a number is refused, naming it. ---
    std::string const sym_text = read_file(sym_xml);
    fs::path const bad_sym = dir / "bad_sym.xml";
    write_file(bad_sym,
               replace_once(sym_text, "symmetry_factor=\"0.5\"", "symmetry_factor=\"0.5x\""));
    check(refuses(bad_sym, "0.5x"), "symmetry_factor=\"0.5x\" must be refused at load");

    fs::path const bad_inline = dir / "bad_inline.xml";
    write_file(bad_inline, replace_once(tfun_text, "yData=\"0,200\"", "yData=\"0,200x\""));
    check(refuses(bad_inline, "200x"), "TFUN yData \"0,200x\" must be refused at load");

    fs::path const bad_dir = dir / "bad_file";
    fs::create_directories(bad_dir);
    write_file(bad_dir / "ft_tfun.xml", read_file(tfun_file_xml));
    write_file(bad_dir / "rates.tfun",
               replace_once(read_file(tfun_file_xml.parent_path() / "rates.tfun"), "2.0", "2.0x"));
    check(refuses(bad_dir / "ft_tfun.xml", "2.0x"),
          "a .tfun value \"2.0x\" must be refused at load");

    // A TFUN read from a file, with a fractional value in each column:
    // rates.tfun's x column is whole numbers.
    fs::path const frac_dir = dir / "frac_file";
    fs::create_directories(frac_dir);
    fs::path const tfun_frac_xml = frac_dir / "ft_tfun.xml";
    write_file(tfun_frac_xml, read_file(tfun_file_xml));
    write_file(
        frac_dir / "rates.tfun",
        replace_once(read_file(tfun_file_xml.parent_path() / "rates.tfun"), "\n20 ", "\n20.5 "));

    // --- The same trajectory under a decimal-comma locale, set either way. ---
    std::vector<fs::path> const models = {sym_xml, tfun_inline_xml, tfun_frac_xml};
    std::vector<std::vector<std::vector<double>>> reference;
    reference.reserve(models.size());
    for (const auto& m : models)
      reference.push_back(trajectory(m));

    std::string const comma = decimal_comma_locale();
    if (comma.empty()) {
      std::fprintf(stdout, "locale_number_test: no decimal-comma locale installed, "
                           "locale half skipped\n");
      skipped = true;
    } else {
      auto same_as_reference = [&](const std::string& how) {
        for (size_t i = 0; i < models.size(); ++i) {
          std::string const what = models[i].filename().string() + " under " + how;
          try {
            check(trajectory(models[i]) == reference[i], what + " differs from the C locale");
          } catch (const std::exception& e) {
            check(false, what + " failed: " + e.what());
          }
        }
      };

      // The C locale only, as setlocale(LC_ALL, "") in a host does.
      set_c_locale(LC_ALL, comma.c_str());
      same_as_reference("setlocale(LC_ALL, \"" + comma + "\")");
      set_c_locale(LC_ALL, "C");

      // The C++ global locale, which iostreams read and which also sets the
      // C locale when it has a name.
      bool have_cxx_locale = true;
      try {
        std::locale::global(std::locale(comma));
      } catch (const std::runtime_error&) {
        have_cxx_locale = false;
        std::fprintf(stdout,
                     "locale_number_test: std::locale(\"%s\") unavailable, that "
                     "case skipped\n",
                     comma.c_str());
      }
      if (have_cxx_locale)
        same_as_reference("std::locale::global(\"" + comma + "\")");
      std::locale::global(std::locale::classic());
      set_c_locale(LC_ALL, "C");
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "FAIL: exception: %s\n", e.what());
    ++g_failures;
  }
  fs::remove_all(dir);

  if (g_failures != 0)
    return 1;
  std::fprintf(stdout, "locale_number_test: OK\n");
  return skipped ? 77 : 0;
}
