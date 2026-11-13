#pragma once

// Issue #1235 — numbers read and written in the "C" locale, whatever the
// process's LC_NUMERIC.
//
// strtod, std::stod, atof and snprintf follow the C locale's decimal point. In a
// process that has called setlocale(LC_ALL, "") on a German or French machine,
// "0.5" read as 0 and 0.5 was written as "0,5". Python leaves LC_NUMERIC at "C"
// unless something calls setlocale, and its own float() reads a point
// whatever the locale, so only this side changed: a homodimer's stat factor
// `0.5*kf1` was refused as an unknown parameter, a group factor `0.5*1,2`
// read 0 with no error, and an inline tfun's 2.5 read 2. ExprTk parses its
// own numbers and is not affected.
//
// Each helper takes the same arguments and gives the same answer as the call
// it replaces, read or written in "C". Where a "C" locale cannot be made (it
// always can on the platforms bngsim builds for), they fall back to the
// process's own.
//
// Code that is not bngsim's to change reads numbers the same way, in-process:
// RuleMonkey's loader (std::stod: a TFUN's 0.5 read 0, with no error), NFsim's
// TFUN parser (strtod: refused), MIR's C front end (strtod: every decimal
// literal refused) and libantimony. ScopedCNumericLocale puts "C" in force
// around each call into them, for the calling thread only.

#include <cerrno>
#include <clocale>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <locale.h>
#include <stdexcept>
#include <string>
#if defined(__APPLE__) || defined(__FreeBSD__)
#include <xlocale.h>
#endif

#if defined(__GNUC__) || defined(__clang__)
#define BNGSIM_PRINTF_FORMAT(fmt, args) __attribute__((format(printf, fmt, args)))
#else
#define BNGSIM_PRINTF_FORMAT(fmt, args)
#endif

namespace bngsim {

namespace detail {
#if defined(_WIN32)
inline _locale_t c_numeric_locale() {
    static const _locale_t loc = _create_locale(LC_ALL, "C");
    return loc;
}
#else
inline locale_t c_numeric_locale() {
    static const locale_t loc = newlocale(LC_ALL_MASK, "C", static_cast<locale_t>(0));
    return loc;
}
#endif
} // namespace detail

// strtod, read in "C".
inline double c_strtod(const char *s, char **end) {
    const auto loc = detail::c_numeric_locale();
    if (!loc)
        return std::strtod(s, end);
#if defined(_WIN32)
    return _strtod_l(s, end, loc);
#else
    return strtod_l(s, end, loc);
#endif
}

// atof, read in "C".
inline double c_atof(const char *s) { return c_strtod(s, nullptr); }

// std::stod, read in "C": the same leading whitespace, the same count of
// characters read in *pos, std::invalid_argument where nothing converts and
// std::out_of_range where strtod reports ERANGE (as libstdc++ and libc++ do).
inline double c_stod(const std::string &s, std::size_t *pos = nullptr) {
    const char *begin = s.c_str();
    char *end = nullptr;
    const int saved = errno;
    errno = 0;
    const double v = c_strtod(begin, &end);
    const int got = errno;
    errno = saved;
    if (end == begin)
        throw std::invalid_argument("stod");
    if (got == ERANGE)
        throw std::out_of_range("stod");
    if (pos)
        *pos = static_cast<std::size_t>(end - begin);
    return v;
}

// snprintf, written in "C". The locale is switched for this thread only, for
// the length of the call. It returns what vsnprintf does, the length the whole
// output needs, on Windows too: _vsnprintf_s_l alone gives -1 when it truncates.
BNGSIM_PRINTF_FORMAT(3, 4)
inline int c_snprintf(char *buf, std::size_t n, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    const auto loc = detail::c_numeric_locale();
    int r;
#if defined(_WIN32)
    if (loc) {
        va_list len_ap;
        va_copy(len_ap, ap);
        r = _vscprintf_l(fmt, loc, len_ap);
        va_end(len_ap);
        if (n > 0)
            _vsnprintf_s_l(buf, n, _TRUNCATE, fmt, loc, ap);
    } else {
        r = std::vsnprintf(buf, n, fmt, ap);
    }
#else
    if (loc) {
        const locale_t old = uselocale(loc);
        r = std::vsnprintf(buf, n, fmt, ap);
        uselocale(old);
    } else {
        r = std::vsnprintf(buf, n, fmt, ap);
    }
#endif
    va_end(ap);
    return r;
}

// "C" in force for LC_NUMERIC on the calling thread while the object lives, and
// the thread's previous setting put back after: uselocale on POSIX, and on
// Windows setlocale under _configthreadlocale(_ENABLE_PER_THREAD_LOCALE). The
// process's own setting is never touched, so no other thread reads in "C", and
// none has its own switch undone when this one ends, as with a process-wide
// setlocale. Nests.
class ScopedCNumericLocale {
  public:
    ScopedCNumericLocale() {
#if defined(_WIN32)
        prev_mode_ = _configthreadlocale(_ENABLE_PER_THREAD_LOCALE);
        if (const char *cur = std::setlocale(LC_NUMERIC, nullptr))
            prev_numeric_ = cur;
        std::setlocale(LC_NUMERIC, "C");
#else
        if (const auto loc = detail::c_numeric_locale())
            prev_ = uselocale(loc);
#endif
    }
    ~ScopedCNumericLocale() {
#if defined(_WIN32)
        if (!prev_numeric_.empty())
            std::setlocale(LC_NUMERIC, prev_numeric_.c_str());
        if (prev_mode_ != -1)
            _configthreadlocale(prev_mode_);
#else
        if (prev_)
            uselocale(prev_);
#endif
    }
    ScopedCNumericLocale(const ScopedCNumericLocale &) = delete;
    ScopedCNumericLocale &operator=(const ScopedCNumericLocale &) = delete;

  private:
#if defined(_WIN32)
    int prev_mode_ = -1;
    std::string prev_numeric_;
#else
    locale_t prev_ = static_cast<locale_t>(0);
#endif
};

} // namespace bngsim

#undef BNGSIM_PRINTF_FORMAT
