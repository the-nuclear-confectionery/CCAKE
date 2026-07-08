#ifndef CCAKE_PHASE_PROFILER_H
#define CCAKE_PHASE_PROFILER_H

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <map>
#include <ostream>
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>

namespace ccake {

/// @class PhaseProfiler
/// @brief Lightweight, opt-in per-phase wall-time profiler for the SPH step loop.
/// @details
/// Two independent instrumentation channels, both added by @ref ScopedPhase:
///  - A Kokkos profiling region (`pushRegion`/`popRegion`) is ALWAYS emitted.
///    These are no-ops unless a Kokkos tool is attached (e.g. via
///    `KOKKOS_TOOLS_LIBS=.../libkp_kernel_timer.so`), so external profilers see
///    the named phases with no measurable overhead otherwise.
///  - Built-in wall-time accumulation per phase, active ONLY when the
///    environment variable `CCAKE_PROFILE` is set to a non-empty, non-"0" value.
///    This lets a normal run print a per-phase breakdown without any external
///    tooling.
/// When `CCAKE_PROFILE` is unset the only added work is the (no-op) region
/// push/pop, so default runs are unchanged. Timing relies on each instrumented
/// phase completing its own Kokkos work before the scope closes (true on the
/// OpenMP backend, and wherever the phase ends in a fence), so NO extra fences
/// are introduced and the measured numbers are accurate on the CPU path.
class PhaseProfiler {
public:
  static PhaseProfiler& instance() {
    static PhaseProfiler p;
    return p;
  }

  bool enabled() const { return enabled_; }

  /// Step interval for interim reports (env CCAKE_PROFILE_EVERY); 0 = end only.
  int report_every() const { return report_every_; }

  /// Accumulate @p seconds against phase @p name (first call registers it).
  void add(const char* name, double seconds) {
    auto it = index_.find(name);
    if (it == index_.end()) {
      index_.emplace(name, entries_.size());
      entries_.push_back(Entry{name, seconds, 1});
    } else {
      Entry& e = entries_[it->second];
      e.total += seconds;
      e.count += 1;
    }
    total_ += seconds;
  }

  /// Print a summary table (sorted by total time, descending). No-op unless
  /// enabled and at least one phase was recorded.
  void report(std::ostream& os) const {
    if (!enabled_ || entries_.empty()) return;
    std::vector<const Entry*> rows;
    rows.reserve(entries_.size());
    for (const auto& e : entries_) rows.push_back(&e);
    std::sort(rows.begin(), rows.end(),
              [](const Entry* a, const Entry* b) { return a->total > b->total; });

    os << "\n  Per-phase wall-time profile (CCAKE_PROFILE)\n";
    os << "  " << std::string(70, '-') << "\n";
    os << "  " << std::left << std::setw(34) << "phase"
       << std::right << std::setw(11) << "total[s]"
       << std::setw(9) << "calls"
       << std::setw(10) << "ms/call"
       << std::setw(7) << "%" << "\n";
    os << "  " << std::string(70, '-') << "\n";
    for (const Entry* e : rows) {
      const double ms = e->count ? 1.0e3 * e->total / static_cast<double>(e->count) : 0.0;
      const double pct = total_ > 0.0 ? 100.0 * e->total / total_ : 0.0;
      os << "  " << std::left << std::setw(34) << e->name
         << std::right << std::fixed << std::setprecision(4) << std::setw(11) << e->total
         << std::setw(9) << e->count
         << std::setprecision(3) << std::setw(10) << ms
         << std::setprecision(1) << std::setw(6) << pct << "%" << "\n";
    }
    os << "  " << std::string(70, '-') << "\n";
    os << "  " << std::left << std::setw(34) << "sum(profiled)"
       << std::right << std::fixed << std::setprecision(4) << std::setw(11) << total_
       << "\n\n";
  }

private:
  struct Entry {
    std::string name;
    double total;
    long count;
  };

  PhaseProfiler() {
    const char* env = std::getenv("CCAKE_PROFILE");
    enabled_ = (env != nullptr && std::strlen(env) > 0 && std::strcmp(env, "0") != 0);
    const char* ev = std::getenv("CCAKE_PROFILE_EVERY");
    report_every_ = (enabled_ && ev != nullptr) ? std::atoi(ev) : 0;
    if (report_every_ < 0) report_every_ = 0;
  }

  bool enabled_ = false;
  int report_every_ = 0;
  double total_ = 0.0;
  std::vector<Entry> entries_;
  std::map<std::string, std::size_t> index_;
};

/// @class ScopedPhase
/// @brief RAII helper: opens a Kokkos profiling region and (when enabled)
/// accumulates the scope's wall time into the singleton @ref PhaseProfiler.
class ScopedPhase {
public:
  explicit ScopedPhase(const char* name) : name_(name) {
    Kokkos::Profiling::pushRegion(name_);
    if (PhaseProfiler::instance().enabled())
      t0_ = std::chrono::high_resolution_clock::now();
  }
  ~ScopedPhase() {
    if (PhaseProfiler::instance().enabled()) {
      const auto t1 = std::chrono::high_resolution_clock::now();
      PhaseProfiler::instance().add(
          name_, std::chrono::duration<double>(t1 - t0_).count());
    }
    Kokkos::Profiling::popRegion();
  }
  ScopedPhase(const ScopedPhase&) = delete;
  ScopedPhase& operator=(const ScopedPhase&) = delete;

private:
  const char* name_;
  std::chrono::high_resolution_clock::time_point t0_;
};

} // namespace ccake

#define CCAKE_PHASE_CONCAT_(a, b) a##b
#define CCAKE_PHASE_CONCAT(a, b) CCAKE_PHASE_CONCAT_(a, b)
/// Time the enclosing scope under @p name (string literal). Pair with braces:
/// `{ CCAKE_PHASE("smooth_fields"); smooth_all_particle_fields(t2); }`
#define CCAKE_PHASE(name) \
  ccake::ScopedPhase CCAKE_PHASE_CONCAT(ccake_phase_, __LINE__)(name)

#endif // CCAKE_PHASE_PROFILER_H
