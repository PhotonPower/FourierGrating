#include "fg/core/prof.hpp"

#include <algorithm>
#include <map>
#include <mutex>
#include <nlohmann/json.hpp>

namespace fg::prof {

namespace detail {
std::atomic<bool> on{false};
}

namespace {
struct Store {
  std::mutex mu;
  std::map<std::string, Entry> entries;
};
Store& store() {
  static Store s;
  return s;
}
}  // namespace

void enable(bool o) {
  detail::on.store(o, std::memory_order_relaxed);
}
bool enabled() {
  return detail::on.load(std::memory_order_relaxed);
}

void reset() {
  Store& s = store();
  const std::lock_guard<std::mutex> lock(s.mu);
  s.entries.clear();
}

void record(const char* name, double seconds) {
  Store& s = store();
  const std::lock_guard<std::mutex> lock(s.mu);
  Entry& e = s.entries[name];
  if (e.name.empty()) e.name = name;
  e.total_s += seconds;
  ++e.count;
  e.max_s = std::max(e.max_s, seconds);
}

std::vector<Entry> snapshot() {
  Store& s = store();
  std::vector<Entry> out;
  {
    const std::lock_guard<std::mutex> lock(s.mu);
    for (const auto& [n, e] : s.entries) out.push_back(e);
  }
  std::stable_sort(out.begin(), out.end(),
                   [](const Entry& a, const Entry& b) { return a.total_s > b.total_s; });
  return out;
}

std::string json() {
  nlohmann::ordered_json j = nlohmann::ordered_json::object();
  for (const Entry& e : snapshot())
    j[e.name] = {{"total_s", e.total_s}, {"count", e.count}, {"max_s", e.max_s}};
  return j.dump(2);
}

}  // namespace fg::prof
