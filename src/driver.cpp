#include "enrico/driver.h"
#include <iomanip>

namespace enrico {

bool Driver::active() const
{
  return comm_.active();
}

void Driver::timer_report()
{
  std::map<std::string, double> times{{"init_step", timer_init_step.elapsed()},
                                      {"solve_step", timer_solve_step.elapsed()},
                                      {"write_step", timer_write_step.elapsed()},
                                      {"finalize_step", timer_finalize_step.elapsed()}};

  for (const auto& t : times) {
    std::stringstream msg;
    msg << "    " << std::setw(24) << std::left << t.first << std::right
        << std::scientific << t.second;
    comm_.message(msg.str());
  }
}

} // namespace enrico
