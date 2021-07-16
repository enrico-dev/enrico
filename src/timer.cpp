//
// Created by Ronald Rahaman on 3/18/21.
//

#include "enrico/timer.h"

namespace enrico {

void Timer::start()
{
  if (comm_.active()) {
    running_ = true;
    comm_.Barrier();
    start_ = MPI_Wtime();
  }
}

void Timer::stop()
{
  if (comm_.active()) {
    elapsed_ = elapsed();
    running_ = false;
  }
}

void Timer::reset()
{
  if (comm_.active()) {
    running_ = false;
    elapsed_ = 0.0;
  }
}

double Timer::elapsed()
{
  if (comm_.active()) {
    if (running_) {
      comm_.Barrier();
      auto diff = MPI_Wtime() - start_;
      return elapsed_ + diff;
    } else {
      return elapsed_;
    }
  } else {
    return 0.0;
  }
}

double TimeAmt::sum_times(const std::vector<TimeAmt>& times)
{
  double tot = 0.0;
  for (const auto& t : times) {
    tot += t.time;
  }
  return tot;
}

double TimeAmt::sum_percent(const std::vector<TimeAmt>& times)
{
  double tot = 0.0;
  for (const auto& t : times) {
    tot += t.percent;
  }
  return tot;
}

void TimeAmt::print_times(const std::string& header_name,
                          const std::vector<TimeAmt>& times,
                          const Comm& comm)
{
  std::ios_base::fmtflags old_flags(std::cout.flags());
  std::stringstream msg;
  msg << "  " << header_name << " time (seconds, percent)";
  comm.message(msg.str());
  for (const auto& t : times) {
    std::stringstream msg;
    msg << "    " << std::setw(22) << std::left << t.name << std::right << std::scientific
        << std::setprecision(4) << t.time << "    " << std::setw(8) << std::fixed
        << std::left << std::right << t.percent * 100.0;
    comm.message(msg.str());
  }
  std::cout.flags(old_flags);
}

}
