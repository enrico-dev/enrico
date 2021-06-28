#ifndef ENRICO_INCLUDE_ENRICO_TIMER_H
#define ENRICO_INCLUDE_ENRICO_TIMER_H

#include "comm.h"
#include <string>
#include <iomanip>

namespace enrico {

class Timer {
public:
  explicit Timer(const Comm& comm)
    : comm_(comm)
  {}

  void start();
  void stop();
  double elapsed();
  void reset();

private:
  const Comm comm_;
  double start_ = 0.0;
  double stop_ = 0.0;
  double elapsed_ = 0.0;
  bool running_ = false;
};

class TimeAmt {
  public:
    explicit TimeAmt(const std::string& name) : name(name) {};
    TimeAmt(const std::string& name, double time) : name(name), time(time) {};
    TimeAmt(const std::string& name, double time, double percent) : name(name), time(time), percent(percent) {};
    const std::string name;
    double time;
    double percent;
};

double sum_times(const std::vector<TimeAmt>& times);
void print_times(const std::string &header_name, const std::vector<TimeAmt>& times, const Comm& comm);

}

#endif // ENRICO_INCLUDE_ENRICO_TIMER_H
