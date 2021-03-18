#ifndef ENRICO_INCLUDE_ENRICO_TIMER_H
#define ENRICO_INCLUDE_ENRICO_TIMER_H

#include "comm.h"

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

}

#endif // ENRICO_INCLUDE_ENRICO_TIMER_H
