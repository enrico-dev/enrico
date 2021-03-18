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

}
