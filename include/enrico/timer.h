#ifndef ENRICO_INCLUDE_ENRICO_TIMER_H
#define ENRICO_INCLUDE_ENRICO_TIMER_H

#include "comm.h"
#include <iomanip>
#include <string>

namespace enrico {

//! Class for measuring and collecting time on a given MPI communicator
class Timer {
public:
  //! Initializes timer for a given MPI communicator
  //! \param comm The MPI communicator for which time is measured
  explicit Timer(const Comm& comm)
    : comm_(comm)
  {}

  //! Begin accumulating elapsed time.
  //!
  //! If the timer has been previously stopped, then restarting it will
  //! not reset elapsed time at 0.
  void start();

  //! Stop accumulating elapsed time.
  //!
  //! This does not reset elapsed time to 0.
  void stop();

  //! Accumulated time between all consecutive calls to start() and stop()
  //!
  //! If the timer is currently running, then this returns the elapsed time
  //! at the time of the function call.
  //!
  //! If the timer is not currently running (i.e., has been stopped), then this returns
  //! the accumulated elapsed time from the time it was stopped.
  //!
  //! \return Elapsed time in seconds.
  double elapsed();

  //! Reset the elapsed time to 0.
  void reset();

private:
  const Comm comm_;      //!< MPI comm for which this instance measures time
  double start_ = 0.0;   //!< Start time at most recent call to start()
  double elapsed_ = 0.0; //!< Time accumulated between all consecutive start/stops
  bool running_ = false; //!< True if started; false if stopped
};

//! Class for storing times associated with an arbitrary label
class TimeAmt {
public:
  explicit TimeAmt(const std::string& name)
    : name(name){};
  TimeAmt(const std::string& name, double time)
    : name(name)
    , time(time){};
  TimeAmt(const std::string& name, double time, double percent)
    : name(name)
    , time(time)
    , percent(percent){};

  //! Get the total time for a vector of TimeAmt
  static double sum_times(const std::vector<TimeAmt>& times);

  //! Get the total percent for a vector of TimeAmt
  static double sum_percent(const std::vector<TimeAmt>& times);

  //! Print the times for a vector TimeAmt in a nicely-formatted way
  //!
  //! \param header_name An arbitrary header that is printed
  //! \param times Each TimeAmt is printed on a separate line
  //! \param comm The root process of this Comm will print the output
  static void print_times(const std::string& header_name,
                          const std::vector<TimeAmt>& times,
                          const Comm& comm);

  const std::string name; //!< Arbitrary label
  double time;            //!< The time in arbitrary units
  double percent;         //!< The percent wrt. a total time
};

}

#endif // ENRICO_INCLUDE_ENRICO_TIMER_H
