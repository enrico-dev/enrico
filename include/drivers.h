#ifndef STREAM_DRIVERS_H
#define STREAM_DRIVERS_H

#include "mpi.h"
#include "omc_interface.h

class SubDriver
{
public:
  MPI_Comm comm;

  SubDriver(MPI_Comm comm) : comm(comm) {};
  virtual ~SubDriver() {};

  // Different from constructor?
  virtual void init();
  virtual void initStep();
  virtual void solve();
  virtual void closeStep();
  // Different from destructor?
  virtual void clean();
};

class NekDriver : public SubDriver
{
  MPI_Comm comm;

  NekDriver(MPI_Comm comm) : comm(comm) {};
  ~NekDriver() {};

  // Different from constructor?
  void init() {};
  void initStep() {};
  void solve() {};
  void closeStep() {};
  // Different from destructor?
  void clean() {};

  // Changes interface of SubDriver.  Maybe we don't need ABC?
  int getElemNumber();
  // Return type?
  double getElemCentroid(int elemId);
};


class OmcDriver : public SubDriver
{
public:
  MPI_Comm comm;

  OmcDriver(MPI_Comm comm) : comm(comm) {};
  ~OmcDriver() {};

  // Different from constructor?
  void init() {};
  void initStep() {};
  void solve() {};
  void closeStep() {};
  // Different from destructor?
  void clean() {};

  // Changes interface of SubDriver.  Maybe subdriver doesn't need to be ABC?
  double coord


};

#endif //STREAM_DRIVERS_H
