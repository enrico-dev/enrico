#ifndef STREAM_DRIVERS_H
#define STREAM_DRIVERS_H

#include "mpi.h"
#include "openmc.h"

// ============================================================================
// Base Classes
// ============================================================================

class ThDriver {
public:
  MPI_Comm comm;

  ThDriver() {};
  ThDriver(MPI_Comm comm) : comm(comm) {};
  virtual ~ThDriver() {};

  virtual void initStep();
  virtual void solveStep();
  virtual void finalizeStep();

};

class NeutronDriver {
public:
  MPI_Comm comm;

  NeutronDriver() {};
  NeutronDriver(MPI_Comm comm) : comm(comm) {};
  virtual ~NeutronDriver() {};

  virtual void initStep();
  virtual void solveStep();
  virtual void finalizeStep();
};

// ============================================================================
// Implementations
// ============================================================================

class OpenmcDriver : public NeutronDriver {
public:
  MPI_Comm comm;

  OpenmcDriver(MPI_Comm comm);
  ~OpenmcDriver();

  void initStep();
  void solveStep();
  void finalizeStep();

  int getCellCentroid(int32_t cellId, double *centroid);
};

class NekDriver : public ThDriver {
public:
  MPI_Comm comm;

  NekDriver(MPI_Comm comm);
  ~NekDriver();

  void initStep();
  void solveStep();
  void finalizeStep();

  int getElemNumber();
  int getElemCentroid(int globalElem);
};

class CoupledDriver {
public:
  MPI_Comm globalComm;

  NeutronDriver neutronDriver;
  ThDriver thDriver;

  CoupledDriver(){};
  CoupledDriver(MPI_Comm globalComm, MPI_Comm neutronComm, MPI_Comm thComm);
  ~CoupledDriver() {};
};

#endif //STREAM_DRIVERS_H
