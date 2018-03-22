#ifndef STREAM_DRIVERS_H
#define STREAM_DRIVERS_H

#include "mpi.h"
#include "openmc.h"

class SubDriver {
public:
  MPI_Comm comm;

  SubDriver() {};
  SubDriver(MPI_Comm comm) : comm(comm) {};
  virtual ~SubDriver() {};

  // ROR: 2018-03-22: Different from constructor/destructor?
  void initDriver();
  void freeDriver();

  // Different from constructor?
  virtual void initStep();
  virtual void solveStep();
  virtual void finalizeStep();
};


class OpenmcDriver : public SubDriver {
public:
  MPI_Comm comm;

  OpenmcDriver(MPI_Comm comm);
  ~OpenmcDriver() {};

  // ROR: 2018-03-22: Different from constructor/destructor?
  void initDriver();
  void freeDriver();

  void initStep();
  void solveStep();
  void finalizeStep();

  // Prefer this?
  //   int getCellCentroid(int32_t cellId, double centroid[3]);
  int getCellCentroid(int32_t cellId, double *centroid);
};


class NekDriver : public SubDriver {
public:
  MPI_Comm comm;

  NekDriver(MPI_Comm comm);
  ~NekDriver() {};

  // ROR: 2018-03-22: Different from constructor/destructor?
  void initDriver();
  void freeDriver();

  void initStep();
  void solveStep();
  void finalizeStep();

  int getElemNumber();
  int getElemCentroid(int globalElem);
};


class CoupledDriver {
public:
  MPI_Comm globalComm;
  OpenmcDriver openmc;
  NekDriver nek;

  CoupledDriver(MPI_Comm globalComm, MPI_Comm openmcComm, MPI_Comm nekComm);
  ~CoupledDriver() {};
};

#endif //STREAM_DRIVERS_H
