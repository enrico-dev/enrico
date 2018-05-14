#ifndef STREAM_DRIVERS_H
#define STREAM_DRIVERS_H

#include <unordered_map>
#include <vector>
#include "mpi.h"
#include "openmc.h"
#include "procinfo.h"

// ============================================================================
// Base Classes
// ============================================================================

class ThDriver {
public:
  ProcInfo procInfo;

  explicit ThDriver(MPI_Comm comm) : procInfo(comm) {};
  ThDriver() {};
  virtual ~ThDriver() {};

  virtual void initStep() {};
  virtual void solveStep() {};
  virtual void finalizeStep() {};
};

class NeutronDriver {
public:
  ProcInfo procInfo;

  explicit NeutronDriver(MPI_Comm comm) : procInfo(comm) {};
  NeutronDriver() {};
  virtual ~NeutronDriver() {};

  virtual void initStep() {};
  virtual void solveStep() {};
  virtual void finalizeStep() {};
};

// ============================================================================
// Implementations
// ============================================================================

class OpenmcDriver : public NeutronDriver {
public:
  explicit OpenmcDriver(int argc, char* argv[], MPI_Comm comm);
  ~OpenmcDriver();

  void initStep();
  void solveStep();
  void finalizeStep();

  int getCellCentroid(int32_t cellId, double *centroid);
private:
  // Map that gives a list of Nek element global indices for a given OpenMC
  // material index
  std::unordered_map<int32_t,std::vector<int>> mats_to_nek_elems;
};

class NekDriver : public ThDriver {
public:
  explicit NekDriver(MPI_Comm comm);
  ~NekDriver();

  void initStep();
  void solveStep();
  void finalizeStep();

  int getElemNumber();
  int getElemCentroid(int globalElem);
};

class CoupledDriver {
public:
  ProcInfo procInfo;

  NeutronDriver neutronDriver;
  ThDriver thDriver;

  explicit CoupledDriver(MPI_Comm coupledComm, MPI_Comm neutronComm, MPI_Comm thComm) :
      procInfo(coupledComm),
      neutronDriver(neutronComm),
      thDriver(thComm) {};
  CoupledDriver(){};
  virtual ~CoupledDriver() {};
};

#endif //STREAM_DRIVERS_H
