#ifndef STREAM_DRIVERS_H
#define STREAM_DRIVERS_H

#include <map>
#include <unordered_map>
#include <vector>
#include "mpi.h"
#include "openmc.h"
#include "procinfo.h"
#include "stream_geom.h"

namespace stream {

// ============================================================================
// Base Classes
// ============================================================================

class HeatFluidsDriver {
public:
  ProcInfo procInfo;

  explicit HeatFluidsDriver(MPI_Comm comm) : procInfo(comm) {};
  HeatFluidsDriver() {};
  virtual ~HeatFluidsDriver() {};

  virtual void initStep() {};
  virtual void solveStep() {};
  virtual void finalizeStep() {};
  bool active() const;
};

class TransportDriver {
public:
  ProcInfo procInfo;

  explicit TransportDriver(MPI_Comm comm) : procInfo(comm) {};
  TransportDriver() {};
  virtual ~TransportDriver() {};

  virtual void initStep() {};
  virtual void solveStep() {};
  virtual void finalizeStep() {};
  bool active() const;
};

class CoupledDriver {
public:
  ProcInfo procInfo;

  TransportDriver transportDriver;
  HeatFluidsDriver heatFluidsDriver;

  explicit CoupledDriver(MPI_Comm coupledComm, MPI_Comm neutronComm, MPI_Comm thComm);
  CoupledDriver(){};
  virtual ~CoupledDriver() {};
};

// ============================================================================
// Implementations
// ============================================================================

class OpenmcDriver : public TransportDriver {
public:
  OpenmcDriver(int argc, char* argv[], MPI_Comm comm);
  ~OpenmcDriver();

  void initStep();
  void solveStep();
  void finalizeStep();

  Position getMatCentroid(int32_t matId) const;
  int32_t getMatId(Position position) const;

  int32_t indexTally;
};

class NekDriver : public HeatFluidsDriver {
public:
  explicit NekDriver(MPI_Comm comm);
  ~NekDriver();

  void initStep();
  void solveStep();
  void finalizeStep();

  Position getGlobalElemCentroid(int32_t globalElem) const;

  int lelg;
  int lelt;
  int lx1;
};

// This is not actually derived from CoupledDriver.  Currently, it is unclear
// how or if the base class will be implemented.  The issue will be revisited
class OpenmcNekDriver {
public:
  ProcInfo procInfo;

  OpenmcDriver openmcDriver;
  NekDriver nekDriver;

  OpenmcNekDriver(int argc, char *argv[], MPI_Comm coupledComm, MPI_Comm openmcComm, MPI_Comm nekComm);
  ~OpenmcNekDriver() {};

private:
  void initMatsToElems();
  void initElemsToMats();
  void initTallies();
  // Map that gives a list of Nek element global indices for a given OpenMC
  // material index
  std::unordered_map<int32_t,std::vector<int32_t>> matsToElems;
  // Map that gives a list of OpenMC material indices for a given Nek global element index
  std::map<int32_t, std::vector<int32_t>> elemsToMats;
};

} // namespace stream

#endif //STREAM_DRIVERS_H
