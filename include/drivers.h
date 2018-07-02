#ifndef STREAM_DRIVERS_H
#define STREAM_DRIVERS_H

#include <map>
#include <unordered_map>
#include <vector>
#include "mpi.h"
#include "openmc.h"
#include "procinfo.h"
#include "stream_geom.h"

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

class CoupledDriver {
public:
  ProcInfo procInfo;

  NeutronDriver neutronDriver;
  ThDriver thDriver;

  explicit CoupledDriver(MPI_Comm coupledComm, MPI_Comm neutronComm, MPI_Comm thComm);
  CoupledDriver(){};
  virtual ~CoupledDriver() {};
};

// ============================================================================
// Implementations
// ============================================================================

class OpenmcDriver : public NeutronDriver {
// ROR 2018-06-02: Declaring this member function as a friend isn't working for
// me.  I must have an error that I can't see.
// friend void OpenmcNekDriver::initMatsToNekElems();
friend class OpenmcNekDriver;
public:
  explicit OpenmcDriver(int argc, char* argv[], MPI_Comm comm);
  ~OpenmcDriver();

  void initStep();
  void solveStep();
  void finalizeStep();

  Position getMatCentroid(const int32_t matId);
  int32_t getMatId(const Position position);
};

class NekDriver : public ThDriver {
// ROR 2018-06-02: Declaring a member function (rather than whole class) as a friend
// isn't working for me.  I must have an error that I can't see.
friend class OpenmcNekDriver;
public:
  explicit NekDriver(MPI_Comm comm);
  ~NekDriver();

  void initStep();
  void solveStep();
  void finalizeStep();

  Position getGlobalElemCentroid(const int32_t globalElem);

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

  explicit OpenmcNekDriver(int argc, char *argv[], MPI_Comm coupledComm, MPI_Comm openmcComm, MPI_Comm nekComm);
  ~OpenmcNekDriver() {};

private:
  void initMatsToElems();
  void initElemsToMats();
  // Map that gives a list of Nek element global indices for a given OpenMC
  // material index
  std::unordered_map<int32_t,std::vector<int32_t>> matsToElems;
  // Map that gives a list of OpenMC material indices for a given Nek global element index
  std::map<int32_t, std::vector<int32_t>> elemsToMats;
};

#endif //STREAM_DRIVERS_H
