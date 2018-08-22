#ifndef STREAM_DRIVERS_H
#define STREAM_DRIVERS_H

#include <map>
#include <unordered_map>
#include <vector>
#include "mpi.h"
#include "openmc/capi.h"
#include "comm.h"
#include "stream_geom.h"
#include "openmc_interface.h"

namespace stream {

constexpr double JOULE_PER_EV = 1.6021766208e-19;

// ============================================================================
// Base Classes
// ============================================================================

class HeatFluidsDriver {
public:
  Comm comm_;

  explicit HeatFluidsDriver(MPI_Comm comm) : comm_(comm) {};
  HeatFluidsDriver() {};
  virtual ~HeatFluidsDriver() {};

  virtual void init_step() {};
  virtual void solve_step() {};
  virtual void finalize_step() {};
  bool active() const;
};

class TransportDriver {
public:
  Comm comm_;

  explicit TransportDriver(MPI_Comm comm) : comm_(comm) {};
  TransportDriver() {};
  virtual ~TransportDriver() {};

  virtual void init_step() {};
  virtual void solve_step() {};
  virtual void finalize_step() {};
  bool active() const;
};

class CoupledDriver {
public:
  Comm comm_;

  TransportDriver transport_driver_;
  HeatFluidsDriver heat_fluids_driver_;

  explicit CoupledDriver(MPI_Comm coupled_comm, MPI_Comm neutron_comm, MPI_Comm heat_fluids_comm);
  CoupledDriver(){};
  virtual ~CoupledDriver() {};
};

// ============================================================================
// Implementations
// ============================================================================

class OpenmcDriver : public TransportDriver {
public:
  // Constructors and destructors
  OpenmcDriver(int argc, char* argv[], MPI_Comm comm);
  ~OpenmcDriver();

  // Methods
  void init_step();
  void solve_step();
  void finalize_step();
  Position get_mat_centroid(int32_t mat_id) const;
  CellInstance* get_cell_instance(int32_t mat_index) const;
  void track_cell_instance(CellInstance c);

  // Data
  int32_t index_tally_;   //!< Index in tallies array for fission tally
  int32_t index_filter_;  //!< Index in filters arrays for material filter
  std::vector<CellInstance> cells_;
private:
  //! Mapping of material indices to cell instances
  std::unordered_map<int32_t, CellInstance*> mat_to_instance_;
};

class NekDriver : public HeatFluidsDriver {
public:
  explicit NekDriver(MPI_Comm comm);
  ~NekDriver();

  void init_step();
  void solve_step();
  void finalize_step();

  Position get_global_elem_centroid(int32_t global_elem) const;

  int lelg_; //!< upper bound on number of mesh elements
  int lelt_; //!< upper bound on number of mesh elements per rank
  int lx1_; //!< polynomial order of the solution
  int nelgt_; //!< total number of mesh elements
  int nelt_; //!< number of local mesh elements
};

// This is not actually derived from CoupledDriver.  Currently, it is unclear
// how or if the base class will be implemented.  The issue will be revisited
class OpenmcNekDriver {
public:
  OpenmcNekDriver(int argc, char *argv[], MPI_Comm coupled_comm, MPI_Comm openmc_comm,
                  MPI_Comm nek_comm, MPI_Comm intranode_comm);
  ~OpenmcNekDriver() {};

  void update_heat_source();
  void update_temperature();

  Comm comm_;
  Comm intranode_comm_;
  OpenmcDriver openmc_driver_;
  NekDriver nek_driver_;
private:
  void init_mappings();
  void init_tallies();

  // Map that gives a list of Nek element global indices for a given OpenMC
  // material index
  std::unordered_map<int32_t,std::vector<int32_t>> mats_to_elems_;
  // Map that gives a list of OpenMC material indices for a given Nek global element index
  std::map<int32_t,int32_t> elems_to_mats_;
  int32_t n_materials_;
};

} // namespace stream

#endif //STREAM_DRIVERS_H
