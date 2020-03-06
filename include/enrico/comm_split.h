#ifndef ENRICO_COMM_SPLIT_H
#define ENRICO_COMM_SPLIT_H

#include "comm.h"

#include <array>
#include <mpi.h>

namespace enrico {

//! Splits a given MPI comunicator into new inter- and intra-node communicators
//!
//! This will create two new comms, sub_comm and intranode_comm, for the calling proc:
//! - sub_comm is an internode comm that will span all nodes in the original comm
//! (super_comm).  It
//!   will contain procs_per_node procs in each node.
//! - intranode_comm is an intranode comm that will span a single node.  It will contain
//! the subset
//!   of super_comm that are in the node of the calling proc.
//!
//! In the above description, "node" refers to a shared-memory region as defined by the
//! the MPI implementation's MPI_COMM_TYPE_SHARED.  In most use cases, this will be a
//! node.
//!
//! For the calling proc, the new comms (sub_comm and intranode_comm) will be either valid
//! or invalid depending on these criteria:
//! - If the calling proc is within the new comm, then the comm will be valid for the
//! calling proc.
//! - If the calling proc is not within the new comm, then the comm will be MPI_COMM_NULL
//! for the
//!   calling proc.
//!
//! The calling proc must then make use of with proper value checks.  For example,
//! subsequent MPI operations on sub_comm may need to check if `subComm != MPI_COMM_NULL`.
//!
//! \param[in] super_comm An existing comm that will be split into sub_comm and
//! intranode_comm \param[in] procs_per_node The number of procs per node in the new
//! communicator (sub_comm) \param[out] sub_comm A new internode comm with the desired
//! number of procs per node \param[out] intranode_comm A new intranode comm wil procs in
//! the same node
void get_node_comms(Comm super_comm,
                    int procs_per_node,
                    Comm& sub_comm,
                    Comm& intranode_comm);

void get_driver_comms(Comm super_comm,
                      std::array<int, 2> num_nodes,
                      std::array<int, 2> procs_per_node,
                      std::array<Comm, 2>& driver_comms,
                      Comm& intranode_comm,
                      Comm& coupling_comm);

}

#endif // ENRICO_COMM_SPLIT_H
