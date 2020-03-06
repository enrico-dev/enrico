#include "enrico/comm_split.h"
#include <initializer_list>

namespace enrico {

void get_node_comms(Comm super_comm,
                    int procs_per_node,
                    Comm& sub_comm,
                    Comm& intranode_comm)
{
  // Each intranode_comm contains all ranks in a given shared memory region (usually node)
  MPI_Comm icomm;
  MPI_Comm_split_type(
    super_comm.comm, MPI_COMM_TYPE_SHARED, super_comm.rank, MPI_INFO_NULL, &icomm);
  intranode_comm = Comm(icomm);

  // Split the specified number of procs_per_node from the super_comm
  // We only want the comm where color == 0.  The second comm is destroyed.
  int color = intranode_comm.rank < procs_per_node ? 0 : 1;
  MPI_Comm scomm;
  MPI_Comm_split(super_comm.comm, color, super_comm.rank, &scomm);
  sub_comm = Comm(scomm);
  if (color != 0) {
    sub_comm.free();
  }
}

void get_driver_comms(Comm super_comm,
                      std::array<int, 2> num_nodes,
                      std::array<int, 2> procs_per_node,
                      std::array<Comm, 2>& driver_comms,
                      Comm& intranode_comm,
                      Comm& coupling_comm)
{
  const int KEEP = 0;
  const int DISCARD = 1;
  MPI_Comm temp_comm;

  // Each intranode_comm contains all ranks in a given shared memory region (usually node)
  MPI_Comm_split_type(
    super_comm.comm, MPI_COMM_TYPE_SHARED, super_comm.rank, MPI_INFO_NULL, &temp_comm);
  intranode_comm = Comm(temp_comm);

  // The coupling comm consists of the root process on each intranode_comm.
  int color = intranode_comm.is_root() ? KEEP : DISCARD;
  MPI_Comm_split(super_comm.comm, color, super_comm.rank, &temp_comm);
  coupling_comm = Comm(temp_comm);
  if (color == DISCARD) {
    coupling_comm.free();
  }

  // Get the number of nodes (inferred as the size of the coupling comm)
  int total_nodes = coupling_comm.size;
  intranode_comm.broadcast(total_nodes);

  // Each rank gets its node index (inferred from the ranks of the coupling comms)
  int node_idx;
  if (coupling_comm.active()) {
    std::vector<int> v;
    if (coupling_comm.is_root()) {
      for (int i = 0; i < total_nodes; ++i) {
        v.push_back(i);
      }
    }
    MPI_Scatter(v.data(), 1, MPI_INT, &node_idx, 1, MPI_INT, 0, coupling_comm.comm);
  }
  intranode_comm.broadcast(node_idx);

  // Get the driver comms. driver_comms[0] gets the left-hand nodes, and
  // driver_comms[1] gets the right-hand nodes, both based on the node_idx
  for (const int i : {0, 1}) {

    // 0 is the null value for num_nodes and procs_per_node.  In that case, we use the
    // maximum number of nodes or procs-per-node
    auto n = num_nodes[i] > 0 ? num_nodes[i] : total_nodes;
    auto ppn = procs_per_node[i] > 0 ? procs_per_node[i] : intranode_comm.size;
    auto& scomm = driver_comms[i];

    int color;
    // Left-hand nodes
    if (i == 0) {
      color = (node_idx < n && intranode_comm.rank < ppn) ? KEEP : DISCARD;
    }
    // Right-hand nodes
    else {
      color = (node_idx >= total_nodes - n && intranode_comm.rank < ppn) ? KEEP : DISCARD;
    }
    MPI_Comm_split(super_comm.comm, color, super_comm.rank, &temp_comm);
    scomm = Comm(temp_comm);
    if (color == DISCARD) {
      scomm.free();
    }
  }
}

}
