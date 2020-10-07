#include "enrico/comm_split.h"

namespace enrico {

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
  int node_idx = coupling_comm.rank;
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

std::vector<int> gather_subcomm_ranks(const Comm& super, const Comm& sub)
{
  std::vector<int> ranks(super.size);
  super.Allgather(&sub.rank, 1, MPI_INT, ranks.data(), 1, MPI_INT);
  auto new_end =
    std::remove_if(ranks.begin(), ranks.end(), [](int i) { return i == MPI_PROC_NULL; });
  ranks.erase(new_end, ranks.end());
  assert(sub.size == ranks.size());
  return ranks;
}

}
