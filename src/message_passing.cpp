#include "stream/message_passing.h"

#include <mpi.h>

namespace stream {

void get_node_comms(MPI_Comm super_comm, int procs_per_node, MPI_Comm* sub_comm,
                    MPI_Comm* intranode_comm)
{

  // super_comm_rank is used as the "key" to retain ordering in the comm splits.
  // This can allow the sub_comm to retain some intent from the super_comm's proc layout
  int super_comm_rank;
  MPI_Comm_rank(super_comm, &super_comm_rank);

  // intranode_comm is an intermediate object.  It is only used to get an intranode_comm_rank,
  // which is used as the "color" in the final comm split.
  MPI_Comm_split_type(super_comm,
                      MPI_COMM_TYPE_SHARED,
                      super_comm_rank,
                      MPI_INFO_NULL,
                      intranode_comm);
  int intranode_comm_rank;
  MPI_Comm_rank(*intranode_comm, &intranode_comm_rank);

  // Finally, split the specified number of procs_per_node from the super_comm
  // We only want the comm where color == 0.  The second comm is destroyed.
  int color = intranode_comm_rank < procs_per_node ? 0 : 1;
  MPI_Comm_split(super_comm, color, super_comm_rank, sub_comm);
  if (color != 0)
    MPI_Comm_free(sub_comm);
}

} // namespace stream
