#include "enrico/comm_split.h"

namespace enrico {

void get_node_comms(Comm super_comm,
                    int procs_per_node,
                    Comm& sub_comm,
                    Comm& intranode_comm)
{
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

}
