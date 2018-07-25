#include "message_passing.h"
#include "mpi.h"

namespace stream {

/**
 * @brief Splits a given MPI comunicator into a new comm with a specified number of procs in every node.
 *
 * The new comm (*subComm*) will span **all** nodes in the original comm (*superComm*).  This will work correctly even
 * if *subComm* spans only a single node.
 *
 * Depending on the calling proc, *subComm* will be one of two values:
 *   - If the calling proc is within the new desired comm (the comm with the given number of *procsPerNode*), then
 *     *subComm* will be a valid comm.
 *   - If the calling proc is not within the new desired comm, then *subComm* will be `MPI_COMM_NULL`.
 * The caller must then make use of *subComm* with proper value checks.  For example, subsequent MPI operations on
 * *subComm* may need to check if `subComm != MPI_COMM_NULL`.
 *
 * @param super_comm An existing communicator that will be split
 * @param procs_per_node The number of MPI procs per node in the new communicator, *subComm*
 * @param sub_comm A new communicator with either the given number of procs per node *or* `MPI_COMM_NULL`, depending on
 *                whether the calling proc is in the desired comm.
 */
void get_internode_sub_comm(MPI_Comm super_comm, int procs_per_node, MPI_Comm *sub_comm) {

  // super_comm_rank is used as the "key" to retain ordering in the comm splits.
  // This can allow the sub_comm to retain some intent from the super_comm's proc layout
  int super_comm_rank = MPI_PROC_NULL;
  MPI_Comm_rank(super_comm, &super_comm_rank);

  // intranode_comm is an intermediate object.  It is only used to get an intranode_comm_rank,
  // which is used as the "color" in the final comm split.
  MPI_Comm intranode_comm = MPI_COMM_NULL;
  MPI_Comm_split_type(super_comm, MPI_COMM_TYPE_SHARED, super_comm_rank, MPI_INFO_NULL, &intranode_comm);
  int intranode_comm_rank = MPI_PROC_NULL;
  MPI_Comm_rank(intranode_comm, &intranode_comm_rank);
  MPI_Comm_free(&intranode_comm);

  // Finally, split the specified number of procs_per_node from the super_comm
  // We only want the comm where color == 0.  The second comm is destroyed.
  int color = intranode_comm_rank < procs_per_node ? 0 : 1;
  MPI_Comm_split(super_comm, color, super_comm_rank, sub_comm);
  if (color != 0)
    MPI_Comm_free(sub_comm);
}

} // namespace stream