#include "message_passing.h"
#include "mpi.h"

/**
 * @brief Splits a given MPI comunicator into a new comm with a specified number of procs in every node.
 *
 * Tne new comm (*subComm*) will span **all** nodes in the original comm (*superComm*).  This will work correctly even
 * if *subComm* spans only a single node.
 *
 * Depending on the calling proc, *subComm* will be one of two values:
 *   - If the calling proc is within the new desired comm (the comm with the given number of *procsPerNode*), then
 *     *subComm* will be a valid comm.
 *   - If the calling proc is not within the new desired comm, then *subComm* will be `MPI_COMM_NULL`.
 * The caller must then make use of *subComm* with proper value checks.  For example, subsequent MPI operations on
 * *subComm* may need to check if `subComm != MPI_COMM_NULL`.
 *
 * @param superComm An existing communicator that will be split
 * @param procsPerNode The number of MPI procs per node in the new communicator, *subComm*
 * @param subComm A new communicator with either the given number of procs per node *or* `MPI_COMM_NULL`, depending on
 *                whether the calling proc is in the desired comm.
 */
void getInternodeSubComm(MPI_Comm superComm, int procsPerNode, MPI_Comm *subComm) {

  // superCommRank is used as the "key" to retain ordering in the comm splits.
  // This can allow the subComm to retain some intent from the superComm's proc layout
  int superCommRank = MPI_PROC_NULL;
  MPI_Comm_rank(superComm, &superCommRank);

  // intranodeComm is an intermediate object.  It is only used to get an intranodeCommRank,
  // which is used as the "color" in the final comm split.
  MPI_Comm intranodeComm = MPI_COMM_NULL;
  MPI_Comm_split_type(superComm, MPI_COMM_TYPE_SHARED, superCommRank, MPI_INFO_NULL, &intranodeComm);
  int intranodeCommRank = MPI_PROC_NULL;
  MPI_Comm_rank(intranodeComm, &intranodeCommRank);
  MPI_Comm_free(&intranodeComm);

  // Finally, split the specified number of procsPerNode from the superComm
  // We only want the comm where color == 0.  The second comm is destroyed.
  int color = intranodeCommRank < procsPerNode ? 0 : 1;
  MPI_Comm_split(superComm, color, superCommRank, subComm);
  if (color != 0)
    MPI_Comm_free(subComm);
}
