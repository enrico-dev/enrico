#include "message_passing.h"
#include "mpi.h"

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
