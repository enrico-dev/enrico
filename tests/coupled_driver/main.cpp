#include "drivers.h"
#include "mpi.h"
#include <algorithm>

struct ProcInfo{
  MPI_Comm comm = MPI_COMM_NULL;
  MPI_Group group = MPI_GROUP_NULL;
  int size = 0;
  int rank = MPI_PROC_NULL;
};

bool isValidRank (int rank) {return rank != MPI_PROC_NULL;}

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  ProcInfo nekGlobal;
  MPI_Comm_dup(MPI_COMM_WORLD, &nekGlobal.comm);
  MPI_Comm_group(nekGlobal.comm, &nekGlobal.group);
  MPI_Comm_rank(nekGlobal.comm, &nekGlobal.rank);
  MPI_Comm_size(nekGlobal.comm, &nekGlobal.size);

  ProcInfo nekNode;
  MPI_Comm_split_type(nekGlobal.comm, MPI_COMM_TYPE_SHARED, nekGlobal.rank, MPI_INFO_NULL, &nekNode.comm);
  MPI_Comm_group(nekNode.comm, &nekNode.group);
  MPI_Comm_rank(nekNode.comm, &nekNode.rank);
  MPI_Comm_size(nekNode.comm, &nekNode.size);

  ProcInfo openmcNode;
  MPI_Group_incl(nekNode.group, 1, {0}, &openmcNode.group);
  MPI_Comm_create(nekNode.comm, openmcNode.group, &openmcNode.comm);
  MPI_Comm_rank(openmcNode.comm, &openmcNode.rank);
  MPI_Comm_size(openmcNode.comm, &openmcNode.size);

  ProcInfo openmcGlobal;
  int transRanksIn[] = {MPI_PROC_NULL}, transRanksOut[] = {MPI_PROC_NULL};
  if (openmcNode.group != MPI_GROUP_NULL) {
    transRanksIn[0] = MPI_PROC_NULL;
    MPI_Group_translate_ranks(openmcNode.group, 1, transRanksIn, nekGlobal.group, transRanksOut);
  }
  auto *translatedRanks = new int[nekGlobal.size];
  MPI_Allgather(transRanksOut, 1, MPI_INT, translatedRanks, 1, MPI_INT, nekGlobal.comm);
  if (openmcNode.group != MPI_GROUP_NULL) {

    // ROR: Can I get this array with an iterator?
    auto numValidRanks = std::count_if(translatedRanks, translatedRanks+nekGlobal.size, isValidRank);
    auto *validRanks = new int[numValidRanks];


  }


  delete []translatedRanks;


  MPI_Finalize();
  return 0;
}
