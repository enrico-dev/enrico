#include "drivers.h"
#include "mpi.h"
#include <algorithm>
#include <vector>
#include "unistd.h"
#include <climits>
#include <string>
#include <iostream>

struct ProcInfo{
  MPI_Comm comm = MPI_COMM_NULL;
  MPI_Group group = MPI_GROUP_NULL;
  int size = 0;
  int rank = MPI_PROC_NULL;
};

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  //===========================================================================
  // ProcInfo world: A dup of MPI_COMM_WORLD
  //===========================================================================

  ProcInfo world;
  MPI_Comm_dup(MPI_COMM_WORLD, &world.comm);
  MPI_Comm_group(world.comm, &world.group);
  MPI_Comm_rank(world.comm, &world.rank);
  MPI_Comm_size(world.comm, &world.size);

  //===========================================================================
  // ProcInfo allProcsInMyShmem: Use MPI_comm_split_type to get a new comm in
  // every shmem region.  All procs in a shmem region will belong to the same
  // comm.
  //
  // The shmem region is intended to be a node; the behavior has been
  // untested for more unusual architectures.
  //===========================================================================

  ProcInfo allProcsInMyShmem;
  MPI_Comm_split_type(world.comm, MPI_COMM_TYPE_SHARED, world.rank, MPI_INFO_NULL,
                      &allProcsInMyShmem.comm);
  MPI_Comm_group(allProcsInMyShmem.comm, &allProcsInMyShmem.group);
  MPI_Comm_rank(allProcsInMyShmem.comm, &allProcsInMyShmem.rank);
  MPI_Comm_size(allProcsInMyShmem.comm, &allProcsInMyShmem.size);

  //===========================================================================
  // ProcInfo oneProcInMyShmem: In each group from allProcsInMyShmem, use
  // MPI_group_incl to select one proc.  There will be a different one of these
  // groups in every shmem domain.  Most procs will not be in these groups.
  //===========================================================================

  ProcInfo oneProcInMyShmem;
  MPI_Group_incl(allProcsInMyShmem.group, 1, {0}, &oneProcInMyShmem.group);
  MPI_Comm_create(allProcsInMyShmem.comm, oneProcInMyShmem.group, &oneProcInMyShmem.comm);
  MPI_Comm_rank(oneProcInMyShmem.comm, &oneProcInMyShmem.rank);
  MPI_Comm_size(oneProcInMyShmem.comm, &oneProcInMyShmem.size);

  //===========================================================================
  // ProcInfo oneProcInAllShmems: Collect all the procs in all groups from
  // oneProcInMyShmem.  Create a single new group with one proc from every
  // shmem.
  //===========================================================================

  // First, translate each rank in a oneProcInMyShmem group to a rank in world.
  // The ranks that aren't in a oneProcInMyShmem group will be translated to
  // MPI_PROC_NULL.
  int myTransRankOut = MPI_PROC_NULL;
  int myTransRankIn = oneProcInMyShmem.group != MPI_GROUP_NULL ? oneProcInMyShmem.rank : MPI_PROC_NULL;
  if (oneProcInMyShmem.group != MPI_GROUP_NULL)
    MPI_Group_translate_ranks(oneProcInMyShmem.group, 1, &myTransRankIn, world.group, &myTransRankOut);

  // Next, gather all the translated (world) ranks from *all* procs.
  std::vector<int> allTransRanks(world.size, MPI_PROC_NULL);
  MPI_Allgather(&myTransRankOut, 1, MPI_INT, allTransRanks.data(), 1, MPI_INT, world.comm);

  // Next, remove the invalid ranks.  This will leave us with only the procs
  // that are in a oneProcInMyShmem group
  auto validRanks = allTransRanks;
  auto validEnd = std::remove_if(validRanks.begin(), validRanks.end(), [](int v){return v == MPI_PROC_NULL;});
  validRanks.erase(validRanks.begin(), validEnd);

  // Finally, make a new group and comm from the valid ranks.
  ProcInfo oneProcInAllShmems;
  MPI_Group_incl(world.comm, static_cast<int>(validRanks.size()), validRanks.data(), &oneProcInAllShmems.group);
  MPI_Comm_create(world.comm, oneProcInAllShmems.group, &oneProcInAllShmems.comm);
  if (oneProcInAllShmems.comm == MPI_COMM_NULL) {
    MPI_Comm_rank(oneProcInAllShmems.comm, &oneProcInAllShmems.rank);
    MPI_Comm_size(oneProcInAllShmems.comm, &oneProcInAllShmems.size);
  }

  //===========================================================================
  // Debug by printing node memberships
  //===========================================================================

  // Can't do:
  //    std::string myHostName(HOST_NAME_MAX, ' ');
  //    char *cMyHostName = myHostName.c_str()
  // since c_str() returns a const pointer and can't be passed to gethostname

  // Get hostname
  auto cMyHostName = new char[HOST_NAME_MAX];
  gethostname(cMyHostName, HOST_NAME_MAX);
  std::string myHostName(cMyHostName);

  // Pring hostnames
  if (world.rank == 0)
    std::cout << "Ranks/hosts in world:" << std::endl;
  MPI_Barrier(world.comm);
  for (auto i; i < world.size; i++) {
    if (world.rank == i)
      std::cout << world.rank << myHostName << std::endl;
  }
  MPI_Barrier(world.comm);

  if (world.rank == 0)
    std::cout << std::endl << "Ranks/hosts in oneProcInAllShmems:" << std::endl;
  MPI_Barrier(world.comm);
  for (auto i = 0; i < oneProcInAllShmems.size; i++) {
    if (oneProcInAllShmems.rank != MPI_PROC_NULL and oneProcInAllShmems.rank == i)
      std::cout << oneProcInAllShmems.rank << myHostName << std::endl;
  }
  MPI_Barrier(world.comm);

  MPI_Finalize();
  return 0;
}
