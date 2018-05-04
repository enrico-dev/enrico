#include <climits>
#include <iostream>
#include <string>
#include "mpi.h"
#include "unistd.h"

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
  world.comm = MPI_COMM_WORLD;
  MPI_Comm_group(world.comm, &world.group);
  MPI_Comm_rank(world.comm, &world.rank);
  MPI_Comm_size(world.comm, &world.size);

  //===========================================================================
  // ProcInfo allProcsInMyShmem: Use MPI_comm_split_type to get a new comm in
  // every shmem region.  All procs in a shmem region will belong to the same
  // comm.
  //
  // The shmem region is intended to be a node; the behavior has not been
  // tested for more unusual architectures.
  //===========================================================================

  ProcInfo allProcsInMyShmem;
  MPI_Comm_split_type(world.comm, MPI_COMM_TYPE_SHARED, world.rank, MPI_INFO_NULL,
                      &allProcsInMyShmem.comm);
  MPI_Comm_group(allProcsInMyShmem.comm, &allProcsInMyShmem.group);
  MPI_Comm_rank(allProcsInMyShmem.comm, &allProcsInMyShmem.rank);
  MPI_Comm_size(allProcsInMyShmem.comm, &allProcsInMyShmem.size);

  //===========================================================================
  // ProcInfo oneProcInAllShmems: In each shmem domain, choose one proc
  // (arbitrarily, we choose one such that allProcsInMyShmem.rank == 0).
  // Split all these procs (from all shmem) off of world into a new comm.
  //===========================================================================

  ProcInfo oneProcInAllShmems;
  int myColor = allProcsInMyShmem.rank == 0 ? 0 : 1;
  MPI_Comm_split(world.comm, myColor, world.rank, &oneProcInAllShmems.comm);
  if (allProcsInMyShmem.rank != 0) {
    MPI_Comm_free(&oneProcInAllShmems.comm);
  }
  else {
    MPI_Comm_group(oneProcInAllShmems.comm, &allProcsInMyShmem.group);
    MPI_Comm_rank(oneProcInAllShmems.comm, &oneProcInAllShmems.rank);
    MPI_Comm_size(oneProcInAllShmems.comm, &oneProcInAllShmems.size);
  }

  //===========================================================================
  // Debug: Print which comms/ranks are on which host
  //===========================================================================

  auto cMyHostName = new char[HOST_NAME_MAX];
  gethostname(cMyHostName, HOST_NAME_MAX);
  std::string myHostName(cMyHostName);

  for (auto i = 0; i < world.size; i++) {
    if (world.rank == i)
      std::cout << "world (rank,host)\t" << world.rank << "\t" << myHostName << std::endl;
    MPI_Barrier(world.comm);
  }

  for (auto i = 0; i < world.size; i++) {
    if (oneProcInAllShmems.rank != MPI_PROC_NULL and world.rank == i)
      std::cout << "oneProcInAllShmems (rank,host)\t" << oneProcInAllShmems.rank << "\t" << myHostName << std::endl;
    MPI_Barrier(world.comm);
  }

  MPI_Finalize();
  return 0;
}
