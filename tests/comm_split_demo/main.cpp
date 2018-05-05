//===========================================================================
// Demonstrates splitting comm world into another comm with a specified
// number of procs per node.  The latter comm is named nProcsInAllShmems
//
// Usage:
//   mpirun -np <nWorldProcs> ./a.out <nProcsPerNode>
//
// Output:
//   Prints ranks and nodes for the two comms
//===========================================================================

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
  // Get procsPerNode from command line
  //===========================================================================

  int procsPerNode = 1;
  if (argc > 1)
    procsPerNode = std::stoi(argv[1]);

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
  // This is an intermediate step on the way to creating nProcsInAllShmems.
  // The latter is the ultimate objective.  
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
  // ProcInfo nProcsInAllShmems: In each shmem domain, choose n procs
  // (arbitrarily, we choose one such that allProcsInMyShmem.rank < nProcsPerNode).
  // Split all these procs (from all shmem) off of world into a new comm.
  //===========================================================================

  ProcInfo nProcsInAllShmems;
  int myColor = allProcsInMyShmem.rank < procsPerNode ? 0 : 1;
  MPI_Comm_split(world.comm, myColor, world.rank, &nProcsInAllShmems.comm);
  if (myColor == 0) {
    MPI_Comm_group(nProcsInAllShmems.comm, &allProcsInMyShmem.group);
    MPI_Comm_rank(nProcsInAllShmems.comm, &nProcsInAllShmems.rank);
    MPI_Comm_size(nProcsInAllShmems.comm, &nProcsInAllShmems.size);
  }
  else
    MPI_Comm_free(&nProcsInAllShmems.comm);

  //===========================================================================
  // Debug: Print which comms/ranks are on which host
  //===========================================================================

  auto cMyHostName = new char[HOST_NAME_MAX];
  gethostname(cMyHostName, HOST_NAME_MAX);
  std::string myHostName(cMyHostName);

  for (int i = 0; i < world.size; i++) {
    if (world.rank == i)
      std::cout << "world (rank,host)\t" << world.rank << "\t" << myHostName << std::endl;
    MPI_Barrier(world.comm);
  }

  for (int i = 0; i < world.size; i++) {
    if (nProcsInAllShmems.rank != MPI_PROC_NULL and world.rank == i)
      std::cout << "nProcsInAllShmems (rank,host)\t" << world.rank << "\t" << myHostName << std::endl;
    MPI_Barrier(world.comm);
  }

  MPI_Finalize();
  return 0;
}
