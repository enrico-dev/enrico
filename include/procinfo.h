#ifndef STREAM_PROCINFO_H
#define STREAM_PROCINFO_H

#include "mpi.h"

class ProcInfo{
public:
  MPI_Comm comm = MPI_COMM_NULL;
  MPI_Group group = MPI_GROUP_NULL;
  int size = 0;
  int rank = MPI_PROC_NULL;

  ProcInfo() {};
  explicit ProcInfo(MPI_Comm comm) : comm(comm) {
    MPI_Comm_group(comm, &group);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
  }
};

#endif //STREAM_PROCINFO_H
