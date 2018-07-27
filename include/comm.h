#ifndef STREAM_COMM_H
#define STREAM_COMM_H

#include "mpi.h"

namespace stream {

class Comm {
public:
  MPI_Comm comm = MPI_COMM_NULL;
  MPI_Group group = MPI_GROUP_NULL;
  int size = 0;
  int rank = MPI_PROC_NULL;

  Comm() {};
  explicit Comm(MPI_Comm comm) : comm(comm) {
    if (comm != MPI_COMM_NULL) {
      MPI_Comm_group(comm, &group);
      MPI_Comm_rank(comm, &rank);
      MPI_Comm_size(comm, &size);
    }
  }
};

} // namespace stream

#endif //STREAM_COMM_H
