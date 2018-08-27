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
  explicit Comm(MPI_Comm comm) : comm(comm)
  {
    if (comm != MPI_COMM_NULL) {
      MPI_Comm_group(comm, &group);
      MPI_Comm_rank(comm, &rank);
      MPI_Comm_size(comm, &size);
    }
  }

  int Bcast(void* buffer, int count, MPI_Datatype datatype, int root = 0)
  {
    return MPI_Bcast(buffer, count, datatype, root, comm);
  }

  int Gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
             void* recvbuf, int recvcount, MPI_Datatype recvtype,
             int root = 0)
  {
    return MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                      recvtype, root, comm);
  }

  int Gatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
              void* recvbuf, const int recvcounts[], const int displs[],
              MPI_Datatype recvtype, int root = 0)
  {
    return MPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts,
                       displs, recvtype, root, comm);
  }
};

} // namespace stream

#endif //STREAM_COMM_H
