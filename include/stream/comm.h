//! \file comm.h
//! Info and function wrappers for a specified MPI communictor.
#ifndef STREAM_COMM_H
#define STREAM_COMM_H

#include "mpi.h"
#include <string>
#include <iostream>

namespace stream {

//! Info and function wrappers for a specified MPI communictor.
class Comm {
public:
  MPI_Comm comm = MPI_COMM_NULL; //!< The MPI communicator described by this instance of Comm.
  MPI_Group group = MPI_GROUP_NULL; //!< The group associated with Comm::comm.
  int size = 0; //!< The size of Comm::comm.
  int rank = MPI_PROC_NULL; //!< The calling process's rank in Comm:comm

  //! Default constructor
  Comm() {};

  //! Retrieves info about a given MPI communicator.
  //! \param comm An exisiting or null MPI communicator.
  explicit Comm(MPI_Comm comm) : comm(comm)
  {
    if (comm != MPI_COMM_NULL) {
      MPI_Comm_group(comm, &group);
      MPI_Comm_rank(comm, &rank);
      MPI_Comm_size(comm, &size);
    }
  }

  //! Block until all processes have reached this call
  //!
  //! \return Error value
  int Barrier() const
  {
    return MPI_Barrier(comm);
  }

  //! Broadcasts a message from the process with rank "root" to all other processes in this comm.
  //!
  //! Currently, a wrapper for MPI_Bcast.
  //!
  //! \param[in,out] buffer Starting address of buffer
  //! \param[in] count Number of entries in buffer
  //! \param[in] datatype Data type of buffer
  //! \param[in] root Rank of broadcast root
  //! \return Error value
  int Bcast(void* buffer, int count, MPI_Datatype datatype, int root = 0) const
  {
    return MPI_Bcast(buffer, count, datatype, root, comm);
  }

  //! Gathers together values from the processes in this comm onto a given root.
  //!
  //! Currently, a wrapper for MPI_Gather.
  //!
  //! \param[in] sendbuf Starting address of send buffer
  //! \param[in] sendcount Number of elements in send buffer
  //! \param[in] sendtype Data type of send buffer elements
  //! \param[out] recvbuf Address of receive buffer
  //! \param[in] recvcount Number of elements for any single receive
  //! \param[in] recvtype Data type of recv buffer elements
  //! \param[in] root Rank of receiving process
  //! \return Error value
  int Gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
             void* recvbuf, int recvcount, MPI_Datatype recvtype,
             int root = 0) const
  {
    return MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                      recvtype, root, comm);
  }

  //! Gathers into specified locations from all processes onto a given root
  //!
  //! Currently, a wrapper for MPI_Gatherv.
  //!
  //! \param[in] sendbuf Starting address of send buffer
  //! \param[in] sendcount Number of elements in send buffer
  //! \param[in] sendtype Data type of send buffer elements
  //! \param[out] recvbuf Address of receive buffer
  //! \param[in] recvcounts Integer array (of length group size) containing the number of elements
  //!                       that are received from each process
  //! \param[in] displs Integer array (of length group size). Entry i specifies the displacement
  //!                   relative to recvbuf at which to place the incoming data from process i.
  //! \param[in] recvtype Data type of recv buffer elements
  //! \param[in] root Rank of receiving process
  //! \return Error value
  int Gatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
              void* recvbuf, const int recvcounts[], const int displs[],
              MPI_Datatype recvtype, int root = 0) const
  {
    return MPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts,
                       displs, recvtype, root, comm);
  }

  //! Gathers data from all tasks and distribute the combined data to all tasks.
  //!
  //! Currently, a wrapper for MPI_Allgather
  //!
  //! \param[in] sendbuf Starting address of send buffer
  //! \param[in] sendcount Number of elements in send buffer
  //! \param[in] sendtype Data type of send buffer elements
  //! \param[out] recvbuf  Starting address of receive buffer
  //! \param[in] recvcount Number of elements received from any process
  //! \param[in] recvtype Data type of receive buffer elements
  //! \return
  int Allgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
                void* recvbuf, int recvcount, MPI_Datatype recvtype) const
  {
    return MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  }

  //! Displays a message from rank 0
  //! \param A message to display
  void message(const std::string& msg)
  {
    if (rank == 0) std::cout << "[STREAM]: " << msg << std::endl;
  }
};

} // namespace stream

#endif //STREAM_COMM_H
