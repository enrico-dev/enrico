//! \file comm.h
//! Info and function wrappers for a specified MPI communictor.
#ifndef ENRICO_COMM_H
#define ENRICO_COMM_H

#include "enrico/mpi_types.h"
#include "xtensor/xtensor.hpp"

#include <mpi.h>

#include <iostream>
#include <string>
#include <vector>

namespace enrico {

//! Info and function wrappers for a specified MPI communictor.
class Comm {
public:
  //! Default constructor
  Comm() = default;

  //! Retrieves info about a given MPI communicator.
  //! \param comm An exisiting or null MPI communicator.
  explicit Comm(MPI_Comm comm)
    : comm(comm)
  {
    if (comm != MPI_COMM_NULL) {
      MPI_Comm_group(comm, &group);
      MPI_Comm_rank(comm, &rank);
      MPI_Comm_size(comm, &size);
    }
  }

  //! Frees the underlying MPI Communicator and nullifies/zeroes the group, rank, and size
  //! \return Error value
  int free()
  {
    auto ierr = MPI_Comm_free(&comm);
    if (ierr == MPI_SUCCESS) {
      group = MPI_GROUP_NULL;
      rank = MPI_PROC_NULL;
      size = 0;
    }
    return ierr;
  };

  //! Queries whether the communicator is active
  //! \return True if the communicator is not MPI_COMM_NULL
  bool active() const { return comm != MPI_COMM_NULL; }

  //! Queries whether the calling rank is the root of this comm
  //! \return True if the calling rank is the root of this comm
  bool is_root() const { return rank == 0; }

  //! Block until all processes have reached this call
  //!
  //! \return Error value
  int Barrier() const { return MPI_Barrier(comm); }

  //! Broadcasts a message from the process with rank "root" to all other processes in
  //! this comm.
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

  //! Broadcast a scalar value across ranks
  //! \param value Value to broadcast (significant at rank 0)
  template<typename T>
  std::enable_if_t<std::is_scalar<std::decay_t<T>>::value> broadcast(T& value) const;

  //! Broadcast a vector across ranks, possibly resizing it
  //! \param values Values to broadcast (significant at rank 0)
  template<typename T>
  void broadcast(std::vector<T>& values) const;

  //! Broadcast an xtensor across ranks, possibly resizing it to match root's shape
  //! \param values Values to broadcast (significant at rank 0)
  template<typename T, size_t N>
  void broadcast(xt::xtensor<T, N>& values) const;

  //! Send a scalar from one rank to another
  //! \param value Value to broadcast (significant at source rank)
  //! \param dest Destination rank
  //! \param source Source rank
  template<typename T>
  std::enable_if_t<std::is_scalar<std::decay_t<T>>::value>
  sendrecv_replace(T& value, int dest, int source) const;

  //! Send a vector from one rank to another, possibly resizing it at destination
  //! \param values Values to broadcast (significant at source rank)
  //! \param dest Destination rank
  //! \param source Source rank
  template<typename T>
  void sendrecv_replace(std::vector<T>& values, int dest, int source) const;

  //! Send an xtensor from one rank to another, possibly resizing it at destination
  //! \param values Values to broadcast (significant at source rank)
  //! \param dest Destination rank
  //! \param source Source rank
  template<typename T, size_t N>
  void sendrecv_replace(xt::xtensor<T, N>& values, int dest, int source) const;

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
  int Gather(const void* sendbuf,
             int sendcount,
             MPI_Datatype sendtype,
             void* recvbuf,
             int recvcount,
             MPI_Datatype recvtype,
             int root = 0) const
  {
    return MPI_Gather(
      sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  }

  //! Gathers into specified locations from all processes onto a given root
  //!
  //! Currently, a wrapper for MPI_Gatherv.
  //!
  //! \param[in] sendbuf Starting address of send buffer
  //! \param[in] sendcount Number of elements in send buffer
  //! \param[in] sendtype Data type of send buffer elements
  //! \param[out] recvbuf Address of receive buffer
  //! \param[in] recvcounts Integer array (of length group size) containing the number of
  //! elements
  //!                       that are received from each process
  //! \param[in] displs Integer array (of length group size). Entry i specifies the
  //! displacement
  //!                   relative to recvbuf at which to place the incoming data from
  //!                   process i.
  //! \param[in] recvtype Data type of recv buffer elements
  //! \param[in] root Rank of receiving process
  //! \return Error value
  int Gatherv(const void* sendbuf,
              int sendcount,
              MPI_Datatype sendtype,
              void* recvbuf,
              const int recvcounts[],
              const int displs[],
              MPI_Datatype recvtype,
              int root = 0) const
  {
    return MPI_Gatherv(
      sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm);
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
  int Allgather(const void* sendbuf,
                int sendcount,
                MPI_Datatype sendtype,
                void* recvbuf,
                int recvcount,
                MPI_Datatype recvtype) const
  {
    return MPI_Allgather(
      sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  }

  //! Displays a message from rank 0
  //! \param A message to display
  void message(const std::string& msg)
  {
    if (rank == 0)
      std::cout << "[ENRICO]: " << msg << std::endl;
  }

  // Data members
  MPI_Comm comm =
    MPI_COMM_NULL; //!< The MPI communicator described by this instance of Comm.
  MPI_Group group = MPI_GROUP_NULL; //!< The group associated with Comm::comm.
  int size = 0;                     //!< The size of Comm::comm.
  int rank = MPI_PROC_NULL;         //!< The calling process's rank in Comm:comm
};

template<typename T>
std::enable_if_t<std::is_scalar<std::decay_t<T>>::value>
Comm::sendrecv_replace(T& value, int dest, int source) const
{
  MPI_Sendrecv_replace(
    &value, 1, get_mpi_type<T>(), dest, 0, source, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
};

template<typename T>
void Comm::sendrecv_replace(std::vector<T>& values, int dest, int source) const
{
  // Send the size of the vector from the source
  auto n = values.size();
  sendrecv_replace(n, dest, source);
  // Resize the vector on dest
  if (rank == dest && values.size() != n) {
    values.resize(n);
  }
  // Send the vector
  MPI_Sendrecv_replace(values.data(),
                       n,
                       get_mpi_type<T>(),
                       dest,
                       0,
                       source,
                       MPI_ANY_TAG,
                       comm,
                       MPI_STATUS_IGNORE);
}

template<typename T, size_t N>
void Comm::sendrecv_replace(xt::xtensor<T, N>& values, int dest, int source) const
{
  // Make sure the shapes match
  const auto& s = values.shape();
  std::vector<size_t> my_shape(s.begin(), s.end());
  std::vector<size_t> root_shape(my_shape);
  sendrecv_replace(root_shape, dest, source);
  if (rank == dest && my_shape != root_shape) {
    values.resize(root_shape);
  }
  // Next, broadcast size
  auto n = values.size();
  sendrecv_replace(n, dest, source);
  // Finally, broadcast data
  MPI_Sendrecv_replace(values.data(),
                       n,
                       get_mpi_type<T>(),
                       dest,
                       0,
                       source,
                       MPI_ANY_TAG,
                       comm,
                       MPI_STATUS_IGNORE);
}

template<typename T>
std::enable_if_t<std::is_scalar<std::decay_t<T>>::value> Comm::broadcast(T& value) const
{
  this->Bcast(&value, 1, get_mpi_type<T>());
}

template<typename T>
void Comm::broadcast(std::vector<T>& values) const
{
  if (this->active()) {
    // First broadcast the size of the vector
    int n = values.size();
    this->broadcast(n);

    // Resize vector (for rank != 0) and broacast data
    if (values.size() != n)
      values.resize(n);
    this->Bcast(values.data(), n, get_mpi_type<T>());
  }
}

template<typename T, size_t N>
void Comm::broadcast(xt::xtensor<T, N>& values) const
{
  if (this->active()) {
    // First, make sure shape of `values` matches root's
    const auto& s = values.shape();
    std::vector<size_t> my_shape(s.begin(), s.end());
    std::vector<size_t> root_shape(my_shape);

    this->broadcast(root_shape);
    if (my_shape != root_shape) {
      values.resize(root_shape);
    }

    // Next, broadcast size
    auto n = values.size();
    this->broadcast(n);

    // Finally, broadcast data
    this->Bcast(values.data(), n, get_mpi_type<T>());
  }
}

} // namespace enrico

#endif // ENRICO_COMM_H
