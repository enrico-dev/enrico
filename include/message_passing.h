#ifndef STREAM_MESSAGE_PASSING_H
#define STREAM_MESSAGE_PASSING_H

#include "mpi.h"

namespace stream {

void get_node_comms(MPI_Comm super_comm, int procs_per_node, MPI_Comm* sub_comm,
                    MPI_Comm* intranode_comm);

} // namespace stream

#endif //STREAM_MESSAGE_PASSING_H
