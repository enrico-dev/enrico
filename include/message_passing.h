#ifndef STREAM_MESSAGE_PASSING_H
#define STREAM_MESSAGE_PASSING_H

#include "mpi.h"

namespace stream {

void getInternodeSubComm(MPI_Comm superComm, int procsPerNode, MPI_Comm *subComm);

} // namespace stream

#endif //STREAM_MESSAGE_PASSING_H
