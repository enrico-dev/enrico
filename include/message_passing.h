#ifndef STREAM_MESSAGE_PASSING_H
#define STREAM_MESSAGE_PASSING_H

#include "mpi.h"

void getInternodeSubComm(MPI_Comm superComm, int procsPerNode, MPI_Comm *subComm);

#endif //STREAM_MESSAGE_PASSING_H
