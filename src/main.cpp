#include "mpi.h"
#include "drivers.h"

int main(int argc, char* argv[])
{
  OmcDriver *o = new OmcDriver(MPI_COMM_WORLD);
  return 0;
}

