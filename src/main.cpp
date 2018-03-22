#include "mpi.h"
#include "drivers.h"

int main(int argc, char* argv[])
{
  OpenmcDriver *o = new OpenmcDriver(MPI_COMM_WORLD);
  return 0;
}

