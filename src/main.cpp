#include "mpi.h"
#include "drivers.h"

int main(int argc, char* argv[])
{
  OpenmcDriver *o = new OpenmcDriver(argc, argv, MPI_COMM_WORLD);
  return 0;
}
