#include "mpi.h"
#include "drivers.h"

int main(int argc, char* argv[])
{
  {
    stream::OpenmcDriver o {argc, argv, MPI_COMM_WORLD};
  }
  return 0;
}
