#include "drivers.h"
#include "mpi.h"

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  OpenmcDriver driver(MPI_COMM_WORLD);
  driver.initStep();
  driver.solveStep();
  driver.finalizeStep();
  MPI_Finalize();
}
