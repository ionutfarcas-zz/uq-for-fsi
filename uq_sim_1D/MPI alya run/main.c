#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>

int main(int argc, char** argv)
{
	int rank = 0;
	int nprocs = 0;

	int run_ok = 0;
	int i = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	run_ok = system("./runs");
	assert(run_ok >= 0);
	
	MPI_Finalize();
	return 0;
}
