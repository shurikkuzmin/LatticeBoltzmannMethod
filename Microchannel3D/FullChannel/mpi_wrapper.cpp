#include <mpi.h>
#include "mpi_wrapper.h"


MPIWrapper::MPIWrapper(int argc, char ** argv)
{
	MPI_Init(&argc, &argv);

	// get this processor's name
	char buf[MPI_MAX_PROCESSOR_NAME];
	int buf_len = 0;
	MPI_Get_processor_name(buf, &buf_len);
	name = std::string(buf);

	// get this processor's number
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// get the number of processors
	MPI_Comm_size(MPI_COMM_WORLD, &size);
}

MPIWrapper::~MPIWrapper()
{
	MPI_Finalize();
}
