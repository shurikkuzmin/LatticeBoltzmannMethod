#include "mpi_singleton.h"

template<> MPIWrapper * MPISingleton::instance = 0;
