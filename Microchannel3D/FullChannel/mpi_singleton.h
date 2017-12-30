#ifndef __MIXER3D_MPI_SINGLETON_H__
#define __MIXER3D_MPI_SINGLETON_H__

#include "singleton.h"
#include "mpi_wrapper.h"

typedef Singleton<MPIWrapper> MPISingleton;

template<> inline MPIWrapper * MPISingleton::create_instance()
{
	return new MPIWrapper(0, 0);
}

inline MPIWrapper & mpi_wrapper()
{
	return MPISingleton::get_instance();
}

#endif
