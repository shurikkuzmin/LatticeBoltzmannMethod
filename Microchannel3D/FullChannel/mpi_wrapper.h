#ifndef __MIXER3D_MPI_WRAPPER_H__
#define __MIXER3D_MPI_WRAPPER_H__

#include <string>

class MPIWrapper
{
public:

	MPIWrapper(int argc, char ** argv);
	~MPIWrapper();

	int get_size()const {return size;}
	int get_rank()const {return rank;}

	const std::string & get_name()const {return name;}

private:

	int size;
	int rank;
	std::string name;
};

#endif
