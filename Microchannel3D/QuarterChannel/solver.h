#ifndef SOLVER_H
#define SOLVER_H

#include "lattice.h"
#include "mpi_singleton.h"

#include "destroyer.h"
#include "params_list.h"

template<typename D>
class Solver
{
    public:
    explicit Solver(ParamsList _params_list);
    ~Solver()
    {
        delete lattice;
        if (rank==0)
        {
            delete[] density;
            delete[] phase;
            delete[] velocity;
        }
    }

    bool checkNAN();

    void load_file(std::string name);

    void putDensity(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double dense);
    void putDensity(double dense);
    void putPhase(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double dense);
    void putPhase(double dense);
    void putVelocity(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double * velocity);
    void putVelocity(double * velocity);

    void getWholeDensity();
    void getWholePhase();
    void getWholeVelocity();

    void writeWholeDensity(std::string name);
    void writeWholePhase(std::string name);
    void writeWholeVelocity(std::string name);
    void writeWholeDensityPhaseVelocity(std::string name);

    void writeTextXZVelocity(std::string name);
    void writeTextXYVelocity(std::string name);
    void writeTextYZVelocity(std::string name);


    void writeLocalDensity(std::string name);
    void writeLocalPhase(std::string name);
    void writeLocalTextDensity(std::string name);
    void writeLocalTextPhase(std::string name);

    void writeTextWholeVelocity(std::string name);
    void writeTextWholePhase(std::string name);
    void writeTextWholeDensity(std::string name);
    void writeVTKWholePhase(std::string name);

    void writeTextXZPhase(std::string name);
    void writeTextXYPhase(std::string name);
    void writeTextYZPhase(std::string name);

    void writeTextXZDensity(std::string name);
    void writeTextXYDensity(std::string name);
    void writeTextYZDensity(std::string name);

    void init();
	void updatePhase();
    void propagate();
    void collide_stream();
    void exchangeMatrices();


    Lattice<D> * getLattice()
    {
        return lattice;
    }


	private:

    const int size;
    const int rank;
    const int NX;
    const int NY;
    const int NZ;


    Lattice<D> * lattice;
    double * density;
    double * phase;
    double * velocity;
};

#endif
