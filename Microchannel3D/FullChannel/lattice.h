#ifndef LATTICE_H
#define LATTICE_H

#include <vtkXMLStructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include "params_list.h"


template<typename D>
class Lattice
{
    public:
    Lattice(int _zbegin,int _zend,int _NX,int _NY,int _NZ,int _NUMLOCAL,ParamsList params):
        zbegin(_zbegin),zend(_zend),
        NX(_NX),NY(_NY),NZ(_NZ),
        NUMLOCAL(_NUMLOCAL),
        aconst(params("aconst").value<double>()),
        gammaconst(params("gammaconst").value<double>()),
        kconst(params("kconst").value<double>()),
        force_x(params("force_x").value<double>()),
        force_y(params("force_y").value<double>()),
        force_z(params("force_z").value<double>()),
        tau_gas(params("tau_gas").value<double>()),
        tau_liq(params("tau_liq").value<double>()),
        tau_phi(params("tau_phi").value<double>()),
        rho_wall(params("rho_wall").value<double>()),
        phase_wall(params("phase_wall").value<double>()),
        phase_gradient(params("phase_gradient").value<double>())
    {
        //Populations
        f = new double[NUMLOCAL*D::NPOP];
        f2= new double[NUMLOCAL*D::NPOP];
        g = new double[NUMLOCAL*D::NPOP];
        g2= new double[NUMLOCAL*D::NPOP];

        //Macro variables
        phase = new double[NUMLOCAL];
        rho = new double[NUMLOCAL];
        ux  = new double[NUMLOCAL];
        uy  = new double[NUMLOCAL];
        uz  = new double[NUMLOCAL];

        //Helpful materials
        phase_top = new double[NX*NY];
        phase_bottom = new double[NX*NY];
        phase_temp = new double[NX*NY];
        layer_top_f = new double[NX*NY*D::NPOP];
		layer_bottom_f = new double[NX*NY*D::NPOP];
		layer_top_g = new double[NX*NY*D::NPOP];
		layer_bottom_g = new double[NX*NY*D::NPOP];
        layer_temp=new double[NX*NY*D::NPOP];

        //Symmetric populations
        symmetricx[0]=0;
        for(int k=1;k<D::NPOP;k++)
            for(int l=1;l<D::NPOP;l++)
                if ((-D::cx[k]==D::cx[l])&&(D::cy[k]==D::cy[l])&&(D::cz[k]==D::cz[l]))
                {
                    symmetricx[k]=l;
                    break;
                }
        symmetricy[0]=0;
        for(int k=1;k<D::NPOP;k++)
            for(int l=1;l<D::NPOP;l++)
                if ((D::cx[k]==D::cx[l])&&(D::cy[k]==-D::cy[l])&&(D::cz[k]==D::cz[l]))
                {
                    symmetricy[k]=l;
                    break;
                }
        symmetricxy[0]=0;
        for(int k=1;k<D::NPOP;k++)
            for(int l=1;l<D::NPOP;l++)
                if ((D::cx[k]==-D::cx[l])&&(D::cy[k]==-D::cy[l])&&(D::cz[k]==D::cz[l]))
                {
                    symmetricxy[k]=l;
                    break;
                }

        for(int k=0;k<D::NPOP;k++)
        {
            std::cout<<"Symmetricx["<<k<<"]="<<symmetricx[k]<<"\n";
            std::cout<<"Symmetricy["<<k<<"]="<<symmetricy[k]<<"\n";
            std::cout<<"Symmetricxy["<<k<<"]="<<symmetricxy[k]<<"\n";
        }
    }

    ~Lattice()
    {
        //Populations
        delete[] f;
        delete[] f2;
        delete[] g;
        delete[] g2;

        //Macro variables
        delete[] phase;
        delete[] rho;
        delete[] ux;
        delete[] uy;
        delete[] uz;

        //Helpful materials
        delete[] phase_top;
        delete[] phase_bottom;
        delete[] phase_temp;
        delete[] layer_top_f;
        delete[] layer_bottom_f;
        delete[] layer_top_g;
        delete[] layer_bottom_g;
        delete[] layer_temp;
      
    }

    int getZbegin()
    {
        return zbegin;
    }
    int getZend()
    {
        return zend;
    }
    int getNUMLOCAL()
    {
        return NUMLOCAL;
    }

    void putDensity(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double _density)
    {
        for(int iZ=iZbegin;iZ<=iZend;iZ++)
            for(int iY=iYbegin;iY<=iYend;iY++)
                for(int iX=iXbegin;iX<=iXend;iX++)
                    rho[iZ*NX*NY+iY*NX+iX]=_density;
    }

    void putPhase(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double _phase)
    {
        for(int iZ=iZbegin;iZ<=iZend;iZ++)
            for(int iY=iYbegin;iY<=iYend;iY++)
                for(int iX=iXbegin;iX<=iXend;iX++)
                    phase[iZ*NX*NY+iY*NX+iX]=_phase;
    }

    void putVelocity(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double * _velocity)
    {
        for(int iZ=iZbegin;iZ<=iZend;iZ++)
            for(int iY=iYbegin;iY<=iYend;iY++)
                for(int iX=iXbegin;iX<=iXend;iX++)
                {
                    ux[iZ*NX*NY+iY*NX+iX]=_velocity[0];
                    uy[iZ*NX*NY+iY*NX+iX]=_velocity[1];
                    uz[iZ*NX*NY+iY*NX+iX]=_velocity[2];
                }
    }

    void updateMacro();

    //Output files
    void writePhase(std::string name);
    void writeDensity(std::string name);

    void writeTextPhase(std::string name);
    void writeTextDensity(std::string name);


    //Initialization from the initial field
    void init();
    void updateWall();
    void updateSymmetric();
    void updateSymmetricPopulations();

    //Collide and stream
    void collide_stream();

    //Prepare the top and the bottom populations to exchange
    void preparePopulations();
    void finishPropagation();

    //Exchange matrices
    void exchangeMatrices();

    private:
        const int NX;
        const int NY;
        const int NZ;
        const int NUMLOCAL;
        const int zbegin;
        const int zend;
        int symmetricx[D::NPOP];
        int symmetricy[D::NPOP];
        int symmetricxy[D::NPOP];
    public:

        double * f;
        double * f2;
        double * g;
        double * g2;
        double * phase;
        double * rho;
        double * ux;
        double * uy;
        double * uz;
        double * phase_top;
        double * phase_bottom;
        double * phase_temp;
   		double * layer_top_f;
		double * layer_bottom_f;

		double * layer_top_g;
		double * layer_bottom_g;
                double * layer_temp;

		//Parameters of the binary liquid model
		const double aconst;
		const double kconst;
		const double gammaconst;
		const double force_x;
		const double force_y;
		const double force_z;
		const double tau_gas;
		const double tau_liq;
		const double tau_phi;
		const double rho_wall;
		const double phase_wall;
		const double phase_gradient;

};
#endif
