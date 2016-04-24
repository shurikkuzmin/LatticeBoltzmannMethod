#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <limits>

//Domain size
const int NY=128;
const int NX=128;

int radius=20;
//Time steps
int N=20000;
int NOUTPUT=100;
int NSIGNAL=100;
double force_y;
double force_x;

//Fields and populations
double f[NX][NY][9], f2[NX][NY][9], g[NX][NY][9], g2[NX][NY][9];
double rho[NX][NY],ux[NX][NY],uy[NX][NY],phase[NX][NY];


//Binary-liquid parameters
double aconst=0.04;
double kconst=0.04;
double gammaconst=1.0;
double wall_gradient=-0.15;

//double tau_gas=0.7;
//double tau_liq=2.5;
double tau_gas=1.0;
double tau_liq=1.0;

//BGK relaxation parameter
double omega_rho=1.0;
double omega_phi=1.0;

//Underlying lattice parameters
double weights[]={4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};
int cx[]={0,1,0,-1,0,1,-1,-1,1};
int cy[]={0,0,1,0,-1,1,1,-1,-1};
int compliment[]={0,3,4,1,2,7,8,5,6};
float wxx[] = {0.0, 1.0/3.0, -1.0/6.0, 1.0/3.0, -1.0/6.0, -1.0/24.0, -1.0/24.0, -1.0/24.0, -1.0/24.0};
float wyy[] = {0.0, -1.0/6.0, 1.0/3.0, -1.0/6.0, 1.0/3.0, -1.0/24.0, -1.0/24.0, -1.0/24.0, -1.0/24.0};
float wxy[] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0/4.0, -1.0/4.0, 1.0/4.0, -1.0/4.0};

float gradstencilx[9]={0.0,4.0/12.0,0.0,-4.0/12.0,0.0,
			          1.0/12.0,-1.0/12.0,-1.0/12.0,1.0/12.0};

float gradstencily[9]={0.0,0.0,4.0/12.0,0.0,-4.0/12.0,
			          1.0/12.0,1.0/12.0,-1.0/12.0,-1.0/12.0};

float laplacestencil[9]={-20.0/6.0,4.0/6.0,4.0/6.0,4.0/6.0,4.0/6.0,
					   1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0};

void writedensity(std::string const & fName)
{
	std::string fullName = fName+ ".dat";
	std::ofstream fout(fullName.c_str());
	fout.precision(10);

	for (int iY=NY-1; iY>=0; --iY)
	{
		for (int iX=0; iX<NX; ++iX)
			fout<<rho[iX][iY]<<" ";
		fout<<"\n";
	}

}

void writephase(std::string const & fName)
{
	std::string fullName = fName+ ".dat";
	std::ofstream fout(fullName.c_str());
	fout.precision(10);

	for (int iY=NY-1; iY>=0; --iY)
	{
		for (int iX=0; iX<NX; ++iX)
			fout<<phase[iX][iY]<<" ";
		fout<<"\n";
	}

}


void writevelocityx(std::string const & fName)
{
	std::string fullName = fName+ ".dat";
	std::ofstream fout(fullName.c_str());
	fout.precision(10);

	for (int iY=NY-1; iY>=0; --iY)
	{
		for (int iX=0; iX<NX; ++iX)
			fout<<ux[iX][iY]<<" ";
		fout<<"\n";
	}

}

void writevelocityy(std::string const & fName)
{
	std::string fullName = fName+ ".dat";
	std::ofstream fout(fullName.c_str());
	fout.precision(10);

	for (int iY=NY-1; iY>=0; --iY)
	{
		for (int iX=0; iX<NX; ++iX)
			fout<<uy[iX][iY]<<" ";
		fout<<"\n";
	}

}

void init()
{

	//Phase initialization prior any equilibrium functions calculations
    for(int iX=0;iX<NX;iX++)
		for(int iY=0; iY<NY; iY++)
			//if ( (iX-(NX-1)/2)*(iX-(NX-1)/2)+(iY-(NY-1)/2)*(iY-(NY-1)/2)<=radius*radius )
            if ((abs(iX-NX/2)<=20)&&(iY<40)&&(iY>=1))
                phase[iX][iY]=1.0;
            else
                phase[iX][iY]=-1.0;
     
    //Walls initialization
    
    for(int iX=0;iX<NX;iX++)
    {
        phase[iX][0]=phase[iX][1]-wall_gradient;
        phase[iX][NY-1]=phase[iX][NY-2]-wall_gradient;
  	}
	//Bulk nodes initialization
	for(int iX=0;iX<NX;iX++)
		for(int iY=0;iY<NY;iY++)
		{


			double laplace_temp=0.0;
			double gradx_temp=0.0;
			double grady_temp=0.0;
			for(int k=0;k<9;k++)
			{
				int iX2=(iX+cx[k]+NX) % NX;
				int iY2=(iY+cy[k]+NY) % NY;
				laplace_temp+=laplacestencil[k]*phase[iX2][iY2];
				gradx_temp+=gradstencilx[k]*phase[iX2][iY2];
				grady_temp+=gradstencily[k]*phase[iX2][iY2];
			}


			//Initialization of the macroscopic fields
			rho[iX][iY]=1.0;
			ux[iX][iY]=0.0;
			uy[iX][iY]=0.0;


			double phase_temp=phase[iX][iY];
			double dense_temp=rho[iX][iY];
			double ux_temp=ux[iX][iY];
			double uy_temp=uy[iX][iY];

			double feq;
			double geq;
			double sum=0.0;
			double sum_phase=0.0;
			double phase_square=phase_temp*phase_temp;
			double pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
			double chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);

			for (int k=1; k<9; k++)
			{
				feq=weights[k]*(3.0*pressure_bulk+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+2.0*ux_temp*uy_temp*cx[k]*cy[k]))
					+kconst*(wxx[k]*gradx_temp*gradx_temp+wyy[k]*grady_temp*grady_temp+wxy[k]*gradx_temp*grady_temp);
				geq=weights[k]*(3.0*chemical_pot+3.0*phase_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*phase_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp
												 +2.0*ux_temp*uy_temp*cx[k]*cy[k]));
				sum+=feq;
				sum_phase+=geq;

				f[iX][iY][k]=feq;
				g[iX][iY][k]=geq;
			}
			f[iX][iY][0]=dense_temp-sum;
			g[iX][iY][0]=phase_temp-sum_phase;

		}

}


void update_bounce_back()
{
	//BB nodes density and velocity specification
	for(int iX=0;iX<NX;iX++)
	{
		int iXtop=(iX+1+NX)%NX;
		int iXbottom=(iX-1+NX)%NX;

		f2[iX][0][2]=f2[iX][1][4];
		f2[iX][0][5]=f2[iXtop][1][7];
		f2[iX][0][6]=f2[iXbottom][1][8];

		f2[iX][NY-1][4]=f2[iX][NY-2][2];
		f2[iX][NY-1][7]=f2[iXbottom][NY-2][5];
		f2[iX][NY-1][8]=f2[iXtop][NY-2][6];


		g2[iX][0][2]=g2[iX][1][4];
		g2[iX][0][5]=g2[iXtop][1][7];
		g2[iX][0][6]=g2[iXbottom][1][8];

		g2[iX][NY-1][4]=g2[iX][NY-2][2];
		g2[iX][NY-1][7]=g2[iXbottom][NY-2][5];
		g2[iX][NY-1][8]=g2[iXtop][NY-2][6];


		rho[iX][0]=1.0;
		rho[iX][NY-1]=1.0;
		ux[iX][0]=0.0;ux[iX][NY-1]=0.0;
		uy[iX][0]=0.0;uy[iX][NY-1]=0.0;
	}

}


void collide_bulk(int counter)
{
    //The phase field should be calculated prior the laplacians
    for(int iX=0;iX<NX;iX++)
        for(int iY=1;iY<NY-1;iY++)
		{
            phase[iX][iY]=0.0;
            for(int iPop=0;iPop<9;iPop++)
   				phase[iX][iY]+=g[iX][iY][iPop];
		}
	for(int iX=0;iX<NX;iX++)
	{
		phase[iX][0]=phase[iX][1]-wall_gradient;
		phase[iX][NY-1]=phase[iX][NY-2]-wall_gradient;
	}
	
    for(int iX=0;iX<NX;iX++)
    {
		    for(int iY=1;iY<NY-1;iY++)
			{

				//Construction equilibrium
				rho[iX][iY]=0.0;
				ux[iX][iY]=0.0;
				uy[iX][iY]=0.0;

				for(int iPop=0;iPop<9;iPop++)
				{
					rho[iX][iY]+=f[iX][iY][iPop];
					ux[iX][iY]+=f[iX][iY][iPop]*cx[iPop];
					uy[iX][iY]+=f[iX][iY][iPop]*cy[iPop];
				}

				ux[iX][iY]=ux[iX][iY]/rho[iX][iY];
				uy[iX][iY]=uy[iX][iY]/rho[iX][iY];


				double laplace_temp=0.0;
				double gradx_temp=0.0;
				double grady_temp=0.0;
				for(int k=0;k<9;k++)
				{
					int iX2=(iX+cx[k]+NX) % NX;
					int iY2=(iY+cy[k]+NY) % NY;
					laplace_temp+=laplacestencil[k]*phase[iX2][iY2];
					gradx_temp+=gradstencilx[k]*phase[iX2][iY2];
					grady_temp+=gradstencily[k]*phase[iX2][iY2];
				}

				double phase_temp=phase[iX][iY];
				double dense_temp=rho[iX][iY];
				double ux_temp=ux[iX][iY];
				double uy_temp=uy[iX][iY];

				double sum=0.0;
				double sum_phase=0.0;
				double phase_square=phase_temp*phase_temp;
				double pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
				double chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);

				double feqeq[9],geqeq[9];

				for (int k=1; k<9; k++)
				{
					feqeq[k]=weights[k]*(3.0*pressure_bulk+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
									+4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+2.0*ux_temp*uy_temp*cx[k]*cy[k]))
					+kconst*(wxx[k]*gradx_temp*gradx_temp+wyy[k]*grady_temp*grady_temp+wxy[k]*gradx_temp*grady_temp);
					geqeq[k]=weights[k]*(3.0*chemical_pot+3.0*phase_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
									+4.5*phase_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+2.0*ux_temp*uy_temp*cx[k]*cy[k]));
					sum+=feqeq[k];
					sum_phase+=geqeq[k];

				}

				feqeq[0]=dense_temp-sum;
				geqeq[0]=phase_temp-sum_phase;

		        double tau_rho=tau_gas+(phase_temp+1.0)/2.0*(tau_liq-tau_gas);
		        omega_rho=1.0/tau_rho;

	

				for(int k=0; k < 9; k++)
				{
					f2[iX][iY][k]=f[iX][iY][k]*(1.0-omega_rho)+omega_rho*feqeq[k];
					g2[iX][iY][k]=g[iX][iY][k]*(1.0-omega_phi)+omega_phi*geqeq[k];
				}


			}
	
	}
}

int main(int argc, char* argv[])
{

    init();


	for(int counter=0;counter<=N;counter++)
	{
        collide_bulk(counter);
     	update_bounce_back();

		//Streaming
		for(int iX=0;iX<NX;iX++)
			for(int iY=1;iY<NY-1;iY++)
				for(int iPop=0;iPop<9;iPop++)
				{
					int iX2=(iX-cx[iPop]+NX)%NX;
					int iY2=(iY-cy[iPop]+NY)%NY;
					f[iX][iY][iPop]=f2[iX2][iY2][iPop];
           	    	g[iX][iY][iPop]=g2[iX2][iY2][iPop];
				}
		
		if (counter%NSIGNAL==0)
		{
			std::cout<<"Time is "<<counter<<"\n";
		}
		//Writing files
		if (counter%NOUTPUT==0)
		{
			if (std::isnan(phase[NX/2][NY/2]))
			{
				std::cout<<"There were some NaN values.\n";
				return 0;
			}
			std::cout<<counter<<"\n";
	
 			std::stringstream filewritephase;
 			std::stringstream filewritevelocityx;
 			std::stringstream filewritevelocityy;
 			std::stringstream counterconvert;
 			counterconvert<<counter;
			filewritephase<<std::fixed;
			filewritevelocityx<<std::fixed;
			filewritevelocityy<<std::fixed;
			

            filewritephase<<"phase"<<std::string(6-counterconvert.str().size(),'0')<<counter;
            filewritevelocityx<<"velocityx"<<std::string(6-counterconvert.str().size(),'0')<<counter;
            filewritevelocityy<<"velocityy"<<std::string(6-counterconvert.str().size(),'0')<<counter;
            

            writephase(filewritephase.str());
            writevelocityx(filewritevelocityx.str());
            writevelocityy(filewritevelocityy.str());
		}


	}
	



	return 0;
}
