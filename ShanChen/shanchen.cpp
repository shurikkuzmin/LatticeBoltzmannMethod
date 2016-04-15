#include <iostream>
#include <cmath>
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>

//Parameters of the grid
const int nx=128;
const int ny=128;
const int npop=9;
const int nsteps=30000;
const int noutput=1000;

//Parameters of the LBM
const int cx[]={0,1,0,-1,0,1,-1,-1,1};
const int cy[]={0,0,1,0,-1,1,1,-1,-1};
const double weights[]={4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};
double tau=1.0;

//Parameters of the Shan-Chen model
double rhol=1.95;
double rhog=0.15;
int radius=20;
double g=-5.0;
//Arrays
double rho[nx*ny];
double u1[nx*ny];
double u2[nx*ny];
double f_mem[nx*ny*npop];
double f2_mem[nx*ny*npop];
double feq[npop];

void writeFile(std::string name, double* f, int n)
{
	name="tmp/"+name;
	std::ofstream fout(name.c_str());
	for(int counterX=0;counterX<nx;counterX++)
	{	
		for(int counterY=0;counterY<ny;counterY++)
			fout<<f[counterX*ny+counterY]<<" ";
		fout<<"\n";
	}		
	fout<<std::endl;
}

void calculateDensities()
{
	double rho_gas=0.0;
	double rho_liq=0.0;
	
	for(int counterX=0;counterX<10;counterX++)
		for(int counterY=0;counterY<10;counterY++)
			rho_gas+=rho[counterX*ny+counterY];
	
	for(int counterX=nx/2-5;counterX<nx/2+5;counterX++)
		for(int counterY=ny/2-5;counterY<ny/2+5;counterY++)
			rho_liq+=rho[counterX*ny+counterY];

	rho_gas=rho_gas/100.0;
	rho_liq=rho_liq/100.0;
	
	double pressure_liq,pressure_gas;
	pressure_liq=rho_liq/3.0+g/6.0*(1.0-exp(-rho_liq))*(1.0-exp(-rho_liq));
	pressure_gas=rho_gas/3.0+g/6.0*(1.0-exp(-rho_gas))*(1.0-exp(-rho_gas));
	std::cout<<"Calculated densities are: Rho_gas="<<rho_gas<<" Rho_liq="<<rho_liq<<"\n";
	std::cout<<"Delta Pressure is:"<<pressure_liq-pressure_gas<<"\n";
}


int main(int argc, char** argv)
{
	
	if (argc==3)
	{
		tau=atof(argv[2]);
		radius=atoi(argv[3]);
	}

	//Memory preparation
	double* f=f_mem;
	double* f2=f2_mem;
	
	//Initialization
	for (int i=0; i < nx*ny; i++)
	{	
		if ((i/nx - ny/2.0)*(i/nx - ny/2.0) + (i%nx - nx/2.0)*(i%nx - nx/2.0) <= radius * radius)
		{
			rho[i]=rhol;
		}
		else 
			rho[i]=rhog;

		double dense,v1,v2;
		
		dense=rho[i];
		v1=v2=u1[i]=u2[i]=0.0;
		float usq = v1*v1 + v2*v2;
		feq[0] = 4.0/9.0 * dense * (1.0 - 1.5 * usq); 
		feq[1] = 1.0/9.0 * dense * (1.0 + 3*v1 + 4.5*v1*v1 - 1.5*usq); 
		feq[2] = 1.0/9.0 * dense * (1.0 + 3*v2 + 4.5*v2*v2 - 1.5*usq); 
		feq[3] = 1.0/9.0 * dense * (1.0 - 3*v1 + 4.5*v1*v1 - 1.5*usq); 
		feq[4] = 1.0/9.0 * dense * (1.0 - 3*v2 + 4.5*v2*v2 - 1.5*usq); 
		feq[5] = 1.0/36.0 * dense * (1.0 + 3*(v1 + v2) + 4.5*(v1 + v2)*(v1 + v2) - 1.5*usq); 
		feq[6] = 1.0/36.0 * dense * (1.0 + 3*(-v1 + v2) + 4.5*(-v1 + v2)*(-v1 + v2) - 1.5*usq);
		feq[7] = 1.0/36.0 * dense * (1.0 + 3*(-v1 - v2) + 4.5*(v1 + v2)*(v1 + v2) - 1.5*usq); 
		feq[8] = 1.0/36.0 * dense * (1.0 + 3*(v1 - v2) + 4.5*(v1 - v2)*(v1 -v2) - 1.5*usq); 
		for (int k=0; k<npop; k++) {
			f[9*i+k]=feq[k];
			f2[9*i+k]=feq[k];
		}
	}
	

	time_t start, finish;
	start = time(NULL);
	//Main loop
	for (int timecounter=0; timecounter<=nsteps;timecounter++) 
	{
		
		//Calculation of the density field
		for (int i=0; i<nx*ny; i++) 
		{
			rho[i]=0; 
			for (int k=0; k<9; k++ )
			{			
				rho[i]+=f[9*i+k]; 
			}		
			
		}

		
		//Collision and streaming
		for (int iY=0; iY<ny; iY++) 
		    for(int iX=0;iX<nx;iX++)
		    {
		    
				int i=iY*nx+iX;
				double dense,v1,v2;
		       	
				dense=rho[i];

				float fx=0.0;
				float fy=0.0;

				for(int k=0;k<9;k++)
				{
					int iX2=(iX+cx[k]+nx) % nx; 
					int iY2=(iY+cy[k]+ny) % ny;
					fx+=weights[k]*cx[k]*(1.0-exp(-rho[nx*iY2+iX2]));
					fy+=weights[k]*cy[k]*(1.0-exp(-rho[nx*iY2+iX2]));
				}
			
				fx=-g*(1.0-exp(-rho[i]))*fx;
				fy=-g*(1.0-exp(-rho[i]))*fy;
			
				v1=u1[i]=(f[9*i+1]-f[9*i+3]+f[9*i+5]-f[9*i+6]-f[9*i+7]+f[9*i+8])/dense+fx/(2.0*dense); 
				v2=u2[i]=(f[9*i+2]-f[9*i+4]+f[9*i+5]+f[9*i+6]-f[9*i+7]-f[9*i+8])/dense+fy/(2.0*dense); 
			
				double fpop[9];
				for(int k=0;k<9;k++)
					fpop[k]=weights[k]*(1-0.5/tau)*((3*(cx[k]-v1)+9*cx[k]*(cx[k]*v1+cy[k]*v2))*fx
		               +(3*(cy[k]-v2)+9*cy[k]*(cx[k]*v1+cy[k]*v2))*fy);
			
				float usq = v1*v1 + v2*v2;	
			
				feq[0] = 4.0/9.0 * dense * (1.0 - 1.5 * usq); 
				feq[1] = 1.0/9.0 * dense * (1.0 + 3*v1 + 4.5*v1*v1 - 1.5*usq); 
				feq[2] = 1.0/9.0 * dense * (1.0 + 3*v2 + 4.5*v2*v2 - 1.5*usq); 
				feq[3] = 1.0/9.0 * dense * (1.0 - 3*v1 + 4.5*v1*v1 - 1.5*usq); 
				feq[4] = 1.0/9.0 * dense * (1.0 - 3*v2 + 4.5*v2*v2 - 1.5*usq); 
				feq[5] = 1.0/36.0 * dense * (1.0 + 3*(v1 + v2) + 4.5*(v1 + v2)*(v1 + v2) - 1.5*usq); 
				feq[6] = 1.0/36.0 * dense * (1.0 + 3*(-v1 + v2) + 4.5*(-v1 + v2)*(-v1 + v2) - 1.5*usq);
				feq[7] = 1.0/36.0 * dense * (1.0 + 3*(-v1 - v2) + 4.5*(v1 + v2)*(v1 + v2) - 1.5*usq);
				feq[8] = 1.0/36.0 * dense * (1.0 + 3*(v1 - v2) + 4.5*(v1 - v2)*(v1 -v2) - 1.5*usq);
			
				for(int k=0; k<9; k++) 
				{  
					int iX2=(iX+cx[k]+nx) % nx; 
					int iY2=(iY+cy[k]+nx) % nx;
					f[9*i+k]+=-1.0/tau*(f[9*i+k]-feq[k])+fpop[k]; 
					f2[9*(nx*iY2+iX2)+k]=f[9*i+k]; 
				}  
			}
				
			std::swap(f,f2);
			//Preparation of the output
			if (timecounter%noutput==0)
			{
				std::stringstream denstream;
				std::stringstream velxstream;
				std::stringstream velystream;
				std::stringstream len;
				len<<timecounter;
		
				denstream << "density"<<std::string(5-len.str().size(),'0')<<timecounter<<".dat";
				velxstream << "velx"<<std::string(5-len.str().size(),'0')<<timecounter<<".dat";
				velystream << "vely"<<std::string(5-len.str().size(),'0')<<timecounter<<".dat";
				
				writeFile(denstream.str(), rho, nx*ny);
				writeFile(velxstream.str(), u1, nx*ny);
				writeFile(velystream.str(), u2, nx*ny);

				std::cout<<"Iteration is "<<timecounter<<"\n";
			}
		
	}

	
	finish = time(NULL);
	
	std::cout<<"Overall time is "<<finish-start<<" sec"<<"\n";
	
    return 0;
}
