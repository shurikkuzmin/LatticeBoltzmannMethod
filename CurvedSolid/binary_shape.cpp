//Simple mass transfer


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>

//Uncomment if you want to clamp tau parameters
//#define CLAMP_TAU

//Domain size
const int NY=301;
const int NX=301;
const int NUM=NX*NY;

//Time steps
int N=20000;
int NOUTPUT=200;
int NSIGNAL=100;

//Other constants
const int NPOP=9;
const int radius=50;
int radius_droplet=10;
double wall_gradient=-0.05;
const int NUM_DROPLETS=3;
const int angle=360/NUM_DROPLETS;
const double pi=4.0*std::atan(1.0);


const double phi_wall=1.5;
const double rho_wall=0.5;
const double force_x=0.0;
const double force_y=0.0;

//BGK relaxation parameter
double omega_rho=1.0;
double omega_phi=1.0;

//Binary-liquid parameters
double aconst=0.04;
double kconst=0.04;
double gammaconst=1.0;
double tau_gas=1.0;
double tau_liq=1.0;




//Fields and populations
double *f;
double *f2;
double *g;
double *g2;
double *phi;
double *rho;
double *ux;
double *uy;
int * geometry;


std::vector<int> bb_nodes;
std::vector<char>* dirs;
std::vector<char> main_dir;

//Underlying lattice parameters
double weights[]={4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};
int cx[]={0,1,0,-1,0,1,-1,-1,1};
int cy[]={0,0,1,0,-1,1,1,-1,-1};
int compliment[]={0,3,4,1,2,7,8,5,6};
double wxx[] = {0.0, 1.0/3.0, -1.0/6.0, 1.0/3.0, -1.0/6.0, -1.0/24.0, -1.0/24.0, -1.0/24.0, -1.0/24.0};
double wyy[] = {0.0, -1.0/6.0, 1.0/3.0, -1.0/6.0, 1.0/3.0, -1.0/24.0, -1.0/24.0, -1.0/24.0, -1.0/24.0};
double wxy[] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0/4.0, -1.0/4.0, 1.0/4.0, -1.0/4.0};

double gradstencilx[9]={0.0,4.0/12.0,0.0,-4.0/12.0,0.0,
			          1.0/12.0,-1.0/12.0,-1.0/12.0,1.0/12.0};

double gradstencily[9]={0.0,0.0,4.0/12.0,0.0,-4.0/12.0,
			          1.0/12.0,1.0/12.0,-1.0/12.0,-1.0/12.0};

double laplacestencil[9]={-20.0/6.0,4.0/6.0,4.0/6.0,4.0/6.0,4.0/6.0,
					   1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0};

int symmetricx[NPOP];
int symmetricy[NPOP];
int symmetricxy_left[NPOP];
int symmetricxy_right[NPOP];

void writephase(std::string const & fname)
{

	std::string filename=fname+".dat";
	std::ofstream fout(filename.c_str());
	fout.precision(10);

	for (int iY=0; iY<NY; iY++)
	{
		for (int iX=0; iX<NX; ++iX)	
		{
			int counter=iY*NX+iX;
			fout<<phi[counter]<<" ";
		}
		fout<<"\n";
	}
}

void writegeometry(std::string const & fname)
{
	std::string filename=fname+".dat";
	std::ofstream fout(filename.c_str());
	fout.precision(10);

	for (int iY=0; iY<NY; iY++)
	{
		for (int iX=0; iX<NX; ++iX)	
		{
			int counter=iY*NX+iX;
			fout<<geometry[counter]<<" ";
		}
		fout<<"\n";
	}

}
void writedensity(std::string const & fname)
{
	std::string filename=fname+".dat";
	std::ofstream fout(filename.c_str());
	fout.precision(10);

	for (int iY=0; iY<NY; iY++)
	{
		for (int iX=0; iX<NX; ++iX)	
		{
			int counter=iY*NX+iX;
			fout<<rho[counter]<<" ";
		}
		fout<<"\n";
	}

}

void writevelocityx(std::string const & fname)
{
	std::string filename=fname+".dat";
	std::ofstream fout(filename.c_str());
	fout.precision(10);

	for (int iY=0; iY<NY; iY++)
	{
		for (int iX=0; iX<NX; ++iX)	
		{
			int counter=iY*NX+iX;
			fout<<ux[counter]<<" ";
		}
		fout<<"\n";
	}

}

void writevelocityy(std::string const & fname)
{
	std::string filename=fname+".dat";
	std::ofstream fout(filename.c_str());
	fout.precision(10);

	for (int iY=0; iY<NY; iY++)
	{
		for (int iX=0; iX<NX; ++iX)	
		{
			int counter=iY*NX+iX;
			fout<<uy[counter]<<" ";
		}
		fout<<"\n";
	}

}

void writevtk(std::string const & fname)
{
	std::string filename=fname+".vtk";
	std::ofstream fout(filename.c_str());
    fout<<"# vtk DataFile Version 3.0\n";
    fout<<"Binary liquid phase field vtk representation\n";
    fout<<"ASCII\n\n";
    fout<<"DATASET STRUCTURED_GRID\n";
    fout<<"DIMENSIONS "<<NX<<" "<<NY<<" "<<1<<"\n";
    fout<<"POINTS "<<NX*NY*1<<" double\n";
    for(int counter=0;counter<NUM;counter++)
    {
    	int iX=counter%NX;
    	int iY=counter/NX;
        fout<<iX<<" "<<iY<<" "<<0<<"\n";
    }
    fout<<"\n";
    fout<<"POINT_DATA "<<NX*NY*1<<"\n";
    
    fout<<"SCALARS phase double\n";
    fout<<"LOOKUP_TABLE phase_table\n";
    for(int counter=0;counter<NUM;counter++)
   		fout<<phi[counter]<<"\n";
    fout<<"SCALARS density double\n";
    fout<<"LOOKUP_TABLE density_table\n";
    for(int counter=0;counter<NUM;counter++)
   		fout<<rho[counter]<<"\n";
 
  	fout<<"SCALARS geometry double\n";
    fout<<"LOOKUP_TABLE geometry_table\n";
    for(int counter=0;counter<NUM;counter++)
    	fout<<geometry[counter]<<"\n";
    fout<<"VECTORS velocity double\n";
    for(int counter=0;counter<NUM;counter++)
    		fout<<ux[counter]<<" "<<uy[counter]<<" 0\n";

    fout.close();
}


void init()
{
	//Creating arrays
	f=new double[NUM*NPOP];
	f2=new double[NUM*NPOP];
	g=new double[NUM*NPOP];
	g2=new double[NUM*NPOP];
	
	//Bulk nodes initialization
	double feq;
	double geq;
	
	for(int iY=0;iY<NY;iY++)
		for(int iX=0;iX<NX;iX++)
		{
			int  counter=iY*NX+iX;
			double dense_temp=rho[counter];
			double ux_temp=ux[counter];
			double uy_temp=uy[counter];
            double phase_temp=phi[counter];
            for (int k=0; k<NPOP; k++)
			{
				feq=weights[k]*(dense_temp+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp
								                +(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp
								                +2.0*ux_temp*uy_temp*cx[k]*cy[k]));
				geq=weights[k]*(phase_temp+3.0*phase_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*phase_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp
								                +(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp
								                +2.0*ux_temp*uy_temp*cx[k]*cy[k]));
								                
								                
								                
				f[counter*NPOP+k]=feq;
				g[counter*NPOP+k]=geq;
			}
			
		}

}

void collide_bulk()
{

    for(int counter=0;counter<NUM;counter++)
    {
		if (geometry[counter]!=1)
			continue;
			
		int iX=counter%NX;
		int iY=counter/NX;
		
		double dense_temp=0.0;
        double phase_temp=phi[counter];
		double ux_temp=0.0;
		double uy_temp=0.0;
		
		for(int k=0;k<NPOP;k++)
		{
			dense_temp+=f[counter*NPOP+k];
			ux_temp+=f[counter*NPOP+k]*cx[k];
			uy_temp+=f[counter*NPOP+k]*cy[k];
		}
	
		ux_temp=(ux_temp+0.5*force_x)/dense_temp;
		uy_temp=(uy_temp+0.5*force_y)/dense_temp;
	
		rho[counter]=dense_temp;
		ux[counter]=ux_temp;
		uy[counter]=uy_temp;

		double laplace_temp=0.0;
		double gradx_temp=0.0;
		double grady_temp=0.0;
		
		
		for(int k=0;k<9;k++)
		{
			int iX2=(iX+cx[k]+NX) % NX;
			int iY2=(iY+cy[k]+NY) % NY;
			int counter2=iY2*NX+iX2;
			laplace_temp+=laplacestencil[k]*phi[counter2];
			gradx_temp+=gradstencilx[k]*phi[counter2];
			grady_temp+=gradstencily[k]*phi[counter2];
		}

		double sum=0.0;
		double sum_phase=0.0;
		double phase_square=phase_temp*phase_temp;
		double pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
		double chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);

		double feqeq[9],geqeq[9],force[9];

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
        #ifdef CLAMP_TAU
        if (phase_temp>1.0)
        	tau_rho=tau_liq;
        if (phase_temp<-1.0)
           	tau_rho=tau_gas;
        #endif   	
        omega_rho=1.0/tau_rho;

		//Obtain force population
        for (int k=0;k<9;k++)
        {
			force[k]=weights[k]*(1.0-0.5*omega_rho)*(3.0*force_x*cx[k]+3.0*force_y*cy[k]+
        	                9.0*((cx[k]*cx[k]-1.0/3.0)*force_x*ux_temp+cx[k]*cy[k]*(force_x*uy_temp+force_y*ux_temp)+
 							(cy[k]*cy[k]-1.0/3.0)*force_y*uy_temp));
        }


		
		for(int k=0; k < NPOP; k++)
		{
			f2[counter*NPOP+k]=f[counter*NPOP+k]*(1.0-omega_rho)+omega_rho*feqeq[k];    	
            g2[counter*NPOP+k]=g[counter*NPOP+k]*(1.0-omega_phi)+omega_phi*geqeq[k];
        }
    }

}

void update_bounce_back()
{
	for(int counter=0;counter<bb_nodes.size();counter++)
	{
		
		//Perform pure bounce back for both phases, probably need to do something else
		for(int k=0;k<dirs[counter].size();k++)
		{
			int dir=dirs[counter][k];
			int counter2=bb_nodes[counter]+cy[dir]*NX+cx[dir];
			f2[bb_nodes[counter]*NPOP+dir]=f2[counter2*NPOP+compliment[dir]];
			g2[bb_nodes[counter]*NPOP+dir]=g2[counter2*NPOP+compliment[dir]];
		}

		//int dir=main_dir[counter];
		//int counter2=bb_nodes[counter]+cy[main_dir[counter]]*NX+cx[main_dir[counter]];
		//f2[bb_nodes[counter]*NPOP+k]=f2[counter2*NPOP+compliment()]
		//if ( (dir==1) || (dir==3))
		//	for(int k=0;k<NPOP;k++)
		//		f2[bb_nodes[counter]*NPOP+k]=f2[counter2*NPOP+symmetricx[k]];
	    
		//else if ( (dir==2) || (dir==4))
		//	for(int k=0;k<NPOP;k++)
		//		f2[bb_nodes[counter]*NPOP+k]=f2[counter2*NPOP+symmetricy[k]];
		//else if ( (dir==5) || (dir==7))
		//	for(int k=0;k<NPOP;k++)
		//		f2[bb_nodes[counter]*NPOP+k]=f2[counter2*NPOP+symmetricxy_right[k]];
		//else 
		//	for(int k=0;k<NPOP;k++)
		//		f2[bb_nodes[counter]*NPOP+k]=f2[counter2*NPOP+symmetricxy_left[k]];

		//density and velocity specification
		
		//double dense_temp=0.0;
		//double ux_temp=0.0;
		//double uy_temp=0.0;
		
		//for(int k=0;k<NPOP;k++)
		//{
		//	dense_temp=dense_temp+f2[bb_nodes[counter]*NPOP+k];
		//	ux_temp=ux_temp+f2[bb_nodes[counter]*NPOP+k]*cx[k];
		//	uy_temp=uy_temp+f2[bb_nodes[counter]*NPOP+k]*cy[k];
		//}
		//ux_temp=ux_temp/dense_temp;
		//uy_temp=uy_temp/dense_temp;
		
		//rho[bb_nodes[counter]]=dense_temp;
		//ux[bb_nodes[counter]]=ux_temp;
		//uy[bb_nodes[counter]]=uy_temp;
	}
	
}

void update_phase()
{
	//Phase need to be calculated before laplacians
    double phase_temp;
	for(int counter=0;counter<NUM;counter++)
	{
        phase_temp=0.0;
        if (geometry[counter]==1)
        {
        	for(int k=0;k<NPOP;k++)
        		phase_temp+=g[counter*NPOP+k];
        	phi[counter]=phase_temp;
		}
	}

	for(int counter=0;counter<bb_nodes.size();counter++)
	{
		int dir=main_dir[counter];
		int counter2=bb_nodes[counter]+cy[dir]*NX+cx[dir];
		phi[bb_nodes[counter]]=phi[counter2]-wall_gradient;
	}
	
}



void initialize_geometry()
{
    geometry=new int[NUM];
    rho=new double[NUM];
    ux=new double[NUM];
    uy=new double[NUM];
    phi=new double[NUM];
    

	//Solid nodes 
	for(int counter=0;counter<NUM;counter++)
	{
		int iX=counter%NX;
		int iY=counter/NX;
	    
	    if ((iX-(NX-1)/2)*(iX-(NX-1)/2)+(iY-(NY-1)/2)*(iY-(NY-1)/2)<radius*radius)
	    	geometry[counter]=-1;
	    else
	    	geometry[counter]=1;  
	}

	for(int counter=0;counter<NUM;counter++)
	{
	    bool flag=false;
	    if (geometry[counter]==-1)
	    {
	    	int iX=counter%NX;
	    	int iY=counter/NX;
	    	for (int iPop=1;iPop<NPOP;iPop++)
	    	{
				int iX2=(iX+cx[iPop]+NX)%NX;
			 	int iY2=(iY+cy[iPop]+NY)%NY;
	    		int counter2=iY2*NX+iX2;
	    		if (geometry[counter2]==1)
	    			flag=true;
	    	}
	    }
	    if (flag)
	    	geometry[counter]=0;
	}

	//Initialization of density
    for(int counter=0;counter<NUM;counter++)
    {
		ux[counter]=0.0;
		uy[counter]=0.0;

		if (geometry[counter]==0)
		{
		    rho[counter]=rho_wall;
			bb_nodes.push_back(counter);
		}
		else if(geometry[counter]==-1)
			rho[counter]=rho_wall;
		else
			rho[counter]=1.0;
	}
	
	//Identifying fluids
	for(int counter=0;counter<NUM;counter++)
	{
		int iX=counter%NX;
		int iY=counter/NX;
		
		bool flag_droplet=false;
		for (int counter_droplet=0;counter_droplet<NUM_DROPLETS;counter_droplet++)
		{	
			double center_x=(NX-1)/2-std::sin(angle*(counter_droplet-1)*pi/180.0)*(radius_droplet+radius);
			double center_y=(NY-1)/2+std::cos(angle*(counter_droplet-1)*pi/180.0)*(radius_droplet+radius);
			if (((iX-center_x)*(iX-center_x)+(iY-center_y)*(iY-center_y)<=radius_droplet*radius_droplet) && (geometry[counter]==1))
		    {	
		    	phi[counter]=1.0;
		    	flag_droplet=true;
		    }
		}
		if (!flag_droplet)
		{
			if (geometry[counter]==1)
				phi[counter]=-1.0;
			else
				phi[counter]=phi_wall;
		
		}
	}
	
	symmetricx[0]=0;
    for(int k=1;k<NPOP;k++)
    	for(int l=1;l<NPOP;l++)
    	    if ((-cx[k]==cx[l])&&(cy[k]==cy[l]))
            {
                symmetricx[k]=l;
                break;
            }
    
     symmetricy[0]=0;
     for(int k=1;k<NPOP;k++)
        for(int l=1;l<NPOP;l++)
            if ((cx[k]==cx[l])&&(cy[k]==-cy[l]))
            {
           		symmetricy[k]=l;
                break;
            }
     
     symmetricxy_left[0]=0;
     symmetricxy_left[1]=2;
     symmetricxy_left[2]=1;
     symmetricxy_left[3]=4;
     symmetricxy_left[4]=3;
     symmetricxy_left[5]=5;
     symmetricxy_left[6]=8;
     symmetricxy_left[7]=7;
     symmetricxy_left[8]=6;
           
     symmetricxy_right[0]=0;
     symmetricxy_right[1]=4;
     symmetricxy_right[2]=3;
     symmetricxy_right[3]=2;
     symmetricxy_right[4]=1;
     symmetricxy_right[5]=7;
     symmetricxy_right[6]=6;
     symmetricxy_right[7]=5;
     symmetricxy_right[8]=8;

     
	//Finding directions for BB nodes
    dirs=new std::vector<char>[bb_nodes.size()];
    for(int counter=0;counter<bb_nodes.size();counter++)
	{
		for(int k=1;k<NPOP;k++)
		{
			int counter2=bb_nodes[counter]+cy[k]*NX+cx[k];
			if (geometry[counter2]==1)
				dirs[counter].push_back(k);
		}
	}
	

	for(int counter=0;counter<bb_nodes.size();counter++)
	{
     	int nx=0;
     	int ny=0;
		bool flag=false;
     	for(int k=1;k<5;k++)
		{
			int counter2=bb_nodes[counter]+cy[k]*NX+cx[k];
			if (geometry[counter2]==1)
			{
				flag=true;
				nx=nx+cx[k];
				ny=ny+cy[k];
			}
			
		}
		if (!flag)
			for(int k=5;k<NPOP;k++)
			{
				int counter2=bb_nodes[counter]+cy[k]*NX+cx[k];
				if (geometry[counter2]==1)
				{
					flag=true;
					nx=nx+cx[k];
					ny=ny+cy[k];
				}
			}
		
		for(int k=1;k<NPOP;k++)
			if ((nx==cx[k])&&(ny==cy[k]))
			{
				main_dir.push_back(k);
			}
			
	}
	std::cout<<"BB size="<<bb_nodes.size()<<"\n";
	std::cout<<"Main size="<<main_dir.size()<<"\n";
	writephase("phase");
	writegeometry("geometry");
}


void finish_simulation()
{
	delete[] geometry;
	delete[] rho;
	delete[] ux;
	delete[] uy;
	delete[] f;
	delete[] f2;
	delete[] dirs;
	delete[] phi;
	delete[] g;
	delete[] g2;
}

void stream()
{
    for(int counter=0;counter<NUM;counter++)
	{
		if (geometry[counter]!=1)
			continue;
		int iX=counter%NX;
		int iY=counter/NX;

		for(int iPop=0;iPop<NPOP;iPop++)
		{
			int iX2=(iX-cx[iPop]+NX)%NX;
			int iY2=(iY-cy[iPop]+NY)%NY;
			int counter2=iY2*NX+iX2;
			f[counter*NPOP+iPop]=f2[counter2*NPOP+iPop];
			g[counter*NPOP+iPop]=g2[counter2*NPOP+iPop];
		}
	}

	
}

void calculate_mass()
{
	double mass=0.0;
	for(int counter=0;counter<NUM;counter++)
		if (geometry[counter]==1)
			mass+=phi[counter];
	std::cout<<"Mass is "<<mass<<"\n";
}

int main(int argc, char* argv[])
{

    if (argc==3)
    {
    	radius_droplet=atoi(argv[1]);
    	wall_gradient=atof(argv[2]);
	}
	std::cout<<"Radius of the droplet is "<<radius_droplet<<"\n";
	std::cout<<"Wall gradient is "<<wall_gradient<<"\n";
	
    initialize_geometry();
    init();
      

	for(int counter=0;counter<=N;counter++)
	{
        update_phase();
        collide_bulk();
        update_bounce_back();
		stream();
        
	    if (counter%NSIGNAL==0)
	    {
	    	std::cout<<"Counter="<<counter<<"\n";
    		calculate_mass();
    	}
		//Writing files
		if (counter%NOUTPUT==0)
		{
			std::stringstream filewritedensity;
  			std::stringstream filewritevelocityx;
  			std::stringstream filewritevelocityy;
 			std::stringstream filevtk;
 			
 			std::stringstream counterconvert;
 			counterconvert<<counter;
 			filewritedensity<<std::fixed;
 			filewritevelocityx<<std::fixed;
 			filewritevelocityy<<std::fixed;
 			filevtk<<std::fixed;

			//filewritedensity<<"den"<<std::string(7-counterconvert.str().size(),'0')<<counter;
			//filewritevelocityx<<"velx"<<std::string(7-counterconvert.str().size(),'0')<<counter;
			//filewritevelocityy<<"vely"<<std::string(7-counterconvert.str().size(),'0')<<counter;
			filevtk<<"vtk"<<std::string(7-counterconvert.str().size(),'0')<<counter;
			
 			writedensity(filewritedensity.str());
 			writevelocityx(filewritevelocityx.str());
 			writevelocityy(filewritevelocityy.str());
 			writevtk(filevtk.str());
		}

	}

    finish_simulation();
   
   	return 0;
}
