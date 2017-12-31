//The code to simulate the deposition: We first initialize the hydrofields and the impose equilibrium phase values.
//Certain flags are indicated

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>

//Uncomment if you want to clamp tau parameters
//#define CLAMP_TAU

//Uncomment it if you want to read the file
//Use it with caution for deposition
#define READ_FILE

//Domain size
int NX,NY,NUM;

//Time steps
int N=100000;
int NHYDRO=10000;
int NOUTPUT=1000;
int NSIGNAL=100;

//Other constants
const int NPOP=9;
const int radius=30;
int radius_droplet=10;
double wall_gradient=-0.35;
double wall_gradient_boundary=0.0;
const int NUM_DROPLETS=3;
const int angle=360/NUM_DROPLETS;
const double pi=4.0*std::atan(1.0);


const double phi_wall=1.5;
const double rho_wall=0.5;

double force_x=0.0;
double force_y=0.0;
double force_x_hydro=0.00001;
double force_y_hydro=0.0;

//BGK relaxation parameter
double omega_rho=1.0;
double omega_phi=1.0;
double omega_hydro=1.0;

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
double *h;
double *h2;
double *phi;
double *rho_f;
double *rho_h;
double *ux_f;
double *ux_h;
double *uy_f;
double *uy_h;
int *geometry;


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
			fout<<rho_f[counter]<<" ";
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
			fout<<ux_f[counter]<<" ";
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
			fout<<uy_f[counter]<<" ";
		}
		fout<<"\n";
	}

}
//Did modifications to add another field
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
   		fout<<rho_f[counter]<<"\n";
 
  	fout<<"SCALARS geometry double\n";
    fout<<"LOOKUP_TABLE geometry_table\n";
    for(int counter=0;counter<NUM;counter++)
    	fout<<geometry[counter]<<"\n";
    fout<<"VECTORS velocity double\n";
    for(int counter=0;counter<NUM;counter++)
    		fout<<ux_f[counter]<<" "<<uy_f[counter]<<" 0\n";
    fout<<"SCALARS density_one_component double\n";
    fout<<"LOOKUP_TABLE density_one_component_table\n";
    for(int counter=0;counter<NUM;counter++)
        fout<<rho_h[counter]<<"\n";
    fout<<"VECTORS velocity_one_component double\n";
    for(int counter=0;counter<NUM;counter++)
        fout<<ux_h[counter]<<" "<<uy_h[counter]<<" 0\n";
    
    fout.close();
        
}


void init()
{
	//Creating arrays
	f=new double[NUM*NPOP];
	f2=new double[NUM*NPOP];
	g=new double[NUM*NPOP];
	g2=new double[NUM*NPOP];
	h=new double[NUM*NPOP];
	h2=new double[NUM*NPOP];
    
	//Bulk nodes initialization
	double feq;
	double geq;
	
	for(int iY=0;iY<NY;iY++)
		for(int iX=0;iX<NX;iX++)
		{
			int  counter=iY*NX+iX;
			double dense_temp=rho_f[counter];
			double ux_temp=ux_f[counter];
			double uy_temp=uy_f[counter];
            double phase_temp=phi[counter];
            for (int k=0; k<NPOP; k++)
			{
				feq=weights[k]*(dense_temp+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp
								                +(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp
								                +2.0*ux_temp*uy_temp*cx[k]*cy[k]));
                f[counter*NPOP+k]=feq;
                h[counter*NPOP+k]=feq;
            	geq=weights[k]*(phase_temp+3.0*phase_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*phase_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp
								                +(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp
								                +2.0*ux_temp*uy_temp*cx[k]*cy[k]));	                
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
		
		if ( (iY==0) || (iY==NY-1) )
            continue;
		double dense_temp=0.0;
        double phase_temp=phi[counter];
		double ux_temp=0.0;
		double uy_temp=0.0;
		
        double dense_temp_hydro=0.0;
        double ux_temp_hydro=0.0;
        double uy_temp_hydro=0.0;
        
		for(int k=0;k<NPOP;k++)
		{
			dense_temp+=f[counter*NPOP+k];
			ux_temp+=f[counter*NPOP+k]*cx[k];
			uy_temp+=f[counter*NPOP+k]*cy[k];
			dense_temp_hydro+=h[counter*NPOP+k];
			ux_temp_hydro+=h[counter*NPOP+k]*cx[k];
			uy_temp_hydro+=h[counter*NPOP+k]*cy[k];
		}
	
		ux_temp=(ux_temp+0.5*force_x)/dense_temp;
		uy_temp=(uy_temp+0.5*force_y)/dense_temp;
		ux_temp_hydro=(ux_temp_hydro+0.5*force_x_hydro)/dense_temp_hydro;
		uy_temp_hydro=(uy_temp_hydro+0.5*force_y_hydro)/dense_temp_hydro;
        
	
		rho_f[counter]=dense_temp;
		ux_f[counter]=ux_temp;
		uy_f[counter]=uy_temp;
		
        rho_h[counter]=dense_temp_hydro;
        ux_h[counter]=ux_temp_hydro;
        uy_h[counter]=uy_temp_hydro;
        
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
        
		double sum_phase=0.0;
		double sum=0.0;
        double sum_hydro=0.0;
        
		double phase_square=phase_temp*phase_temp;
		double pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
		double chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);
        double geqeq[9];
      	double feqeq[9],feqeq_hydro[9],force[9],force_hydro[9];
      	
		double pressure_bulk_hydro=dense_temp_hydro/3.0;
		for (int k=1; k<9; k++)
		{
			feqeq_hydro[k]=weights[k]*(3.0*pressure_bulk_hydro+3.0*dense_temp_hydro*(cx[k]*ux_temp_hydro+cy[k]*uy_temp_hydro)
                +4.5*dense_temp_hydro*((cx[k]*cx[k]-1.0/3.0)*ux_temp_hydro*ux_temp_hydro
                +(cy[k]*cy[k]-1.0/3.0)*uy_temp_hydro*uy_temp_hydro+2.0*ux_temp_hydro*uy_temp_hydro*cx[k]*cy[k]));
            sum_hydro+=feqeq_hydro[k];
            
            feqeq[k]=weights[k]*(3.0*pressure_bulk+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
                +4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+2.0*ux_temp*uy_temp*cx[k]*cy[k]))
                +kconst*(wxx[k]*gradx_temp*gradx_temp+wyy[k]*grady_temp*grady_temp+wxy[k]*gradx_temp*grady_temp);
			
			sum+=feqeq[k];
            
			geqeq[k]=weights[k]*(3.0*chemical_pot+3.0*phase_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
							+4.5*phase_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+2.0*ux_temp*uy_temp*cx[k]*cy[k]));
			sum_phase+=geqeq[k];
		}

		feqeq[0]=dense_temp-sum;
        feqeq_hydro[0]=dense_temp_hydro-sum_hydro;
        geqeq[0]=phase_temp-sum_phase;
        
        double tau_rho=tau_gas+(phase_temp+1.0)/2.0*(tau_liq-tau_gas);
        #ifdef CLAMP_TAU
        if (phase_temp>1.0)
            tau_rho=tau_liq;
        if (phase_temp<-1.0)
            tau_rho=tau_gas;
        #endif   	
        omega_rho=1.0/tau_rho;
        omega_hydro=1.0/tau_gas;
		//Obtain force population
        for (int k=0;k<9;k++)
        {
			force[k]=weights[k]*(1.0-0.5*omega_rho)*(3.0*force_x*cx[k]+3.0*force_y*cy[k]+
        	                9.0*((cx[k]*cx[k]-1.0/3.0)*force_x*ux_temp+cx[k]*cy[k]*(force_x*uy_temp+force_y*ux_temp)+
 							(cy[k]*cy[k]-1.0/3.0)*force_y*uy_temp));
    		force_hydro[k]=weights[k]*(1.0-0.5*omega_hydro)*(3.0*force_x_hydro*cx[k]+3.0*force_y_hydro*cy[k]+
        	                9.0*((cx[k]*cx[k]-1.0/3.0)*force_x_hydro*ux_temp_hydro
        	                +cx[k]*cy[k]*(force_x_hydro*uy_temp_hydro+force_y_hydro*ux_temp_hydro)+
    							(cy[k]*cy[k]-1.0/3.0)*force_y_hydro*uy_temp_hydro));

        }
		
		for(int k=0; k < NPOP; k++)
		{
			f2[counter*NPOP+k]=f[counter*NPOP+k]*(1.0-omega_rho)+omega_rho*feqeq[k]+force[k];    	
            h2[counter*NPOP+k]=h[counter*NPOP+k]*(1.0-omega_hydro)+omega_hydro*feqeq_hydro[k]+force_hydro[k];    	
            g2[counter*NPOP+k]=g[counter*NPOP+k]*(1.0-omega_phi)+omega_phi*geqeq[k];
        }
    }

}

void update_bounce_back()
{
	for(int counter=0;counter<bb_nodes.size();counter++)
	{
		
		for(int k=0;k<dirs[counter].size();k++)
		{
			int dir=dirs[counter][k];
			int counter2=bb_nodes[counter]+cy[dir]*NX+cx[dir];
			f2[bb_nodes[counter]*NPOP+dir]=f2[counter2*NPOP+compliment[dir]];
			g2[bb_nodes[counter]*NPOP+dir]=g2[counter2*NPOP+compliment[dir]];
            h2[bb_nodes[counter]*NPOP+dir]=h2[counter2*NPOP+compliment[dir]];
		}
	}
	
	//Perform bounce back on walls
    for(int iX=0;iX<NX;iX++)
    {
        int iX_top=(iX+1+NX)%NX;
        int iX_bottom=(iX-1+NX)%NX;
        
        f2[iX*NPOP+2]=f2[(NX+iX)*NPOP+4];
        f2[iX*NPOP+5]=f2[(NX+iX_top)*NPOP+7];
        f2[iX*NPOP+6]=f2[(NX+iX_bottom)*NPOP+8];
        f2[(NX*(NY-1)+iX)*NPOP+4]=f2[(NX*(NY-2)+iX)*NPOP+2];
        f2[(NX*(NY-1)+iX)*NPOP+7]=f2[(NX*(NY-2)+iX_bottom)*NPOP+5];
        f2[(NX*(NY-1)+iX)*NPOP+8]=f2[(NX*(NY-2)+iX_top)*NPOP+6];
        
        //Phase populations
        g2[iX*NPOP+2]=g2[(NX+iX)*NPOP+4];
        g2[iX*NPOP+5]=g2[(NX+iX_top)*NPOP+7];
        g2[iX*NPOP+6]=g2[(NX+iX_bottom)*NPOP+8];
        g2[(NX*(NY-1)+iX)*NPOP+4]=g2[(NX*(NY-2)+iX)*NPOP+2];
        g2[(NX*(NY-1)+iX)*NPOP+7]=g2[(NX*(NY-2)+iX_bottom)*NPOP+5];
        g2[(NX*(NY-1)+iX)*NPOP+8]=g2[(NX*(NY-2)+iX_top)*NPOP+6];

        //One-component flow
        h2[iX*NPOP+2]=h2[(NX+iX)*NPOP+4];
        h2[iX*NPOP+5]=h2[(NX+iX_top)*NPOP+7];
        h2[iX*NPOP+6]=h2[(NX+iX_bottom)*NPOP+8];
        h2[(NX*(NY-1)+iX)*NPOP+4]=h2[(NX*(NY-2)+iX)*NPOP+2];
        h2[(NX*(NY-1)+iX)*NPOP+7]=h2[(NX*(NY-2)+iX_bottom)*NPOP+5];
        h2[(NX*(NY-1)+iX)*NPOP+8]=h2[(NX*(NY-2)+iX_top)*NPOP+6];


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
	
	//Update boundary walls
    for(int iX=0;iX<NX;iX++)
    {
        phi[iX]=phi[iX+NX]-wall_gradient_boundary;
        phi[NX*(NY-1)+iX]=phi[NX*(NY-2)+iX]-wall_gradient_boundary;
    }
	
}



void initialize_geometry()
{
    NY=201;
    NX=601;
    NUM=NX*NY;
    geometry=new int[NUM];
    rho_f=new double[NUM];
    ux_f=new double[NUM];
    uy_f=new double[NUM];
    rho_h=new double[NUM];
    ux_h=new double[NUM];
    uy_h=new double[NUM];
    
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

	for(int counter=0;counter<NUM;counter++)
		if (geometry[counter]==0)
			bb_nodes.push_back(counter);
    
	
	     
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
	
    //Determination of the normal
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
	writegeometry("geometry");
}

void initialize_hydro()
{
 	//Initialization of density
    for(int counter=0;counter<NUM;counter++)
    {
		ux_f[counter]=0.0;
		uy_f[counter]=0.0;
        ux_h[counter]=0.0;
        uy_h[counter]=0.0;

		if (geometry[counter]==1) 
		{
            rho_f[counter]=1.0;
            rho_h[counter]=1.0;
		}
		else
		{
		    rho_f[counter]=rho_wall;
            rho_h[counter]=rho_wall;
		}	
	}
    //Initialization of density walls
    for(int iX=0;iX<NX;iX++)
    {
        rho_f[iX]=rho_wall;
        rho_f[NX*(NY-1)+iX]=rho_wall;
        rho_h[iX]=rho_wall;
        rho_h[NX*(NY-1)+iX]=rho_wall;
    }
    
}

void initialize_phase()
{
	//Identifying fluids
	for(int counter=0;counter<NUM;counter++)
		if (geometry[counter]==1)
			phi[counter]=-1.0;
		else 
			phi[counter]=phi_wall;

    int xcenter=20;
    int ycenter=(NY-1)/2;

    for(int counter=0;counter<NUM;counter++)
	{
		int iX=counter%NX;
		int iY=counter/NX;
		
		if ((iX-xcenter)*(iX-xcenter)+(iY-ycenter)*(iY-ycenter)<=radius_droplet*radius_droplet)
            phi[counter]=1;
    }    
    writephase("phase");
}

void finish_simulation()
{
	delete[] geometry;
	delete[] rho_f;
	delete[] ux_f;
	delete[] uy_f;
	delete[] rho_h;
	delete[] ux_h;
	delete[] uy_h;	
	delete[] f;
	delete[] f2;
	delete[] dirs;
	delete[] phi;
	delete[] g;
	delete[] g2;
    delete[] h;
    delete[] h2;
}

void stream()
{
    for(int counter=0;counter<NUM;counter++)
	{
		if (geometry[counter]!=1)
			continue;
		int iX=counter%NX;
		int iY=counter/NX;

        if ((iY==0) || (iY==NY-1))
            continue;
		for(int iPop=0;iPop<NPOP;iPop++)
		{
			int iX2=(iX-cx[iPop]+NX)%NX;
			int iY2=(iY-cy[iPop]+NY)%NY;
			int counter2=iY2*NX+iX2;
			f[counter*NPOP+iPop]=f2[counter2*NPOP+iPop];
            h[counter*NPOP+iPop]=h2[counter2*NPOP+iPop];
			g[counter*NPOP+iPop]=g2[counter2*NPOP+iPop];
		}
	}

	
}

void calculate_mass()
{
    double mass;
    mass=0.0;
    for(int counter=0;counter<NUM;counter++)
	{   
        int iY=counter/NX;
		if ( (geometry[counter]==1) && (iY!=0) && (iY!=NY-1) )
			mass+=rho_f[counter];	
	}	
	std::cout<<"Mass_density is "<<mass<<"\n";
    mass=0.0;
    for(int counter=0;counter<NUM;counter++)
	{   
        int iY=counter/NX;
		if ( (geometry[counter]==1) && (iY!=0) && (iY!=NY-1) )
			mass+=phi[counter];	
	}
    std::cout<<"Mass_phase is "<<mass<<"\n";
}

int main(int argc, char* argv[])
{

    if (argc==4)
    {
    	radius_droplet=atoi(argv[1]);
    	wall_gradient=atof(argv[2]);
        force_x_hydro=atof(argv[3]);
	}
	std::cout<<"Radius of the droplet is "<<radius_droplet<<"\n";
	std::cout<<"Wall gradient is "<<wall_gradient<<"\n";
	
	initialize_geometry();
    initialize_hydro();
    initialize_phase();
    init();
    bool flag_hydro=true;
	for(int counter=0;counter<=N;counter++)
	{
        if ((counter>NHYDRO) && (flag_hydro))
        {
            force_x=force_x_hydro;
            flag_hydro=false;
            std::cout<<"Size of the array"<<sizeof(f)<<"\n";
            memcpy(f,h,sizeof(f));
            memcpy(f2,h2,sizeof(f2));
            memcpy(rho_f,rho_h,sizeof(rho_f));
            memcpy(ux_f,ux_h,sizeof(ux_f));
            memcpy(uy_f,uy_h,sizeof(uy_f));
        }
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
            std::cout<<"Output "<<counter<<"\n";
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
    
   	return 0;
}
