#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>

const int NY=3;
//const int NX=(NY-2)*16;
const int NX=40;
const int N=800;
const double pi=3.141592653589793238462643383279502884197;

double f[NX][NY][9],f2[NX][NY][9];
double rho[NX][NY],ux[NX][NY],uy[NX][NY];


void writedensity(std::string const & fName)
{
	std::string fullName = "./tmp/" + fName+ ".dat";
	std::ofstream fout(fullName.c_str());
	fout.precision(15);

	for (int iY=NY-1; iY>=0; --iY)
	{
		for (int iX=0; iX<NX; ++iX)
			fout<<rho[iX][iY]<<" ";
		fout<<"\n";
	}

}

void writevelocity(std::string const & fName)
{
	std::string fullName = "./tmp/" + fName+ ".dat";
	std::ofstream fout(fullName.c_str());
	fout.precision(15);

	for (int iY=NY-1; iY>=0; --iY)
	{
		for (int iX=0; iX<NX; ++iX)
			fout<<ux[iX][iY]<<" ";
		fout<<"\n";
	}

}



int main(int argc, char* argv[])
{

	double omega=double(atoi(argv[1]))/10.0;
	double omega_bulk=double(atoi(argv[2]))/10.0;

	//std::stringstream additional_name;
	//if (argc>2) 
	//	additional_name<<argv[3];
	std::cout<<"Omega="<<omega<<"\n";
	std::cout<<"Omega_bulk="<<omega_bulk<<"\n";
	double omegaginzburg=8.0*(2.0-omega)/(8.0-omega);
	double omegamat[]={1.0,omega_bulk,omega,1.0,omega,1.0,omega,omega,omega};
	double weights[]={4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};

	double weightsmat[]={1.0/9.0,1.0/36.0,1.0/36.0,1.0/6.0,1.0/12.0,1.0/6.0,1.0/12.0,1.0/4.0,1.0/4.0};
	int cx[]={0,1,0,-1,0,1,-1,-1,1};
	int cy[]={0,0,1,0,-1,1,1,-1,-1};
	double M[9][9];
	double invM[9][9];
	for (int iCoor=0;iCoor<9;iCoor++)
	{
		M[0][iCoor]=1.0;
		M[1][iCoor]=-4.0+3.0*(cx[iCoor]*cx[iCoor]+cy[iCoor]*cy[iCoor]);
		M[2][iCoor]=4.0-21.0/2.0*(cx[iCoor]*cx[iCoor]+cy[iCoor]*cy[iCoor])+9.0/2.0*(cx[iCoor]*cx[iCoor]+cy[iCoor]*cy[iCoor])*(cx[iCoor]*cx[iCoor]+cy[iCoor]*cy[iCoor]);
		M[3][iCoor]=cx[iCoor];
		M[4][iCoor]=(-5.0+3.0*(cx[iCoor]*cx[iCoor]+cy[iCoor]*cy[iCoor]))*cx[iCoor];
		M[5][iCoor]=cy[iCoor];
		M[6][iCoor]=(-5.0+3.0*(cx[iCoor]*cx[iCoor]+cy[iCoor]*cy[iCoor]))*cy[iCoor];
		M[7][iCoor]=cx[iCoor]*cx[iCoor]-cy[iCoor]*cy[iCoor];
		M[8][iCoor]=cx[iCoor]*cy[iCoor];

		invM[iCoor][0]=weightsmat[0]*1.0;
		invM[iCoor][1]=weightsmat[1]*(-4.0+3.0*(cx[iCoor]*cx[iCoor]+cy[iCoor]*cy[iCoor]));
		invM[iCoor][2]=weightsmat[2]*(4.0-21.0/2.0*(cx[iCoor]*cx[iCoor]+cy[iCoor]*cy[iCoor])
			+9.0/2.0*(cx[iCoor]*cx[iCoor]+cy[iCoor]*cy[iCoor])*(cx[iCoor]*cx[iCoor]+cy[iCoor]*cy[iCoor]));
		invM[iCoor][3]=weightsmat[3]*cx[iCoor];
		invM[iCoor][4]=weightsmat[4]*(-5.0+3.0*(cx[iCoor]*cx[iCoor]+cy[iCoor]*cy[iCoor]))*cx[iCoor];
		invM[iCoor][5]=weightsmat[5]*cy[iCoor];
		invM[iCoor][6]=weightsmat[6]*(-5.0+3.0*(cx[iCoor]*cx[iCoor]+cy[iCoor]*cy[iCoor]))*cy[iCoor];
		invM[iCoor][7]=weightsmat[7]*(cx[iCoor]*cx[iCoor]-cy[iCoor]*cy[iCoor]);
		invM[iCoor][8]=weightsmat[8]*cx[iCoor]*cy[iCoor];
	}


	//Initialization of initial conditions

	double ampl=0.001;
	double force_ampl=0.002;
	for(int iX=0;iX<NX;iX++)
		for(int iY=0;iY<NY;iY++)
		{
			double forcex=force_ampl*cos(double(iX)/double(NX)*2.0*pi);
			double forcey=0.0;
			rho[iX][iY]=1.0;
			ux[iX][iY]=ampl*cos(double(iX)/double(NX)*2.0*pi)-forcex/2.0;
			uy[iX][iY]=0.0;

			double fluxx=3.0*ux[iX][iY];
			double fluxy=3.0*uy[iX][iY];

			double qxx=4.5*ux[iX][iY]*ux[iX][iY];
			double qxy=9.0*ux[iX][iY]*uy[iX][iY];
			double qyy=4.5*uy[iX][iY]*uy[iX][iY];

			for(int iPop=0;iPop<9;iPop++)
			{
				f[iX][iY][iPop]=rho[iX][iY]*weights[iPop]*(1.0+fluxx*cx[iPop]+fluxy*cy[iPop]+
								qxx*(cx[iPop]*cx[iPop]-1.0/3.0)+qxy*cx[iPop]*cy[iPop]+qyy*(cy[iPop]*cy[iPop]-1.0/3.0));
			}
		}


    //MRT initialization
//    for(int iX=0;iX<NX;iX++)
//		for(int iY=0;iY<NY;iY++)
//		{
//			double forcex=force_ampl*cos(double(iX)/double(NX)*2.0*pi);
//			double forcey=0.0;
//			rho[iX][iY]=1.0;
//			ux[iX][iY]=ampl*cos(double(iX)/double(NX)*2.0*pi)-forcex/2.0;
//			uy[iX][iY]=0.0;
//
//			for(int iPop=0;iPop<9;iPop++)
//			{
//				f[iX][iY][iPop]=rho[iX][iY]*weights[iPop]*(1.0+fluxx*cx[iPop]+fluxy*cy[iPop]+
//								qxx*(cx[iPop]*cx[iPop]-1.0/3.0)+qxy*cx[iPop]*cy[iPop]+qyy*(cy[iPop]*cy[iPop]-1.0/3.0));
//			}
//		}





	for(int counter=0;counter<=N;counter++)
	{
		for(int iX=0;iX<NX;iX++)
			for(int iY=0;iY<NY;iY++)
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

				double dense=rho[iX][iY];

				ux[iX][iY]/=dense;
				uy[iX][iY]/=dense;

 				double forcex=force_ampl*cos(2.0*pi*double(iX)/double(NX));
 				double forcey=0;

 				ux[iX][iY]+=forcex/(2.0*dense);
 				ux[iX][iY]+=forcey/(2.0*dense);

				double uxeq=ux[iX][iY];
				double uyeq=uy[iX][iY];

				//Construction of the equilibrium moments
				double eq[9];
				eq[0]=dense;
				eq[1]=-2.0*dense+3.0*dense*(uxeq*uxeq+uyeq*uyeq);
				eq[2]=dense-3.0*dense*(uxeq*uxeq+uyeq*uyeq);
				eq[3]=dense*uxeq;
				eq[4]=-dense*uxeq;
				eq[5]=dense*uyeq;
				eq[6]=-dense*uyeq;
				eq[7]=dense*(uxeq*uxeq-uyeq*uyeq);
				eq[8]=dense*uxeq*uyeq;

				double momforce[9];
				momforce[0]=0.0;
				momforce[1]=(1.0-omegamat[1]/2.0)*6.0*(uxeq*forcex+uyeq*forcey);
				momforce[2]=0.0;
				momforce[3]=(1.0-omegamat[3]/2.0)*forcex;
				momforce[4]=0.0;
				momforce[5]=(1.0-omegamat[5]/2.0)*forcey;
				momforce[6]=0.0;
				momforce[7]=(1.0-omegamat[7]/2.0)*2.0*(uxeq*forcex-uyeq*forcey);
				momforce[8]=(1.0-omegamat[8]/2.0)*(uxeq*forcey+uyeq*forcex);

				double feqeq[9];
				for(int iPop=0;iPop<9; iPop++)
					feqeq[iPop]=weights[iPop]*dense*(1.0+3.0*cx[iPop]*uxeq+3.0*cy[iPop]*uyeq+
					4.5*(cx[iPop]*cx[iPop]-1.0/3.0)*uxeq*uxeq+9.0*cx[iPop]*cy[iPop]*uxeq*uyeq+4.5*(cy[iPop]*cy[iPop]-1.0/3.0)*uyeq*uyeq);

				double add[9];
				double addit;
				for(int iPop=0;iPop < 9; iPop++)
				{
					add[iPop]=0.0;
					for(int k=0; k < 9; k++)
						add[iPop]=add[iPop]+M[iPop][k]*(-f[iX][iY][k]);
					add[iPop]+=eq[iPop];
					add[iPop]*=omegamat[iPop];
					add[iPop]+=momforce[iPop];
				}

				for(int k=0; k < 9; k++)
				{
//   					feqforce=(1.0-0.5*omega)*weights[k]*(forcex*(3.0*(cx[k]-uxeq)+9.0*cx[k]*(cx[k]*uxeq+cy[k]*uyeq))+
// 							forcey*(3.0*(cy[k]-uyeq)+9.0*cy[k]*(cx[k]*uxeq+cy[k]*uyeq)));
// //
/*						feqforce=(1.0-0.5*omega)*weights[k]*(forcex*(3.0*(cx[k]-uxeq)+0.0*cx[k]*(cx[k]*uxeq+cy[k]*uyeq))+
							forcey*(3.0*(cy[k]-uyeq)+0.0*cy[k]*(cx[k]*uxeq+cy[k]*uyeq)));*/



					//feqforce=3.0*weights[k]*cx[k]*forcex;

 					//feqforce=weights[k]*(3.0*(1.0-0.5*omegamat[1])*forcex*cx[k]+3.0*(1.0-0.5*omegamat[2])*forcey*cy[k]+
 			//				9.0*((1.0-0.5*omegamat[3])*(cx[k]*cx[k]-1.0/3.0)*forcex*uxeq+
// 							(1.0-0.5*omegamat[4])*cx[k]*cy[k]*(forcex*uyeq+forcey*uxeq)+
// 							(1.0-0.5*omegamat[5])*(cy[k]*cy[k]-1.0/3.0)*forcey*uyeq));
//

					addit=0.0;
					for(int m=0; m < 9; m++)
						addit=addit+invM[k][m]*add[m];
					f2[iX][iY][k]=f[iX][iY][k]+addit;
					//f2[iX][iY][k]=f[iX][iY][k]*(1.0-omega)+omega*feqeq[k];
				}


			}


		//Streaming
		for(int iX=0;iX<NX;iX++)
			for(int iY=0;iY<NY;iY++)
				for(int iPop=0;iPop<9;iPop++)
				{
					int iX2=(iX-cx[iPop]+NX)%NX;
					int iY2=(iY-cy[iPop]+NY)%NY;
					f[iX][iY][iPop]=f2[iX2][iY2][iPop];
				}

		//Writing files
		if (counter%1==0)
		{
			std::cout<<counter<<"\n";

			std::stringstream filewritedensity;
 			std::stringstream filewritevelocity;
 			std::stringstream counterconvert;
 			counterconvert<<counter;
 			filewritedensity<<std::fixed;
			filewritevelocity<<std::fixed;

			filewritedensity<<"proba"<<std::string(6-counterconvert.str().size(),'0')<<counter;
			filewritevelocity<<"phase"<<std::string(6-counterconvert.str().size(),'0')<<counter;

 			writedensity(filewritedensity.str());
			writevelocity(filewritevelocity.str());
		}


	}

	return 0;
}
