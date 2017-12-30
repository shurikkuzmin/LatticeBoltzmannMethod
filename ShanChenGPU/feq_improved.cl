#define BLOCK_SIZE 64
__kernel void CalculateDensity(int nx,int ny,
					__global float* rho,
					__global float* fr0,
					__global float* fe0,
					__global float* fn0,
					__global float* fw0,
					__global float* fs0,
					__global float* fne0,
					__global float* fnw0,
					__global float* fsw0,
					__global float* fse0)
{
	int gx=get_global_id(0); 
	int lx=get_local_id(0);
	int gy=get_global_id(1);

	int k=gx+gy*nx;
	
	float fr = fr0[k];
	float fe = fe0[k];
	float fs = fs0[k];
	float fn = fn0[k];
	float fw = fw0[k];
	float fse = fse0[k];
	float fsw = fsw0[k];
	float fne = fne0[k];
	float fnw = fnw0[k];

	rho[k]=fr+fe+fs+fn+fw+fse+fsw+fne+fnw;

	barrier(CLK_GLOBAL_MEM_FENCE);

}


__kernel void CalculateForce(int nx,int ny,
					__global float* fx,
					__global float* fy,
					__global float* rho)
{
	int gx=get_global_id(0); 
	int lx=get_local_id(0);
	int gy=get_global_id(1);

	int k=gx+gy*nx;

   	float g=-5.0, tau = 1.0f;
	
	__local float shared_rhon[BLOCK_SIZE];
	__local float shared_rhos[BLOCK_SIZE];
	__local float shared_rho[BLOCK_SIZE];
	
	
	int top=(gy+1+ny)%ny;
	int bottom=(gy-1+ny)%ny;
	
	int first=(gx+nx-1)%nx;
	int last=(gx+nx+1)%nx;

	int up=gx+top*nx;
	int down=gx+bottom*nx;

	float rhofirstup=rho[first+top*nx];
	float rhofirst=rho[first+gy*nx];
	float rhofirstdown=rho[first+bottom*nx];
	float rholastup=rho[last+top*nx];
	float rholast=rho[last+gy*nx];
	float rholastdown=rho[last+bottom*nx];
	
	shared_rhon[lx]=rho[up];
	shared_rhos[lx]=rho[down];
	shared_rho[lx]=rho[k];
	barrier(CLK_LOCAL_MEM_FENCE);


	int left=lx-1;
	int right=lx+1;

	if(lx==0)
	{
		fx[k]=1.0/9.0*((1.0-exp(-shared_rho[1]))-(1.0-exp(-rhofirst)))
				+1.0/36.0*((1.0-exp(-shared_rhon[1]))-(1.0-exp(-rhofirstup))-(1.0-exp(-rhofirstdown))+(1.0-exp(-shared_rhos[1])));
		fy[k]=1.0/9.0*((1.0-exp(-shared_rhon[lx]))-(1.0-exp(-shared_rhos[lx])))
				+1.0/36.0*((1.0-exp(-shared_rhon[1]))+(1.0-exp(-rhofirstup))-(1.0-exp(-rhofirstdown))-(1.0-exp(-shared_rhos[1])));
	}
	else if(lx==BLOCK_SIZE-1)
	{
		fx[k]=1.0/9.0*((1.0-exp(-rholast))-(1.0-exp(-shared_rho[BLOCK_SIZE-2])))
				+1.0/36.0*((1.0-exp(-rholastup))-(1.0-exp(-shared_rhon[BLOCK_SIZE-2]))-(1.0-exp(-shared_rhos[BLOCK_SIZE-2]))+(1.0-exp(-rholastdown)));
		fy[k]=1.0/9.0*((1.0-exp(-shared_rhon[lx]))-(1.0-exp(-shared_rhos[lx])))
				+1.0/36.0*((1.0-exp(-rholastup))+(1.0-exp(-shared_rhon[BLOCK_SIZE-2]))-(1.0-exp(-shared_rhos[BLOCK_SIZE-2]))-(1.0-exp(-rholastdown)));
	}
	else
	{
		fx[k]=1.0/9.0*((1.0-exp(-shared_rho[right]))-(1.0-exp(-shared_rho[left])))
				+1.0/36.0*((1.0-exp(-shared_rhon[right]))-(1.0-exp(-shared_rhon[left]))-(1.0-exp(-shared_rhos[left]))+(1.0-exp(-shared_rhos[right])));
		fy[k]=1.0/9.0*((1.0-exp(-shared_rhon[lx]))-(1.0-exp(-shared_rhos[lx])))
				+1.0/36.0*((1.0-exp(-shared_rhon[right]))+(1.0-exp(-shared_rhon[left]))-(1.0-exp(-shared_rhos[left]))-(1.0-exp(-shared_rhos[right])));
	}
	
	fx[k]*=-g*(1.0-exp(-rho[k]));
	fy[k]*=-g*(1.0-exp(-rho[k]));
	
	barrier(CLK_LOCAL_MEM_FENCE);
}

__kernel void CollisionPropogate(int nx,int ny,
					__global float* fx,
					__global float* fy,
					__global float* fr0,
					__global float* fe0,
					__global float* fn0,
					__global float* fw0,
					__global float* fs0,
					__global float* fne0,
					__global float* fnw0,
					__global float* fsw0,
					__global float* fse0,
					__global float* fr1,
					__global float* fe1,
					__global float* fn1,
					__global float* fw1,
					__global float* fs1,
					__global float* fne1,
					__global float* fnw1,
					__global float* fsw1,
					__global float* fse1)
{
    int lx=get_global_id(0); 
	int ix=get_local_id(0);
	int by=get_global_id(1);

	int k=lx+by*nx;

   	float g=-5.0f, tau = 1.0f; 

	float fr = fr0[k];
	float fe = fe0[k];
	float fs = fs0[k];
	float fn = fn0[k];
	float fw = fw0[k];
	float fse = fse0[k];
	float fsw = fsw0[k];
	float fne = fne0[k];
	float fnw = fnw0[k];

	//Probably it's better to calculate the force with sum of populations instead of addressing global memory
	float rhocurr=fr+fe+fn+fw+fs+fne+fnw+fsw+fse; 
	float u1=(fe-fw+fne-fnw-fsw+fse)/rhocurr+fx[k]/(tau*rhocurr); 
	float u2=(fn-fs+fne+fnw-fse-fsw)/rhocurr+fy[k]/(tau*rhocurr); 

	float usq = u1*u1 + u2*u2;
	float feq0 = 4.0/9.0 * rhocurr * (1.0 - 1.5 * usq);
	float feq1 = 1.0/9.0 * rhocurr * (1.0 + 3*u1 + 4.5*u1*u1 - 1.5*usq); 
	float feq2 = 1.0/9.0 * rhocurr * (1.0 + 3*u2 + 4.5*u2*u2 - 1.5*usq); 
	float feq3 = 1.0/9.0 * rhocurr * (1.0 - 3*u1 + 4.5*u1*u1 - 1.5*usq); 
	float feq4 = 1.0/9.0 * rhocurr * (1.0- 3*u2 + 4.5*u2*u2 - 1.5*usq); 
	float feq5 = 1.0/36.0 * rhocurr * (1.0 + 3*(u1 + u2) + 4.5*(u1 + u2)*(u1 + u2) - 1.5*usq); 
	float feq6 = 1.0/36.0 * rhocurr * (1.0 + 3*(-u1 + u2) + 4.5*(-u1 + u2)*(-u1 + u2) - 1.5*usq); 
	float feq7 = 1.0/36.0 * rhocurr * (1.0 + 3*(-u1 - u2) + 4.5*(u1 + u2)*(u1 + u2) - 1.5*usq); 
	float feq8 = 1.0/36.0 * rhocurr * (1.0 + 3*(u1 - u2) + 4.5*(u1 - u2)*(u1 -u2) - 1.5*usq); 

	fr+=-(fr-feq0)/tau;
	fe+=-(fe-feq1)/tau;
	fn+=-(fn-feq2)/tau;
	fw+=-(fw-feq3)/tau;
	fs+=-(fs-feq4)/tau;
	fne+=-(fne-feq5)/tau;
	fnw+=-(fnw-feq6)/tau;
	fsw+=-(fsw-feq7)/tau;
	fse+=-(fse-feq8)/tau;

	__local float shared_fe[BLOCK_SIZE];
	__local float shared_fw[BLOCK_SIZE];
	__local float shared_fne[BLOCK_SIZE];
	__local float shared_fnw[BLOCK_SIZE];
	__local float shared_fsw[BLOCK_SIZE];
	__local float shared_fse[BLOCK_SIZE];
	
	int up=(ix+1+BLOCK_SIZE)%BLOCK_SIZE;
	int down=(ix-1+BLOCK_SIZE)%BLOCK_SIZE;
	shared_fe[up]=fe;
	shared_fne[up]=fne;
	shared_fse[up]=fse;

	shared_fw[down]=fw;
	shared_fnw[down]=fnw;
	shared_fsw[down]=fsw;
	
	
	
	barrier(CLK_LOCAL_MEM_FENCE);
	
	//Exchange in the shared memory location
	fr1[k]=fr;
	fe1[k]=shared_fe[ix];
	fw1[k]=shared_fw[ix];
	
	//Exchange for the upper and bottom limit limit
	up=nx*((by+1+ny)%ny)+lx;
	down=nx*((by-1+ny)%ny)+lx;
	fn1[up]=fn;
	fne1[up]=shared_fne[ix];
	fnw1[up]=shared_fnw[ix];
	fs1[down]=fs;
	fse1[down]=shared_fse[ix];
	fsw1[down]=shared_fsw[ix];
	

	barrier(CLK_LOCAL_MEM_FENCE);
}

__kernel void FinishPropogate(int nx,int ny,
					__global float* fr0,
					__global float* fe0,
					__global float* fn0,
					__global float* fw0,
					__global float* fs0,
					__global float* fne0,
					__global float* fnw0,
					__global float* fsw0,
					__global float* fse0,
					__global float* fr1,
					__global float* fe1,
					__global float* fn1,
					__global float* fw1,
					__global float* fs1,
					__global float* fne1,
					__global float* fnw1,
					__global float* fsw1,
					__global float* fse1)
{
	int bx=nx/BLOCK_SIZE;
	int ystart = get_global_id(1);

	float tempw=fw1[BLOCK_SIZE-1 + nx*ystart];
	float tempnw=fnw1[BLOCK_SIZE-1 + nx*ystart];
	float tempsw=fsw1[BLOCK_SIZE-1 + nx*ystart];
	float tempe=fe1[nx-BLOCK_SIZE + nx*ystart];
	float tempne=fne1[nx-BLOCK_SIZE + nx*ystart];
	float tempse=fse1[nx-BLOCK_SIZE + nx*ystart];

	int xtarget,xstart,ktarget,kstart;

	for(int i=1;i<bx;i++)
	{
		xtarget=i*BLOCK_SIZE-1;
		xstart=(i+1)*BLOCK_SIZE-1;
		ktarget=xtarget+nx*ystart;
		kstart=xstart+nx*ystart;
		
		fw1[ktarget]=fw1[kstart];
		fnw1[ktarget]=fnw1[kstart];
		fsw1[ktarget]=fsw1[kstart];
	}

	for(int i=bx-1;i>0;i--)
	{
		xtarget=i*BLOCK_SIZE;
		xstart=(i-1)*BLOCK_SIZE;
		ktarget=xtarget+nx*ystart;
		kstart=xstart+nx*ystart;
		
		fe1[ktarget]=fe1[kstart];
		fne1[ktarget]=fne1[kstart];
		fse1[ktarget]=fse1[kstart];
	}

	fw1[nx-1 + nx*ystart]=tempw;
	fnw1[nx-1 + nx*ystart]=tempnw;
	fsw1[nx-1 + nx*ystart]=tempsw;
	fe1[0 + nx*ystart]=tempe;
	fne1[0 + nx*ystart]=tempne;
	fse1[0 + nx*ystart]=tempse;

}
