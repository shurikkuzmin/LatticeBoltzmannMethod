/*
 *  main.cpp
 *  Shan-Chen GPU code
 *
 *  Created by Alexandr Kuzmin on 09-09-24.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

////////////////////////////////////////////////////////////////////////////////

#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <CL/cl.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

#define BLOCK_SIZE 64
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
static char *
load_program_source(const char *filename)
{
    struct stat statbuf;
    FILE        *fh;
    char        *source;
	
    fh=fopen(filename, "r");
	if (fh == 0)
	{
		return 0;
	}
    stat(filename, &statbuf);
    source = (char *) malloc(statbuf.st_size + 1);
    fread(source, statbuf.st_size, 1, fh);
    source[statbuf.st_size] = '\0';
	
    return source;
}

void writeFile(std::string name, float* f, int nx,int ny)
{
	std::ofstream fout(name.c_str());
	for(int counterX=0;counterX<nx;counterX++)
	{
		for(int counterY=0;counterY<ny;counterY++)
			fout<<f[counterX*ny+counterY]<<" ";
		fout<<"\n";
	}
	fout<<std::endl;
}

int main(int argc, char** argv)
{
    int err;                            // error code returned from api calls
	int nx = 256;
	int ny=256;
	double rhol=1.8;
	double rhog=0.2;
	double g=-5.0;
	double tau=1.0;
	int radius=40;
	
	float host_fr0[nx*ny];
	float host_fe0[nx*ny];
	float host_fn0[nx*ny];
	float host_fw0[nx*ny];
	float host_fs0[nx*ny];
	float host_fne0[nx*ny];
	float host_fnw0[nx*ny];
	float host_fsw0[nx*ny];
	float host_fse0[nx*ny];
	float host_fr1[nx*ny];
	float host_fe1[nx*ny];
	float host_fn1[nx*ny];
	float host_fw1[nx*ny];
	float host_fs1[nx*ny];
	float host_fne1[nx*ny];
	float host_fnw1[nx*ny];
	float host_fsw1[nx*ny];
	float host_fse1[nx*ny];

	float host_rho[nx*ny];
	float host_fx[nx*ny];
	float host_fy[nx*ny];
	
	
	const unsigned int iNumElements = nx*ny;
	int i;
	for (i=0; i<nx*ny; i++)
	{		
		if ((i/nx-double(ny)/2.0)*(i/nx-double(ny)/2.0)+(i%nx-double(nx)/2.0)*(i%nx-double(nx)/2.0)<=radius*radius)
		{
			host_rho[i]=rhol;
		}
		else 
			 host_rho[i]=rhog;
		double rho=host_rho[i];
		double u1=0;
		double u2=0;
		
		float usq = u1*u1 + u2*u2;
		host_fr1[i] = host_fr0[i] =  4.0f/9.0f * rho * (1.0 - 1.5f * usq);
		host_fe1[i] = host_fe0[i] = 1.0f/9.0f * rho * (1.0 + 3*u1 + 4.5f*u1*u1 - 1.5f*usq); 
		host_fn1[i] = host_fn0[i] = 1.0f/9.0f * rho * (1.0 + 3*u2 + 4.5f*u2*u2 - 1.5f*usq); 
		host_fw1[i] = host_fw0[i] = 1.0f/9.0f * rho * (1.0- 3*u1 + 4.5f*u1*u1 - 1.5f*usq); 
		host_fs1[i] = host_fs0[i] = 1.0f/9.0f * rho * (1.0 - 3*u2 + 4.5f*u2*u2 - 1.5f*usq); 
		host_fne1[i] = host_fne0[i] = 1.0f/36.0f * rho * (1.0 + 3*(u1 + u2) + 4.5f*(u1 + u2)*(u1 + u2) - 1.5f*usq); 
		host_fnw1[i] = host_fnw0[i] = 1.0f/36.0f * rho * (1.0+ 3*(-u1 + u2) + 4.5f*(-u1 + u2)*(-u1 + u2) - 1.5f*usq); 
		host_fsw1[i] = host_fsw0[i] = 1.0f/36.0f * rho * (1.0 + 3*(-u1 - u2) + 4.5f*(u1 + u2)*(u1 + u2) - 1.5f*usq); 
		host_fse1[i] = host_fse0[i] = 1.0f/36.0f * rho * (1.0 + 3*(u1 - u2) + 4.5f*(u1 - u2)*(u1 -u2) - 1.5f*usq); 
	}


    size_t global[2];						// global domain size for our calculation
    size_t local[2];						// local domain size for our calculation

    cl_platform_id platform;
    cl_device_id device_id;					// compute device id
    cl_context context;						// compute context
    cl_command_queue command_queue;         // compute command queue
    cl_program program;						// compute program
    cl_kernel kernelCalculateDensity;
    cl_kernel kernelCollProp;               // compute kernel
	cl_kernel kernelFinishProp;				// compute kernel
	cl_kernel kernelCalculateForce;			// compute kernel

			 
    int gpu = 1;
	clGetPlatformIDs(1, &platform, NULL);

	std::cout << "Connecting to a compute device" << std::endl;
    err = clGetDeviceIDs(platform, gpu ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, 1, &device_id, NULL);
    if (err != CL_SUCCESS)
    {
		std::cout<<"Error: Failed to create a device group!\n";
        return EXIT_FAILURE;
    }

    // Create a compute context
    std::cout << "Creating a compute context" << std::endl;
    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
    if (!context)
    {
        printf("Error: Failed to create a compute context!\n");
        return EXIT_FAILURE;
    }

    // Create a command command_queue
    std::cout << "Creating a command queue" << std::endl;
    command_queue = clCreateCommandQueue(context, device_id, 0, &err);
    if (!command_queue)
    {
        printf("Error: Failed to create a command command_queue!\n");
        return EXIT_FAILURE;
    }

    // Create the compute program from the source buffer
    char *source = load_program_source("feq_improved.cl");
	if(!source)
    {
        printf("Error: Failed to load compute program from file!\n");
        return EXIT_FAILURE;    
    }
	
    // Create the compute program from the source buffer
    program = clCreateProgramWithSource(context, 1, (const char **) & source, NULL, &err);
	if (!program)
    {
        printf("Error: Failed to create compute program!\n");
        return EXIT_FAILURE;
    }

    // Build the program executable
    std::cout << "Building the program executable" << std::endl;
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[8192*8];
        printf("Error: Failed to build program executable!\n");
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("%s\n", buffer);
        exit(1);
    }

    // Create the compute kernel in the program we wish to run
	std::cout << "Creating kernel" << std::endl;
	
	kernelCalculateDensity = clCreateKernel(program, "CalculateDensity", &err);
    if (!kernelCalculateDensity || err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel!\n");
        exit(1);
    }

	
	kernelCollProp = clCreateKernel(program, "CollisionPropogate", &err);
    if (!kernelCollProp || err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel!\n");
        exit(1);
    }
	
	kernelFinishProp = clCreateKernel(program, "FinishPropogate", &err);
    if (!kernelCollProp || err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel!\n");
        exit(1);
    }
			 
	kernelCalculateForce = clCreateKernel(program, "CalculateForce", &err);
	if (!kernelCollProp || err != CL_SUCCESS)
	{
		printf("Error: Failed to create compute kernel!\n");
		exit(1);
	}		 
	
	

    // Create the input and output arrays in device memory for our calculation
	std::cout << "Creating buffers" << std::endl;
	cl_mem dev_rho=clCreateBuffer(context,
								CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_rho,NULL);
	 cl_mem dev_fx=clCreateBuffer(context,
							   CL_MEM_READ_WRITE, sizeof(cl_float)*nx*ny, NULL,NULL);
	 cl_mem dev_fy=clCreateBuffer(context,
							   CL_MEM_READ_WRITE, sizeof(cl_float)*nx*ny, NULL,NULL);

	cl_mem dev_fr0 = clCreateBuffer(context,
								  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fr0,NULL);
    cl_mem dev_fe0 = clCreateBuffer(context,
								   CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fe0,NULL);
    cl_mem dev_fn0 = clCreateBuffer(context,
								   CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fn0, NULL);
    cl_mem dev_fw0 = clCreateBuffer(context,
									CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fw0, NULL);
	cl_mem dev_fs0 = clCreateBuffer(context,
								  CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fs0, NULL);
	cl_mem dev_fne0 = clCreateBuffer(context,
									CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fne0, NULL);
    cl_mem dev_fnw0 = clCreateBuffer(context,
									CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fnw0, NULL);
	cl_mem dev_fsw0 = clCreateBuffer(context,
									CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fsw0, NULL);
	cl_mem dev_fse0 = clCreateBuffer(context,
									 CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fse0, NULL);
	cl_mem dev_fr1 = clCreateBuffer(context,
									CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fr1,NULL);
    cl_mem dev_fe1 = clCreateBuffer(context,
									CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fe1,NULL);
    cl_mem dev_fn1 = clCreateBuffer(context,
									CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fn1, NULL);
    cl_mem dev_fw1 = clCreateBuffer(context,
									CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fw1, NULL);
	cl_mem dev_fs1 = clCreateBuffer(context,
									CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fs1, NULL);
	cl_mem dev_fne1 = clCreateBuffer(context,
									 CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fne1, NULL);
    cl_mem dev_fnw1 = clCreateBuffer(context,
									 CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fnw1, NULL);
	cl_mem dev_fsw1 = clCreateBuffer(context,
									 CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fsw1, NULL);
	cl_mem dev_fse1 = clCreateBuffer(context,
									 CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*nx*ny, host_fse1, NULL);
	

	// Set the arguments to our compute kernel
	
	
	
	std::cout << "Setting kernel arguments" << std::endl;
	
	err = 0;
	err |= clSetKernelArg(kernelCalculateDensity, 0, sizeof(int), &nx);
	err |= clSetKernelArg(kernelCalculateDensity, 1, sizeof(int), &ny);
	err |= clSetKernelArg(kernelCalculateDensity, 2, sizeof(cl_mem), &dev_rho);
	err |= clSetKernelArg(kernelCalculateDensity, 3, sizeof(cl_mem), &dev_fr0);
	err |= clSetKernelArg(kernelCalculateDensity, 4, sizeof(cl_mem), &dev_fe0);
	err |= clSetKernelArg(kernelCalculateDensity, 5, sizeof(cl_mem), &dev_fn0);
	err |= clSetKernelArg(kernelCalculateDensity, 6, sizeof(cl_mem), &dev_fw0);
	err |= clSetKernelArg(kernelCalculateDensity, 7, sizeof(cl_mem), &dev_fs0);
	err |= clSetKernelArg(kernelCalculateDensity, 8, sizeof(cl_mem), &dev_fne0);
	err |= clSetKernelArg(kernelCalculateDensity, 9, sizeof(cl_mem), &dev_fnw0);
	err |= clSetKernelArg(kernelCalculateDensity, 10, sizeof(cl_mem), &dev_fsw0);
	err |= clSetKernelArg(kernelCalculateDensity, 11, sizeof(cl_mem), &dev_fse0);			 
	
	err |= clSetKernelArg(kernelCalculateForce, 0, sizeof(int), &nx);
	err |= clSetKernelArg(kernelCalculateForce, 1, sizeof(int), &ny);
	err |= clSetKernelArg(kernelCalculateForce, 2, sizeof(cl_mem), &dev_fx);
	err |= clSetKernelArg(kernelCalculateForce, 3, sizeof(cl_mem), &dev_fy);
	err |= clSetKernelArg(kernelCalculateForce, 4, sizeof(cl_mem), &dev_rho);
			
	err |= clSetKernelArg(kernelCollProp, 0, sizeof(int), &nx);
    err |= clSetKernelArg(kernelCollProp, 1, sizeof(int), &ny);
	err |= clSetKernelArg(kernelCollProp, 2, sizeof(cl_mem), &dev_fx);
	err |= clSetKernelArg(kernelCollProp, 3, sizeof(cl_mem), &dev_fy);
	err |= clSetKernelArg(kernelCollProp, 4, sizeof(cl_mem), &dev_fr0);
	err |= clSetKernelArg(kernelCollProp, 5, sizeof(cl_mem), &dev_fe0);
	err |= clSetKernelArg(kernelCollProp, 6, sizeof(cl_mem), &dev_fn0);
    err |= clSetKernelArg(kernelCollProp, 7, sizeof(cl_mem), &dev_fw0);
    err |= clSetKernelArg(kernelCollProp, 8, sizeof(cl_mem), &dev_fs0);
  	err |= clSetKernelArg(kernelCollProp, 9, sizeof(cl_mem), &dev_fne0);
	err |= clSetKernelArg(kernelCollProp, 10, sizeof(cl_mem), &dev_fnw0);
	err |= clSetKernelArg(kernelCollProp, 11, sizeof(cl_mem), &dev_fsw0);
	err |= clSetKernelArg(kernelCollProp, 12, sizeof(cl_mem), &dev_fse0);
  	err |= clSetKernelArg(kernelCollProp, 13, sizeof(cl_mem), &dev_fr1);
	err |= clSetKernelArg(kernelCollProp, 14, sizeof(cl_mem), &dev_fe1);
	err |= clSetKernelArg(kernelCollProp, 15, sizeof(cl_mem), &dev_fn1);
    err |= clSetKernelArg(kernelCollProp, 16, sizeof(cl_mem), &dev_fw1);
    err |= clSetKernelArg(kernelCollProp, 17, sizeof(cl_mem), &dev_fs1);
  	err |= clSetKernelArg(kernelCollProp, 18, sizeof(cl_mem), &dev_fne1);
	err |= clSetKernelArg(kernelCollProp, 19, sizeof(cl_mem), &dev_fnw1);
	err |= clSetKernelArg(kernelCollProp, 20, sizeof(cl_mem), &dev_fsw1);
	err |= clSetKernelArg(kernelCollProp, 21, sizeof(cl_mem), &dev_fse1);
	
	err |= clSetKernelArg(kernelFinishProp, 0, sizeof(int), &nx);
    err |= clSetKernelArg(kernelFinishProp, 1, sizeof(int), &ny);
	err |= clSetKernelArg(kernelFinishProp, 2, sizeof(cl_mem), &dev_fr0);
	err |= clSetKernelArg(kernelFinishProp, 3, sizeof(cl_mem), &dev_fe0);
	err |= clSetKernelArg(kernelFinishProp, 4, sizeof(cl_mem), &dev_fn0);
    err |= clSetKernelArg(kernelFinishProp, 5, sizeof(cl_mem), &dev_fw0);
    err |= clSetKernelArg(kernelFinishProp, 6, sizeof(cl_mem), &dev_fs0);
  	err |= clSetKernelArg(kernelFinishProp, 7, sizeof(cl_mem), &dev_fne0);
	err |= clSetKernelArg(kernelFinishProp, 8, sizeof(cl_mem), &dev_fnw0);
	err |= clSetKernelArg(kernelFinishProp, 9, sizeof(cl_mem), &dev_fsw0);
	err |= clSetKernelArg(kernelFinishProp, 10, sizeof(cl_mem), &dev_fse0);
  	err |= clSetKernelArg(kernelFinishProp, 11, sizeof(cl_mem), &dev_fr1);
	err |= clSetKernelArg(kernelFinishProp, 12, sizeof(cl_mem), &dev_fe1);
	err |= clSetKernelArg(kernelFinishProp, 13, sizeof(cl_mem), &dev_fn1);
    err |= clSetKernelArg(kernelFinishProp, 14, sizeof(cl_mem), &dev_fw1);
    err |= clSetKernelArg(kernelFinishProp, 15, sizeof(cl_mem), &dev_fs1);
  	err |= clSetKernelArg(kernelFinishProp, 16, sizeof(cl_mem), &dev_fne1);
	err |= clSetKernelArg(kernelFinishProp, 17, sizeof(cl_mem), &dev_fnw1);
	err |= clSetKernelArg(kernelFinishProp, 18, sizeof(cl_mem), &dev_fsw1);
	err |= clSetKernelArg(kernelFinishProp, 19, sizeof(cl_mem), &dev_fse1);
	
	if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set kernel arguments! %d\n", err);
        exit(1);
    }

	time_t start,finish;
	start = time(NULL);
	for (int timecounter=0; timecounter<3000;timecounter++)
	{
		global[0]=nx;
		global[1]=ny;
		local[0]=BLOCK_SIZE;
		local[1]=1;
		err = clEnqueueNDRangeKernel(command_queue,kernelCalculateDensity,2,NULL,global,local,0,NULL,NULL);
		clFinish(command_queue);
		
		err = clEnqueueNDRangeKernel(command_queue,kernelCalculateForce,2,NULL,global,local,0,NULL,NULL);
		clFinish(command_queue);
	
/*		err = clEnqueueReadBuffer( command_queue, dev_fx, CL_TRUE, 0, sizeof(float) * nx*ny, host_fx, 0, NULL, NULL );
		clFinish(command_queue);
		
		std::stringstream imagestream;
		std::stringstream len;
		len<<timecounter;
		imagestream << "force"<<std::string(5-len.str().size(),'0')<<timecounter<<".dat";
		writeFile(imagestream.str(), host_fx, nx,ny);*/
		
		
		
		err = clEnqueueNDRangeKernel(command_queue, kernelCollProp, 2, NULL, global, local, 0, NULL, NULL);
		clFinish(command_queue);

		global[0]=1; 
		global[1]=ny;
		local[0]=1;
		local[1]=BLOCK_SIZE;
		
		err|=clEnqueueNDRangeKernel(command_queue, kernelFinishProp, 2, NULL, global, local, 0, NULL, NULL);
		clFinish(command_queue);
		
		if (err)
		{
			printf("Error: Failed to execute kernel!\n");
			return EXIT_FAILURE;
		}

		if (timecounter & 1) 
		{
	
			clSetKernelArg(kernelCalculateDensity, 3, sizeof(cl_mem), &dev_fr0);
			clSetKernelArg(kernelCalculateDensity, 4, sizeof(cl_mem), &dev_fe0);
			clSetKernelArg(kernelCalculateDensity, 5, sizeof(cl_mem), &dev_fn0);
			clSetKernelArg(kernelCalculateDensity, 6, sizeof(cl_mem), &dev_fw0);
			clSetKernelArg(kernelCalculateDensity, 7, sizeof(cl_mem), &dev_fs0);
			clSetKernelArg(kernelCalculateDensity, 8, sizeof(cl_mem), &dev_fne0);
			clSetKernelArg(kernelCalculateDensity, 9, sizeof(cl_mem), &dev_fnw0);
			clSetKernelArg(kernelCalculateDensity, 10, sizeof(cl_mem), &dev_fsw0);
			clSetKernelArg(kernelCalculateDensity, 11, sizeof(cl_mem), &dev_fse0);

			
			clSetKernelArg(kernelCollProp, 4, sizeof(cl_mem), &dev_fr0);
			clSetKernelArg(kernelCollProp, 5, sizeof(cl_mem), &dev_fe0);
			clSetKernelArg(kernelCollProp, 6, sizeof(cl_mem), &dev_fn0);
			clSetKernelArg(kernelCollProp, 7, sizeof(cl_mem), &dev_fw0);
			clSetKernelArg(kernelCollProp, 8, sizeof(cl_mem), &dev_fs0);
			clSetKernelArg(kernelCollProp, 9, sizeof(cl_mem), &dev_fne0);
			clSetKernelArg(kernelCollProp, 10, sizeof(cl_mem), &dev_fnw0);
			clSetKernelArg(kernelCollProp, 11, sizeof(cl_mem), &dev_fsw0);
			clSetKernelArg(kernelCollProp, 12, sizeof(cl_mem), &dev_fse0);
			clSetKernelArg(kernelCollProp, 13, sizeof(cl_mem), &dev_fr1);
			clSetKernelArg(kernelCollProp, 14, sizeof(cl_mem), &dev_fe1);
			clSetKernelArg(kernelCollProp, 15, sizeof(cl_mem), &dev_fn1);
			clSetKernelArg(kernelCollProp, 16, sizeof(cl_mem), &dev_fw1);
			clSetKernelArg(kernelCollProp, 17, sizeof(cl_mem), &dev_fs1);
			clSetKernelArg(kernelCollProp, 18, sizeof(cl_mem), &dev_fne1);
			clSetKernelArg(kernelCollProp, 19, sizeof(cl_mem), &dev_fnw1);
			clSetKernelArg(kernelCollProp, 20, sizeof(cl_mem), &dev_fsw1);
			clSetKernelArg(kernelCollProp, 21, sizeof(cl_mem), &dev_fse1);
		
			clSetKernelArg(kernelFinishProp, 2, sizeof(cl_mem), &dev_fr0);
			clSetKernelArg(kernelFinishProp, 3, sizeof(cl_mem), &dev_fe0);
			clSetKernelArg(kernelFinishProp, 4, sizeof(cl_mem), &dev_fn0);
			clSetKernelArg(kernelFinishProp, 5, sizeof(cl_mem), &dev_fw0);
			clSetKernelArg(kernelFinishProp, 6, sizeof(cl_mem), &dev_fs0);
			clSetKernelArg(kernelFinishProp, 7, sizeof(cl_mem), &dev_fne0);
			clSetKernelArg(kernelFinishProp, 8, sizeof(cl_mem), &dev_fnw0);
			clSetKernelArg(kernelFinishProp, 9, sizeof(cl_mem), &dev_fsw0);
			clSetKernelArg(kernelFinishProp, 10, sizeof(cl_mem), &dev_fse0);
			clSetKernelArg(kernelFinishProp, 11, sizeof(cl_mem), &dev_fr1);
			clSetKernelArg(kernelFinishProp, 12, sizeof(cl_mem), &dev_fe1);
			clSetKernelArg(kernelFinishProp, 13, sizeof(cl_mem), &dev_fn1);
			clSetKernelArg(kernelFinishProp, 14, sizeof(cl_mem), &dev_fw1);
			clSetKernelArg(kernelFinishProp, 15, sizeof(cl_mem), &dev_fs1);
			clSetKernelArg(kernelFinishProp, 16, sizeof(cl_mem), &dev_fne1);
			clSetKernelArg(kernelFinishProp, 17, sizeof(cl_mem), &dev_fnw1);
			clSetKernelArg(kernelFinishProp, 18, sizeof(cl_mem), &dev_fsw1);
			clSetKernelArg(kernelFinishProp, 19, sizeof(cl_mem), &dev_fse1);
		} 
		else 
		{

			clSetKernelArg(kernelCalculateDensity, 3, sizeof(cl_mem), &dev_fr1);
			clSetKernelArg(kernelCalculateDensity, 4, sizeof(cl_mem), &dev_fe1);
			clSetKernelArg(kernelCalculateDensity, 5, sizeof(cl_mem), &dev_fn1);
			clSetKernelArg(kernelCalculateDensity, 6, sizeof(cl_mem), &dev_fw1);
			clSetKernelArg(kernelCalculateDensity, 7, sizeof(cl_mem), &dev_fs1);
			clSetKernelArg(kernelCalculateDensity, 8, sizeof(cl_mem), &dev_fne1);
			clSetKernelArg(kernelCalculateDensity, 9, sizeof(cl_mem), &dev_fnw1);
			clSetKernelArg(kernelCalculateDensity, 10, sizeof(cl_mem), &dev_fsw1);
			clSetKernelArg(kernelCalculateDensity, 11, sizeof(cl_mem), &dev_fse1);
			
			
			clSetKernelArg(kernelCollProp, 4, sizeof(cl_mem), &dev_fr1);
			clSetKernelArg(kernelCollProp, 5, sizeof(cl_mem), &dev_fe1);
			clSetKernelArg(kernelCollProp, 6, sizeof(cl_mem), &dev_fn1);
			clSetKernelArg(kernelCollProp, 7, sizeof(cl_mem), &dev_fw1);
			clSetKernelArg(kernelCollProp, 8, sizeof(cl_mem), &dev_fs1);
			clSetKernelArg(kernelCollProp, 9, sizeof(cl_mem), &dev_fne1);
			clSetKernelArg(kernelCollProp, 10, sizeof(cl_mem), &dev_fnw1);
			clSetKernelArg(kernelCollProp, 11, sizeof(cl_mem), &dev_fsw1);
			clSetKernelArg(kernelCollProp, 12, sizeof(cl_mem), &dev_fse1);
			clSetKernelArg(kernelCollProp, 13, sizeof(cl_mem), &dev_fr0);
			clSetKernelArg(kernelCollProp, 14, sizeof(cl_mem), &dev_fe0);
			clSetKernelArg(kernelCollProp, 15, sizeof(cl_mem), &dev_fn0);
			clSetKernelArg(kernelCollProp, 16, sizeof(cl_mem), &dev_fw0);
			clSetKernelArg(kernelCollProp, 17, sizeof(cl_mem), &dev_fs0);
			clSetKernelArg(kernelCollProp, 18, sizeof(cl_mem), &dev_fne0);
			clSetKernelArg(kernelCollProp, 19, sizeof(cl_mem), &dev_fnw0);
			clSetKernelArg(kernelCollProp, 20, sizeof(cl_mem), &dev_fsw0);
			clSetKernelArg(kernelCollProp, 21, sizeof(cl_mem), &dev_fse0);
		
			clSetKernelArg(kernelFinishProp, 2, sizeof(cl_mem), &dev_fr1);
			clSetKernelArg(kernelFinishProp, 3, sizeof(cl_mem), &dev_fe1);
			clSetKernelArg(kernelFinishProp, 4, sizeof(cl_mem), &dev_fn1);
			clSetKernelArg(kernelFinishProp, 5, sizeof(cl_mem), &dev_fw1);
			clSetKernelArg(kernelFinishProp, 6, sizeof(cl_mem), &dev_fs1);
			clSetKernelArg(kernelFinishProp, 7, sizeof(cl_mem), &dev_fne1);
			clSetKernelArg(kernelFinishProp, 8, sizeof(cl_mem), &dev_fnw1);
			clSetKernelArg(kernelFinishProp, 9, sizeof(cl_mem), &dev_fsw1);
			clSetKernelArg(kernelFinishProp, 10, sizeof(cl_mem), &dev_fse1);
			clSetKernelArg(kernelFinishProp, 11, sizeof(cl_mem), &dev_fr0);
			clSetKernelArg(kernelFinishProp, 12, sizeof(cl_mem), &dev_fe0);
			clSetKernelArg(kernelFinishProp, 13, sizeof(cl_mem), &dev_fn0);
			clSetKernelArg(kernelFinishProp, 14, sizeof(cl_mem), &dev_fw0);
			clSetKernelArg(kernelFinishProp, 15, sizeof(cl_mem), &dev_fs0);
			clSetKernelArg(kernelFinishProp, 16, sizeof(cl_mem), &dev_fne0);
			clSetKernelArg(kernelFinishProp, 17, sizeof(cl_mem), &dev_fnw0);
			clSetKernelArg(kernelFinishProp, 18, sizeof(cl_mem), &dev_fsw0);
			clSetKernelArg(kernelFinishProp, 19, sizeof(cl_mem), &dev_fse0);
		}
		
		
		if (timecounter%10==0) 
		{
			if((timecounter & 1) == 0)
			{
				err = clEnqueueReadBuffer( command_queue, dev_fr0, CL_TRUE, 0, sizeof(float) * nx*ny, host_fr1, 0, NULL, NULL );
				err|= clEnqueueReadBuffer( command_queue, dev_fe0, CL_TRUE, 0, sizeof(float) * nx*ny, host_fe1, 0, NULL, NULL );
				err|= clEnqueueReadBuffer( command_queue, dev_fn0, CL_TRUE, 0, sizeof(float) * nx*ny, host_fn1, 0, NULL, NULL );
				err|= clEnqueueReadBuffer( command_queue, dev_fw0, CL_TRUE, 0, sizeof(float) * nx*ny, host_fw1, 0, NULL, NULL );
				err|= clEnqueueReadBuffer( command_queue, dev_fs0, CL_TRUE, 0, sizeof(float) * nx*ny, host_fs1, 0, NULL, NULL );
				err|= clEnqueueReadBuffer( command_queue, dev_fne0, CL_TRUE, 0, sizeof(float) * nx*ny, host_fne1, 0, NULL, NULL );
				err|= clEnqueueReadBuffer( command_queue, dev_fnw0, CL_TRUE, 0, sizeof(float) * nx*ny, host_fnw1, 0, NULL, NULL );
				err|= clEnqueueReadBuffer( command_queue, dev_fsw0, CL_TRUE, 0, sizeof(float) * nx*ny, host_fsw1, 0, NULL, NULL );
				err|= clEnqueueReadBuffer( command_queue, dev_fse0, CL_TRUE, 0, sizeof(float) * nx*ny, host_fse1, 0, NULL, NULL );
				clFinish(command_queue);
			} else {
				err = clEnqueueReadBuffer( command_queue, dev_fr1, CL_TRUE, 0, sizeof(float) * nx*ny, host_fr1, 0, NULL, NULL );
				err|= clEnqueueReadBuffer( command_queue, dev_fe1, CL_TRUE, 0, sizeof(float) * nx*ny, host_fe1, 0, NULL, NULL );
				err|= clEnqueueReadBuffer( command_queue, dev_fn1, CL_TRUE, 0, sizeof(float) * nx*ny, host_fn1, 0, NULL, NULL );
				err|= clEnqueueReadBuffer( command_queue, dev_fw1, CL_TRUE, 0, sizeof(float) * nx*ny, host_fw1, 0, NULL, NULL );
				err|= clEnqueueReadBuffer( command_queue, dev_fs1, CL_TRUE, 0, sizeof(float) * nx*ny, host_fs1, 0, NULL, NULL );
				err|= clEnqueueReadBuffer( command_queue, dev_fne1, CL_TRUE, 0, sizeof(float) * nx*ny, host_fne1, 0, NULL, NULL );
				err|= clEnqueueReadBuffer( command_queue, dev_fnw1, CL_TRUE, 0, sizeof(float) * nx*ny, host_fnw1, 0, NULL, NULL );
				err|= clEnqueueReadBuffer( command_queue, dev_fsw1, CL_TRUE, 0, sizeof(float) * nx*ny, host_fsw1, 0, NULL, NULL );
				err|= clEnqueueReadBuffer( command_queue, dev_fse1, CL_TRUE, 0, sizeof(float) * nx*ny, host_fse1, 0, NULL, NULL );
				clFinish(command_queue);
			}
			
			
			for(int iX=0;iX<nx;iX++)
				for(int iY=0;iY<ny;iY++)
				{
					int coor=iX*ny+iY;
					host_rho[coor]=host_fr1[coor]+host_fe1[coor]+host_fn1[coor]+host_fw1[coor]+host_fs1[coor]+
					host_fne1[coor]+host_fnw1[coor]+host_fsw1[coor]+host_fse1[coor];
				}
			std::stringstream imagestream;
			std::stringstream len;
			len<<timecounter;
			imagestream << "height"<<std::string(5-len.str().size(),'0')<<timecounter<<".dat";
			writeFile(imagestream.str(), host_rho, nx,ny);
			
		}

	}

	finish = time(NULL);

	std::cout<<"Time spent is "<<finish-start<<" sec"<<"\n";
			clReleaseMemObject(dev_rho);
			clReleaseMemObject(dev_fx);
			clReleaseMemObject(dev_fy);
	clReleaseMemObject(dev_fr0);
	clReleaseMemObject(dev_fe0);
	clReleaseMemObject(dev_fn0);
	clReleaseMemObject(dev_fw0);
	clReleaseMemObject(dev_fne0);
	clReleaseMemObject(dev_fnw0);
	clReleaseMemObject(dev_fsw0);
	clReleaseMemObject(dev_fse0);
	
	clReleaseMemObject(dev_fr1);
	clReleaseMemObject(dev_fe1);
	clReleaseMemObject(dev_fn1);
	clReleaseMemObject(dev_fw1);
	clReleaseMemObject(dev_fne1);
	clReleaseMemObject(dev_fnw1);
	clReleaseMemObject(dev_fsw1);
	clReleaseMemObject(dev_fse1);
	

    clReleaseProgram(program);
    clReleaseKernel(kernelCollProp);
	clReleaseKernel(kernelFinishProp);
    clReleaseCommandQueue(command_queue);
    clReleaseContext(context);

    return 0;
}


