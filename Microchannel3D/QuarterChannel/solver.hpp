#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include "solver.h"
#include "lattice.h"
#include <mpi.h>
#include <iostream>
#include <sstream>
#include <string>
#include "descriptor.h"

template<typename D>
Solver<D>::Solver(ParamsList params):
	NX(params("NX").value<int>()),
	NY(params("NY").value<int>()),
	NZ(params("NZ").value<int>()),
	rank(mpi_wrapper().get_rank()),
	size(mpi_wrapper().get_size())
{
	//Initialization of ranges
	int zstep=NZ/size;
	int zbegin=rank*zstep;
	int zend=rank*zstep+zstep-1;
	if(rank==size-1)
		zend=NZ-1;
	int NUMLOCAL=NX*NY*(zend-zbegin+1);

	//Creation the necessary objects
	lattice = new Lattice<D>(zbegin,zend,NX,NY,NZ,NUMLOCAL,params);

	//Allocation of memory
	if (rank==0)
	{
		density=new double[NX*NY*NZ];
		phase = new double[NX*NY*NZ];
        velocity=new double[3*NX*NY*NZ];
	}

    //Preparation of the geometry types above and at the bottom
	int zbottom=(zbegin-1+NZ)%NZ;
	int ztop=(zend+1+NZ)%NZ;

}

template<typename D> bool Solver<D>::checkNAN()
{
    if (mpi_wrapper().get_rank()==0)
        return isnan(phase[NZ/2*NX*NY+NY*NX/2+NX/2]);
    else
        return false;
}


template<typename D>
void Solver<D>::putDensity(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double dense)
{
    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    if (iZbegin>zend)
        return;
    if (iZend<zbegin)
        return;
    if (iZend>zend) iZend=zend;
    if (iZbegin<zbegin) iZbegin=zbegin;

    //cout<<"zend="<<zend<<" zbegin="<<zbegin;
    lattice->putDensity(iXbegin,iYbegin,iZbegin-zbegin,iXend,iYend,iZend-zbegin,dense);
}
template<typename D> void Solver<D>::putDensity(double dense)
{
    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    lattice->putDensity(0,0,0,NX-1,NY-1,zend-zbegin,dense);
}

template<typename D> void Solver<D>::putPhase(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double phase)
{
    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    if (iZbegin>zend)
        return;
    if (iZend<zbegin)
        return;
    if (iZend>zend) iZend=zend;
    if (iZbegin<zbegin) iZbegin=zbegin;

    lattice->putPhase(iXbegin,iYbegin,iZbegin-zbegin,iXend,iYend,iZend-zbegin,phase);
}
template<typename D> void Solver<D>::putPhase(double phase)
{
    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    lattice->putPhase(0,0,0,NX-1,NY-1,zend-zbegin,phase);
}

template<typename D> void Solver<D>::putVelocity(int iXbegin,int iYbegin,int iZbegin,int iXend,int iYend,int iZend,double * velocity)
{
    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    if (iZbegin>zend)
        return;
    if (iZend<zbegin)
        return;
    if (iZend>zend) iZend=zend;
    if (iZbegin<zbegin) iZbegin=zbegin;

    lattice->putVelocity(iXbegin,iYbegin,iZbegin-zbegin,iXend,iYend,iZend-zbegin,velocity);
}
template<typename D> void Solver<D>::putVelocity(double * velocity)
{
    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    lattice->putVelocity(0,0,0,NX-1,NY-1,zend-zbegin,velocity);
}


template<typename D> void Solver<D>::getWholeDensity()
{
	if (size==1)
	{
		for(int iCount=0;iCount<NX*NY*NZ;iCount++)
			density[iCount]=lattice->rho[iCount];
		return;
	}
    else
	{
		if (rank==0)
		{
       		for(int iCount=0;iCount<lattice->getNUMLOCAL();iCount++)
                density[iCount]=lattice->rho[iCount];

			int offset=lattice->getNUMLOCAL();
			for (int nproc=1;nproc<size;nproc++)
			{

				MPI_Status status;
				int sizes[2];

				MPI_Recv(sizes,2,MPI_INT,nproc,2,MPI_COMM_WORLD,&status);

				double * temp_array=new double[(sizes[1]-sizes[0]+1)*NX*NY];
				MPI_Recv(temp_array,(sizes[1]-sizes[0]+1)*NX*NY,MPI_DOUBLE,nproc,3,MPI_COMM_WORLD,&status);

 				for(int iCount=0;iCount<(sizes[1]-sizes[0]+1)*NX*NY;iCount++)
					density[iCount+offset]=temp_array[iCount];

				offset=offset+(sizes[1]-sizes[0]+1)*NX*NY;

				delete[] temp_array;
			}
		}
		else
		{
			int sizes[2]={lattice->getZbegin(),lattice->getZend()};
			MPI_Send(sizes,2,MPI_INT,0,2,MPI_COMM_WORLD);
			MPI_Send(lattice->rho,(lattice->getZend()-lattice->getZbegin()+1)*NX*NY,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);

	}
}

template<typename D> void Solver<D>::getWholePhase()
{
	if (size==1)
	{
		for(int iCount=0;iCount<NX*NY*NZ;iCount++)
			phase[iCount]=lattice->phase[iCount];
		return;
	}
    else
	{
		if (rank==0)
		{
       		for(int iCount=0;iCount<lattice->getNUMLOCAL();iCount++)
                phase[iCount]=lattice->phase[iCount];

			int offset=lattice->getNUMLOCAL();
			for (int nproc=1;nproc<size;nproc++)
			{

				MPI_Status status;
				int sizes[2];

				MPI_Recv(sizes,2,MPI_INT,nproc,2,MPI_COMM_WORLD,&status);

				double * temp_array=new double[(sizes[1]-sizes[0]+1)*NX*NY];
				MPI_Recv(temp_array,(sizes[1]-sizes[0]+1)*NX*NY,MPI_DOUBLE,nproc,3,MPI_COMM_WORLD,&status);

 				for(int iCount=0;iCount<(sizes[1]-sizes[0]+1)*NX*NY;iCount++)
					phase[iCount+offset]=temp_array[iCount];

				offset=offset+(sizes[1]-sizes[0]+1)*NX*NY;

				delete[] temp_array;
			}
		}
		else
		{
			int sizes[2]={lattice->getZbegin(),lattice->getZend()};
			MPI_Send(sizes,2,MPI_INT,0,2,MPI_COMM_WORLD);
			MPI_Send(lattice->phase,(lattice->getZend()-lattice->getZbegin()+1)*NX*NY,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);

	}
}

template<typename D> void Solver<D>::getWholeVelocity()
{
	if (size==1)
	{
		for(int iCount=0;iCount<NX*NY*NZ;iCount++)
        {
			velocity[3*iCount]=lattice->ux[iCount];
            velocity[3*iCount+1]=lattice->uy[iCount];
            velocity[3*iCount+2]=lattice->uz[iCount];
        }
		return;
	}
    else
	{
		if (rank==0)
		{
       		for(int iCount=0;iCount<lattice->getNUMLOCAL();iCount++)
            {
                velocity[3*iCount]=lattice->ux[iCount];
                velocity[3*iCount+1]=lattice->uy[iCount];
                velocity[3*iCount+2]=lattice->uz[iCount];
            }

			int offset=3*lattice->getNUMLOCAL();
			for (int nproc=1;nproc<size;nproc++)
			{

				MPI_Status status;
				int sizes[2];

				MPI_Recv(sizes,2,MPI_INT,nproc,2,MPI_COMM_WORLD,&status);

				double * temp_array_x=new double[(sizes[1]-sizes[0]+1)*NX*NY];
				double * temp_array_y=new double[(sizes[1]-sizes[0]+1)*NX*NY];
				double * temp_array_z=new double[(sizes[1]-sizes[0]+1)*NX*NY];
				MPI_Recv(temp_array_x,(sizes[1]-sizes[0]+1)*NX*NY,MPI_DOUBLE,nproc,3,MPI_COMM_WORLD,&status);
				MPI_Recv(temp_array_y,(sizes[1]-sizes[0]+1)*NX*NY,MPI_DOUBLE,nproc,3,MPI_COMM_WORLD,&status);
				MPI_Recv(temp_array_z,(sizes[1]-sizes[0]+1)*NX*NY,MPI_DOUBLE,nproc,3,MPI_COMM_WORLD,&status);


 				for(int iCount=0;iCount<(sizes[1]-sizes[0]+1)*NX*NY;iCount++)
				{
					velocity[3*iCount+offset]=temp_array_x[iCount];
					velocity[3*iCount+1+offset]=temp_array_y[iCount];
					velocity[3*iCount+2+offset]=temp_array_z[iCount];
				}
				offset=offset+3*(sizes[1]-sizes[0]+1)*NX*NY;

				delete[] temp_array_x;
				delete[] temp_array_y;
				delete[] temp_array_z;
			}
		}
		else
		{
			int sizes[2]={lattice->getZbegin(),lattice->getZend()};
			MPI_Send(sizes,2,MPI_INT,0,2,MPI_COMM_WORLD);
			MPI_Send(lattice->ux,(lattice->getZend()-lattice->getZbegin()+1)*NX*NY,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
			MPI_Send(lattice->uy,(lattice->getZend()-lattice->getZbegin()+1)*NX*NY,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
			MPI_Send(lattice->uz,(lattice->getZend()-lattice->getZbegin()+1)*NX*NY,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);

	}
}


template<typename D>  void Solver<D>::writeWholeDensity(std::string name)
{
	getWholeDensity();

	if(rank==0)
	{
        name=name+".vts";

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                for(int counterX=0;counterX<NX;counterX++)
                    points->InsertNextPoint(counterX,counterY,counterZ);


        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfComponents(1);
        data->SetName("Density");



        for(vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
        {
            double * temp_data=&density[i];
            data->InsertNextTupleValue(temp_data);
        }

        vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
        structuredGrid->SetDimensions(NX,NY,NZ);
        structuredGrid->SetPoints(points);
        structuredGrid->GetPointData()->AddArray(data);

        vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
        writer->SetFileName(name.c_str());
        //writer.SetPoints();
        writer->SetInput(structuredGrid);
        writer->Write();
	}
}

template<typename D> void Solver<D>::writeWholePhase(std::string name)
{
    Solver::getWholePhase();
    if (rank==0)
    {
        name=name+".vts";

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                for(int counterX=0;counterX<NX;counterX++)
                    points->InsertNextPoint(counterX,counterY,counterZ);


        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfComponents(1);
        data->SetName("Phase");



        for(vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
        {
            double * temp_data=&phase[i];
            data->InsertNextTupleValue(temp_data);
        }

        vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
        structuredGrid->SetDimensions(NX,NY,NZ);
        structuredGrid->SetPoints(points);
        structuredGrid->GetPointData()->AddArray(data);

        vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
        writer->SetFileName(name.c_str());
        //writer.SetPoints();
        writer->SetInput(structuredGrid);
        writer->Write();
    }
}

template<typename D> void Solver<D>::writeWholeVelocity(std::string name)
{
    Solver::getWholeVelocity();
    if (rank==0)
    {
        name=name+".vts";

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                for(int counterX=0;counterX<NX;counterX++)
                    points->InsertNextPoint(counterX,counterY,counterZ);


        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        data->SetNumberOfComponents(3);
        data->SetName("Velocity");



        for(vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
        {
            double * temp_data=&velocity[3*i];
            data->InsertNextTupleValue(temp_data);
        }

        vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
        structuredGrid->SetDimensions(NX,NY,NZ);
        structuredGrid->SetPoints(points);
        structuredGrid->GetPointData()->AddArray(data);

        vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
        writer->SetFileName(name.c_str());
        writer->SetInput(structuredGrid);
        writer->Write();
    }
}

template<typename D> void Solver<D>::load_file(std::string name)
{
    if (size==1)
	{

        vtkSmartPointer<vtkXMLStructuredGridReader> reader=vtkSmartPointer<vtkXMLStructuredGridReader>::New();
        name=name+".vts";
        reader->SetFileName(name.c_str());
        reader->Update();

        vtkStructuredGrid* data = vtkStructuredGrid::New();
        data->ShallowCopy(reader->GetOutput());

        vtkPointData* point_data = data->GetPointData();
        vtkDataArray* density_data=point_data->GetScalars("Density");
        vtkDataArray* phase_data=point_data->GetScalars("Phase");
        vtkDataArray* velocity_data=point_data->GetVectors("Velocity");

        for (int counter=0;counter<density_data->GetNumberOfTuples();counter++)
        {
            double* density_value=density_data->GetTuple(counter);
            double* phase_value=phase_data->GetTuple(counter);
            double* velocity_value=velocity_data->GetTuple(counter);
            lattice->rho[counter]=density_value[0];
            lattice->phase[counter]=phase_value[0];
            lattice->ux[counter]=velocity_value[0];
            lattice->uy[counter]=velocity_value[1];
            lattice->uz[counter]=velocity_value[2];
        }

      	return;
    }
    else
	{
		if (rank==0)
		{

            vtkSmartPointer<vtkXMLStructuredGridReader> reader=vtkSmartPointer<vtkXMLStructuredGridReader>::New();
            name=name+".vts";
            reader->SetFileName(name.c_str());
            reader->Update();

            vtkStructuredGrid* data = vtkStructuredGrid::New();
            data->ShallowCopy(reader->GetOutput());

            vtkPointData* point_data = data->GetPointData();
            vtkDataArray* density_data=point_data->GetScalars("Density");
            vtkDataArray* phase_data=point_data->GetScalars("Phase");
            vtkDataArray* velocity_data=point_data->GetVectors("Velocity");

       		for(int counter=0;counter<lattice->getNUMLOCAL();counter++)
            {
                double* density_value=density_data->GetTuple(counter);
                double* phase_value=phase_data->GetTuple(counter);
                double* velocity_value=velocity_data->GetTuple(counter);
                lattice->rho[counter]=density_value[0];
                lattice->phase[counter]=phase_value[0];
                lattice->ux[counter]=velocity_value[0];
                lattice->uy[counter]=velocity_value[1];
                lattice->uz[counter]=velocity_value[2];
            }

			int offset=lattice->getNUMLOCAL();
			for (int nproc=1;nproc<size;nproc++)
			{
				MPI_Status status;
				int sizes[2];

				MPI_Recv(sizes,2,MPI_INT,nproc,2,MPI_COMM_WORLD,&status);

				double * temp_density=new double[(sizes[1]-sizes[0]+1)*NX*NY];
				double * temp_phase  =new double[(sizes[1]-sizes[0]+1)*NX*NY];
				double * temp_array_x=new double[(sizes[1]-sizes[0]+1)*NX*NY];
				double * temp_array_y=new double[(sizes[1]-sizes[0]+1)*NX*NY];
				double * temp_array_z=new double[(sizes[1]-sizes[0]+1)*NX*NY];

 				for(int iCount=0;iCount<(sizes[1]-sizes[0]+1)*NX*NY;iCount++)
				{
					double* density_value=density_data->GetTuple(iCount+offset);
                    double* phase_value=phase_data->GetTuple(iCount+offset);
                    double* velocity_value=velocity_data->GetTuple(iCount+offset);

					temp_density[iCount]=density_value[0];
					temp_phase[iCount]=phase_value[0];
					temp_array_x[iCount]=velocity_value[0];
					temp_array_y[iCount]=velocity_value[1];
					temp_array_z[iCount]=velocity_value[2];
				}

				MPI_Send(temp_density,(sizes[1]-sizes[0]+1)*NX*NY,MPI_DOUBLE,nproc,3,MPI_COMM_WORLD);
				MPI_Send(temp_phase,(sizes[1]-sizes[0]+1)*NX*NY,MPI_DOUBLE,nproc,3,MPI_COMM_WORLD);
				MPI_Send(temp_array_x,(sizes[1]-sizes[0]+1)*NX*NY,MPI_DOUBLE,nproc,3,MPI_COMM_WORLD);
				MPI_Send(temp_array_y,(sizes[1]-sizes[0]+1)*NX*NY,MPI_DOUBLE,nproc,3,MPI_COMM_WORLD);
				MPI_Send(temp_array_z,(sizes[1]-sizes[0]+1)*NX*NY,MPI_DOUBLE,nproc,3,MPI_COMM_WORLD);


				offset=offset+(sizes[1]-sizes[0]+1)*NX*NY;

				delete[] temp_density;
				delete[] temp_phase;
				delete[] temp_array_x;
				delete[] temp_array_y;
				delete[] temp_array_z;
			}
		}
		else
		{
			MPI_Status status;
			int sizes[2]={lattice->getZbegin(),lattice->getZend()};
			MPI_Send(sizes,2,MPI_INT,0,2,MPI_COMM_WORLD);

			MPI_Recv(lattice->rho,(lattice->getZend()-lattice->getZbegin()+1)*NX*NY,MPI_DOUBLE,0,3,MPI_COMM_WORLD,&status);
			MPI_Recv(lattice->phase,(lattice->getZend()-lattice->getZbegin()+1)*NX*NY,MPI_DOUBLE,0,3,MPI_COMM_WORLD,&status);
			MPI_Recv(lattice->ux,(lattice->getZend()-lattice->getZbegin()+1)*NX*NY,MPI_DOUBLE,0,3,MPI_COMM_WORLD,&status);
			MPI_Recv(lattice->uy,(lattice->getZend()-lattice->getZbegin()+1)*NX*NY,MPI_DOUBLE,0,3,MPI_COMM_WORLD,&status);
			MPI_Recv(lattice->uz,(lattice->getZend()-lattice->getZbegin()+1)*NX*NY,MPI_DOUBLE,0,3,MPI_COMM_WORLD,&status);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

}


template<typename D> void Solver<D>::writeWholeDensityPhaseVelocity(std::string name)
{

    Solver::getWholeDensity();
    Solver::getWholePhase();
    Solver::getWholeVelocity();
    if (rank==0)
    {
        name=name+".vts";

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                for(int counterX=0;counterX<NX;counterX++)
                    points->InsertNextPoint(counterX,counterY,counterZ);

        vtkSmartPointer<vtkDoubleArray> data_phase = vtkSmartPointer<vtkDoubleArray>::New();
        vtkSmartPointer<vtkDoubleArray> data_density = vtkSmartPointer<vtkDoubleArray>::New();
        vtkSmartPointer<vtkDoubleArray> data_velocity = vtkSmartPointer<vtkDoubleArray>::New();

        data_phase->SetNumberOfComponents(1);
        data_phase->SetName("Phase");
        data_density->SetNumberOfComponents(1);
        data_density->SetName("Density");
        data_velocity->SetNumberOfComponents(3);
        data_velocity->SetName("Velocity");



        for(vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
        {
            double * temp_velocity=&velocity[3*i];
            double * temp_density=&density[i];
            double * temp_phase=&phase[i];
            data_velocity->InsertNextTupleValue(temp_velocity);
            data_density->InsertNextTupleValue(temp_density);
            data_phase->InsertNextTupleValue(temp_phase);
        }

        vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
        structuredGrid->SetDimensions(NX,NY,NZ);
        structuredGrid->SetPoints(points);
        structuredGrid->GetPointData()->AddArray(data_velocity);
        structuredGrid->GetPointData()->AddArray(data_phase);
        structuredGrid->GetPointData()->AddArray(data_density);

        vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
        writer->SetFileName(name.c_str());
        writer->SetInput(structuredGrid);
        writer->Write();
    }

}

template<typename D> void Solver<D>::writeTextXZVelocity(std::string name)
{
    Solver::getWholeVelocity();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterX=0;counterX<NX;counterX++)
                fout << counterX << " " << counterZ << " "<<velocity[3*(counterZ*NX*NY+NY/2*NX+counterX)]<<" "<<
                velocity[3*(counterZ*NX*NY+NY/2*NX+counterX)+1]<<" "<<velocity[3*(counterZ*NX*NY+NY/2*NX+counterX)+2]<<"\n";
    }
}

template<typename D> void Solver<D>::writeTextXYVelocity(std::string name)
{
    Solver::getWholeVelocity();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterY=0;counterY<NY;counterY++)
            for(int counterX=0;counterX<NX;counterX++)
                fout << counterX << " " << counterY << " "<<velocity[3*(NZ/2*NX*NY+counterY*NX+counterX)]<<
                " "<<velocity[3*(NZ/2*NX*NY+counterY*NX+counterX)+1]<<" "<<velocity[3*(NZ/2*NX*NY+counterY*NX+counterX)+2]<<"\n";
    }
}

template<typename D> void Solver<D>::writeTextYZVelocity(std::string name)
{
    Solver::getWholeVelocity();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                fout << counterY << " " << counterZ <<" "<<velocity[3*(counterZ*NX*NY+counterY*NX+NX/2)]<<
                " "<<velocity[3*(counterZ*NX*NY+counterY*NX+NX/2)+1]<<" "<<velocity[3*(counterZ*NX*NY+counterY*NX+NX/2)+2]<<"\n";
    }
}


template<typename D> void Solver<D>::writeTextWholePhase(std::string name)
{
    Solver::getWholePhase();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                for(int counterX=0;counterX<NX;counterX++)
                    fout << counterX << " " << counterY << " " << counterZ << " "<<phase[counterZ*NX*NY+counterY*NX+counterX]<<"\n";
    }
}

template<typename D> void Solver<D>::writeTextWholeDensity(std::string name)
{
    Solver::getWholeDensity();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        fout.precision(10);
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                for(int counterX=0;counterX<NX;counterX++)
                    fout << counterX << " " << counterY << " " << counterZ << " "<<density[counterZ*NX*NY+counterY*NX+counterX]<<"\n";
    }
}



//Write planes of phase
template<typename D> void Solver<D>::writeTextXYPhase(std::string name)
{
    Solver::getWholePhase();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterY=0;counterY<NY;counterY++)
            for(int counterX=0;counterX<NX;counterX++)
                fout << counterX << " " << counterY << " "<<phase[NZ/2*NX*NY+counterY*NX+counterX]<<"\n";
    }
}


template<typename D> void Solver<D>::writeTextXZPhase(std::string name)
{
    Solver::getWholePhase();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterX=0;counterX<NX;counterX++)
                fout << counterX << " " << counterZ << " "<<phase[counterZ*NX*NY+NY/2*NX+counterX]<<"\n";
    }
}


template<typename D> void Solver<D>::writeTextYZPhase(std::string name)
{
    Solver::getWholePhase();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                fout << counterY << " " << counterZ <<" "<<phase[counterZ*NX*NY+counterY*NX+NX/2]<<"\n";
    }
}


//Write planes of density
template<typename D> void Solver<D>::writeTextXYDensity(std::string name)
{
    Solver::getWholeDensity();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterY=0;counterY<NY;counterY++)
            for(int counterX=0;counterX<NX;counterX++)
                fout << counterX << " " << counterY << " "<<density[NZ/2*NX*NY+counterY*NX+counterX]<<"\n";
    }
}



template<typename D> void Solver<D>::writeTextXZDensity(std::string name)
{
    Solver::getWholeDensity();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterX=0;counterX<NX;counterX++)
                fout << counterX << " " << counterZ << " "<<density[counterZ*NX*NY+NY/2*NX+counterX]<<"\n";
    }
}



template<typename D> void Solver<D>::writeTextYZDensity(std::string name)
{
    Solver::getWholeDensity();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                fout << counterY << " " << counterZ <<" "<<density[counterZ*NX*NY+counterY*NX+NX/2]<<"\n";
    }
}


//Write velocity
template<typename D> void Solver<D>::writeTextWholeVelocity(std::string name)
{
    Solver::getWholeVelocity();
    if (rank==0)
    {
    	name=name+".dat";
        std::ofstream fout(name.c_str());
        fout.precision(12);
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                for(int counterX=0;counterX<NX;counterX++)
                    fout << counterX << " " << counterY << " " << counterZ << " "<<velocity[3*(counterZ*NX*NY+counterY*NX+counterX)]<<" "
                        <<velocity[3*(counterZ*NX*NY+counterY*NX+counterX)+1]<<" "<<velocity[3*(counterZ*NX*NY+counterY*NX+counterX)+2]<<"\n";
    }
}


template<typename D> void Solver<D>::writeVTKWholePhase(std::string name)
{
    Solver::getWholePhase();
    if (rank==0)
    {
   		name=name+".vtk";
        std::ofstream fout(name.c_str());
        fout<<"# vtk DataFile Version 3.0\n";
        fout<<"Binary liquid phase field vtk representation\n";
        fout<<"ASCII\n\n";
        fout<<"DATASET STRUCTURED_GRID\n";
        fout<<"DIMENSIONS "<<NX<<" "<<NY<<" "<<NZ<<"\n";
        //fout<<"ORIGIN "<<floor(NX/2)<<" "<<floor(NY/2)<<" "<<floor(NZ/2)<<"\n";
        //fout<<"SPACING 1 1 1\n";
        fout<<"POINTS "<<NX*NY*NZ<<" double\n";
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                for(int counterX=0;counterX<NX;counterX++)
                    fout<<counterX<<" "<<counterY<<" "<<counterZ<<"\n";
        fout<<"\n";
        fout<<"POINT_DATA "<<NX*NY*NZ<<"\n";
        fout<<"SCALARS phase double\n";
        fout<<"LOOKUP_TABLE phase_table\n";
        for(int counterZ=0;counterZ<NZ;counterZ++)
            for(int counterY=0;counterY<NY;counterY++)
                for(int counterX=0;counterX<NX;counterX++)
                    fout<<phase[counterZ*NX*NY+counterY*NX+counterX]<<"\n";

    }
}


template<typename D> void Solver<D>::writeLocalDensity(std::string name)
{
    std::stringstream process;
    process<<"_process"<<rank;
    name=name+process.str();
    lattice->writeDensity(name);
}
template<typename D> void Solver<D>::writeLocalPhase(std::string name)
{
    std::stringstream process;
    process<<"_process"<<rank;
    name=name+process.str();
    lattice->writePhase(name);
}

template<typename D> void Solver<D>::writeLocalTextDensity(std::string name)
{
    std::stringstream process;
    process<<"_process"<<rank;
    name=name+process.str();
    lattice->writeTextDensity(name);
}

template<typename D> void Solver<D>::writeLocalTextPhase(std::string name)
{
    std::stringstream process;
    process<<"_process"<<rank;
    name=name+process.str();
    lattice->writeTextPhase(name);
}

template<typename D> void Solver<D>::updatePhase()
{
    //Exchange phase fields to calculate it properly
    int zbegin=lattice->getZbegin();
    int zend=lattice->getZend();
    MPI_Status status;

	for (int iCount=0; iCount<NX*NY; iCount++)
	{
		lattice->phase_top[iCount]=lattice->phase[(zend-zbegin)*NX*NY+iCount];
		lattice->phase_bottom[iCount]=lattice->phase[iCount];
	}

	if (size!=1)
	{
		if(rank==0)
		{
			MPI_Sendrecv_replace(lattice->phase_top,NX*NY,MPI_DOUBLE,1,1,1,1,MPI_COMM_WORLD,&status);
			MPI_Sendrecv_replace(lattice->phase_bottom,NX*NY,MPI_DOUBLE,mpi_wrapper().get_size()-1,1,mpi_wrapper().get_size()-1,1,MPI_COMM_WORLD,&status);
		}
        else
		{
			MPI_Sendrecv_replace(lattice->phase_bottom,NX*NY,MPI_DOUBLE,mpi_wrapper().get_rank()-1,1,mpi_wrapper().get_rank()-1,1,MPI_COMM_WORLD,&status);
			MPI_Sendrecv_replace(lattice->phase_top,NX*NY,MPI_DOUBLE,(mpi_wrapper().get_rank()+1)%mpi_wrapper().get_size(),1,(mpi_wrapper().get_rank()+1)%mpi_wrapper().get_size(),1,MPI_COMM_WORLD,&status);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else
	{
		for (int iCount=0; iCount<NX*NY; iCount++)
		{
			lattice->phase_top[iCount]=lattice->phase[iCount];
			lattice->phase_bottom[iCount]=lattice->phase[(zend-zbegin)*NX*NY+iCount];
		}

	}


}


template<typename D> void Solver<D>::init()
{
    lattice->updateSymmetric();
    lattice->updateWall();
	updatePhase();
    lattice->init();
	lattice->updateMacro();

}

template<typename D> void Solver<D>::collide_stream()
{
    lattice->updateMacro();
    updatePhase();
    lattice->collide_stream();
    if (size!=1)
    {
        propagate();
    }
    exchangeMatrices();
}

template<typename D> void Solver<D>::propagate()
{
    lattice->preparePopulations();

    MPI_Status status;
	if(mpi_wrapper().get_rank()==0)
	{
		MPI_Sendrecv_replace(lattice->layer_top_f,NX*NY*D::NPOP,MPI_DOUBLE,1,1,1,1,MPI_COMM_WORLD,&status);
		MPI_Sendrecv_replace(lattice->layer_bottom_f,NX*NY*D::NPOP,MPI_DOUBLE,mpi_wrapper().get_size()-1,1,mpi_wrapper().get_size()-1,1,MPI_COMM_WORLD,&status);

		MPI_Sendrecv_replace(lattice->layer_top_g,NX*NY*D::NPOP,MPI_DOUBLE,1,1,1,1,MPI_COMM_WORLD,&status);
		MPI_Sendrecv_replace(lattice->layer_bottom_g,NX*NY*D::NPOP,MPI_DOUBLE,mpi_wrapper().get_size()-1,1,mpi_wrapper().get_size()-1,1,MPI_COMM_WORLD,&status);
	}
	else
	{
		MPI_Sendrecv_replace(lattice->layer_bottom_f,NX*NY*D::NPOP,MPI_DOUBLE,mpi_wrapper().get_rank()-1,1,mpi_wrapper().get_rank()-1,1,MPI_COMM_WORLD,&status);
		MPI_Sendrecv_replace(lattice->layer_top_f,NX*NY*D::NPOP,MPI_DOUBLE,(mpi_wrapper().get_rank()+1)%mpi_wrapper().get_size(),1,(mpi_wrapper().get_rank()+1)%mpi_wrapper().get_size(),1,MPI_COMM_WORLD,&status);

		MPI_Sendrecv_replace(lattice->layer_bottom_g,NX*NY*D::NPOP,MPI_DOUBLE,mpi_wrapper().get_rank()-1,1,mpi_wrapper().get_rank()-1,1,MPI_COMM_WORLD,&status);
		MPI_Sendrecv_replace(lattice->layer_top_g,NX*NY*D::NPOP,MPI_DOUBLE,(mpi_wrapper().get_rank()+1)%mpi_wrapper().get_size(),1,(mpi_wrapper().get_rank()+1)%mpi_wrapper().get_size(),1,MPI_COMM_WORLD,&status);
	}

    MPI_Barrier(MPI_COMM_WORLD);

    lattice->finishPropagation();
}

template<typename D> void Solver<D>::exchangeMatrices()
{
    lattice->exchangeMatrices();
}
