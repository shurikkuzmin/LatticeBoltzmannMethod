#include "lattice.h"
#include "mpi_singleton.h"

template<typename D>
void Lattice<D>::updateMacro()
{
    //First of all we do update for the bulk
    for (int counter=0;counter<NUMLOCAL;counter++)
    {

        double phase_temp=0.0, dense_temp=0.0, ux_temp=0.0, uy_temp=0.0, uz_temp=0.0;
        double * f_temp=&f[counter*D::NPOP];
        double * g_temp=&g[counter*D::NPOP];

        dense_temp=f_temp[0]+f_temp[1]+f_temp[2]+f_temp[3]+f_temp[4]
                +f_temp[5]+f_temp[6]+f_temp[7]+f_temp[8]+f_temp[9]
                +f_temp[10]+f_temp[11]+f_temp[12]+f_temp[13]+f_temp[14]+
                +f_temp[15]+f_temp[16]+f_temp[17]+f_temp[18];
        phase_temp=g_temp[0]+g_temp[1]+g_temp[2]+g_temp[3]+g_temp[4]
                +g_temp[5]+g_temp[6]+g_temp[7]+g_temp[8]+g_temp[9]
                +g_temp[10]+g_temp[11]+g_temp[12]+g_temp[13]+g_temp[14]+
                +g_temp[15]+g_temp[16]+g_temp[17]+g_temp[18];

        ux_temp=f_temp[1]-f_temp[2]+f_temp[7]-f_temp[8]+f_temp[9]-f_temp[10]+f_temp[15]-f_temp[16]+f_temp[17]-f_temp[18];
        uy_temp=f_temp[3]-f_temp[4]+f_temp[7]+f_temp[8]-f_temp[9]-f_temp[10]+f_temp[11]-f_temp[12]+f_temp[13]-f_temp[14];
        uz_temp=f_temp[5]-f_temp[6]+f_temp[11]+f_temp[12]-f_temp[13]-f_temp[14]+f_temp[15]+f_temp[16]-f_temp[17]-f_temp[18];
        //for (int k=0;k<D::NPOP;k++)
        //{
            //dense_temp+=lattice->f[counter*D::NPOP+k];
            //phase_temp+=lattice->g[counter*D::NPOP+k];
            //ux_temp+=lattice->f[counter*D::NPOP+k]*D::cx[k];
            //uy_temp+=lattice->f[counter*D::NPOP+k]*D::cy[k];
            //uz_temp+=lattice->f[counter*D::NPOP+k]*D::cz[k];
        //}

        ux_temp=ux_temp/dense_temp;
        uy_temp=uy_temp/dense_temp;
        uz_temp=uz_temp/dense_temp;

        phase[counter]=phase_temp;
        rho[counter]=dense_temp;
        ux[counter]=ux_temp;
        uy[counter]=uy_temp;
        uz[counter]=uz_temp;
    }

    //Update Symmetric nodes - unnecessary if everything is OK with other nodes
//    for(int iZ=0;iZ<zend-zbegin+1;iZ++)
//    {
//        for(int iX=1;iX<NX-1;iX++)
//        {
//            int counter=iZ*NX*NY+iX;
//            int counter2=iZ*NX*NY+NX+iX;
//
//            rho[counter]=rho[counter2];
//            phase[counter]=phase[counter2];
//        }
//
//        for(int iY=1;iY<NY-1;iY++)
//        {
//            int counter=iZ*NX*NY+iY*NX;
//            int counter2=iZ*NX*NY+iY*NX+1;
//            rho[counter]=rho[counter2];
//            phase[counter]=phase[counter2];
//        }
//
//        int counter=iZ*NX*NY;
//        int counter2=iZ*NX*NY+NX+1;
//        rho[counter]=rho[counter2];
//        phase[counter]=phase[counter2];
//    }


    //Update BB nodes
    for(int iZ=0;iZ<zend-zbegin+1;iZ++)
    {
        for(int iX=0;iX<NX-1;iX++)
        {
            int counter=iZ*NX*NY+(NY-1)*NX+iX;
            int counter2=iZ*NX*NY+(NY-2)*NX+iX;
            rho[counter]=rho_wall;
            phase[counter]=phase[counter2]-phase_gradient;
            ux[counter]=0.0;
            uy[counter]=0.0;
            uz[counter]=0.0;
        }

        for(int iY=0;iY<NY-1;iY++)
        {
            int counter=iZ*NX*NY+iY*NX+NX-1;
            int counter2=iZ*NX*NY+iY*NX+NX-2;
            rho[counter]=rho_wall;
            phase[counter]=phase[counter2]-phase_gradient;
            ux[counter]=0.0;
            uy[counter]=0.0;
            uz[counter]=0.0;
        }

        int counter=iZ*NX*NY+(NY-1)*NX+NX-1;
        int counter2=iZ*NX*NY+(NY-2)*NX+NX-2;
        rho[counter]=rho_wall;
        phase[counter]=phase[counter2]-phase_gradient;
        ux[counter]=0.0;
        uy[counter]=0.0;
        uz[counter]=0.0;
    }

}


template<typename D>
void Lattice<D>::writePhase(std::string name)
{

    name=name+".vts";

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    for(int counterZ=0;counterZ<zend-zbegin+1;counterZ++)
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
    structuredGrid->SetDimensions(NX,NY,zend-zbegin+1);
    structuredGrid->SetPoints(points);
    structuredGrid->GetPointData()->AddArray(data);

    vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetFileName(name.c_str());
    writer->SetInput(structuredGrid);
    writer->Write();
}

template<typename D>
void Lattice<D>::writeTextPhase(std::string name)
{

    name=name+".dat";
    std::ofstream fout(name.c_str());
    for(int counterZ=0;counterZ<zend-zbegin+1;counterZ++)
        for(int counterY=0;counterY<NY;counterY++)
            for(int counterX=0;counterX<NX;counterX++)
                fout << counterX << " " << counterY << " " << counterZ << " "<<phase[counterZ*NX*NY+counterY*NX+counterX]<<"\n";

}

template<typename D>
void Lattice<D>::writeTextDensity(std::string name)
{

    name=name+".dat";
    std::ofstream fout(name.c_str());
    for(int counterZ=0;counterZ<zend-zbegin+1;counterZ++)
        for(int counterY=0;counterY<NY;counterY++)
            for(int counterX=0;counterX<NX;counterX++)
                fout << counterX << " " << counterY << " " << counterZ << " "<<rho[counterZ*NX*NY+counterY*NX+counterX]<<"\n";

}

template<typename D>
void Lattice<D>::writeDensity(std::string name)
{

    name=name+".vts";

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    for(int counterZ=0;counterZ<zend-zbegin+1;counterZ++)
        for(int counterY=0;counterY<NY;counterY++)
            for(int counterX=0;counterX<NX;counterX++)
                points->InsertNextPoint(counterX,counterY,counterZ);


    vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
    data->SetNumberOfComponents(1);
    data->SetName("Density");



    for(vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
    {
        double * temp_data=&rho[i];
        data->InsertNextTupleValue(temp_data);
    }

    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
    structuredGrid->SetDimensions(NX,NY,zend-zbegin+1);
    structuredGrid->SetPoints(points);
    structuredGrid->GetPointData()->AddArray(data);

    vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetFileName(name.c_str());
    //writer.SetPoints();
    writer->SetInput(structuredGrid);
    writer->Write();
}

template<typename D>
void Lattice<D>::init()
{
 	//Initializing bulk nodes
	for(int counter=0;counter<NUMLOCAL;counter++)
	{
        int iZ=counter/(NX*NY);
        int iY=(counter%(NX*NY))/NX;
        int iX=(counter%(NX*NY))%NX;

        if ((iX==0)||(iX==NX-1)||(iY==0)||(iY==NY-1))
            continue;

        double gradx_temp=0.0,grady_temp=0.0,gradz_temp=0.0,laplace_temp=0.0;
        for (int k=0;k<D::NPOP;k++)
        {
            double phase_temp;
            int iZ2=iZ+D::cz[k];
            int iY2=(iY+D::cy[k]+NY)%NY;
            int iX2=(iX+D::cx[k]+NX)%NX;
            if (iZ2==-1)
                phase_temp=phase_bottom[iY2*NX+iX2];
            else if (iZ2==zend-zbegin+1)
				phase_temp=phase_top[iY2*NX+iX2];
			else
				phase_temp=phase[iZ2*NX*NY+iY2*NX+iX2];

            gradx_temp+=D::gradx_stencil[k]*phase_temp;
            grady_temp+=D::grady_stencil[k]*phase_temp;
            gradz_temp+=D::gradz_stencil[k]*phase_temp;
            laplace_temp+=D::laplace_stencil[k]*phase_temp;
        }

        double phase_temp=phase[counter];
        double dense_temp=rho[counter];
        double ux_temp=ux[counter];
        double uy_temp=uy[counter];
        double uz_temp=uz[counter];

        double feq;
        double geq;
        double sum=0.0;
        double sum_phase=0.0;
        double phase_square=phase_temp*phase_temp;
        double pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
        double chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);

        for (int k=1; k<D::NPOP; k++)
        {
            feq=D::weights[k]*(pressure_bulk+dense_temp*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp)
						+1.5*dense_temp*((D::cx[k]*D::cx[k]-1.0/3.0)*ux_temp*ux_temp+(D::cy[k]*D::cy[k]-1.0/3.0)*uy_temp*uy_temp+(D::cz[k]*D::cz[k]-1.0/3.0)*uz_temp*uz_temp
										 +2.0*(ux_temp*uy_temp*D::cx[k]*D::cy[k]+ux_temp*uz_temp*D::cx[k]*D::cz[k]+uy_temp*uz_temp*D::cy[k]*D::cz[k])))
                +kconst*(D::wxx[k]*gradx_temp*gradx_temp+D::wyy[k]*grady_temp*grady_temp+D::wzz[k]*gradz_temp*gradz_temp
                        +D::wxy[k]*gradx_temp*grady_temp+D::wyz[k]*grady_temp*gradz_temp+D::wzx[k]*gradx_temp*gradz_temp);
            geq=D::weights[k]*(chemical_pot+phase_temp*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp)
						+1.5*phase_temp*((D::cx[k]*D::cx[k]-1.0/3.0)*ux_temp*ux_temp+(D::cy[k]*D::cy[k]-1.0/3.0)*uy_temp*uy_temp+(D::cz[k]*D::cz[k]-1.0/3.0)*uz_temp*uz_temp
										 +2.0*(ux_temp*uy_temp*D::cx[k]*D::cy[k]+ux_temp*uz_temp*D::cx[k]*D::cz[k]+uy_temp*uz_temp*D::cy[k]*D::cz[k])));
            sum+=feq;
            sum_phase+=geq;
            f[counter*D::NPOP+k]=feq;
            g[counter*D::NPOP+k]=geq;
        }
        f[counter*D::NPOP]=dense_temp-sum;
        g[counter*D::NPOP]=phase_temp-sum_phase;
	}

    //Initialize Symmetric nodes
    for(int iZ=0;iZ<zend-zbegin+1;iZ++)
    {
        for(int iX=1;iX<NX-1;iX++)
        {
            int counter=iZ*NX*NY+iX;
            int counter2=iZ*NX*NY+NX+iX;
            for(int k=0;k<D::NPOP;k++)
            {
                f[counter*D::NPOP+symmetricy[k]]=f[counter2*D::NPOP+k];
                g[counter*D::NPOP+symmetricy[k]]=g[counter2*D::NPOP+k];
            }

        }

        for(int iY=1;iY<NY-1;iY++)
        {
            int counter=iZ*NX*NY+iY*NX;
            int counter2=iZ*NX*NY+iY*NX+1;
            for(int k=0;k<D::NPOP;k++)
            {
                f[counter*D::NPOP+symmetricx[k]]=f[counter2*D::NPOP+k];
                g[counter*D::NPOP+symmetricx[k]]=g[counter2*D::NPOP+k];
            }
        }

        int counter=iZ*NX*NY;
        int counter2=iZ*NX*NY+NX+1;

        for(int k=0;k<D::NPOP;k++)
        {
            f[counter*D::NPOP+symmetricxy[k]]=f[counter2*D::NPOP+k];
            g[counter*D::NPOP+symmetricxy[k]]=g[counter2*D::NPOP+k];
        }

    }




    //Initializing BB nodes
    for(int iZ=0;iZ<zend-zbegin+1;iZ++)
    {
        for(int iX=0;iX<NX-1;iX++)
        {
            int counter=iZ*NX*NY+(NY-1)*NX+iX;
            double sum=0.0;
            double sum_phase=0.0;

            for(int k=1;k<D::NPOP;k++)
            {
                f[counter*D::NPOP+k]=D::weights[k]*rho[counter]/3.0;
                g[counter*D::NPOP+k]=D::weights[k]*phase[counter]/3.0;
                sum=sum+f[counter*D::NPOP+k];
                sum_phase=sum_phase+g[counter*D::NPOP+k];
            }
            f[counter*D::NPOP]=rho[counter]-sum;
            g[counter*D::NPOP]=phase[counter]-sum_phase;
        }

        for(int iY=0;iY<NY-1;iY++)
        {
            int counter=iZ*NX*NY+iY*NX+NX-1;
            double sum=0.0;
            double sum_phase=0.0;

            for(int k=1;k<D::NPOP;k++)
            {
                f[counter*D::NPOP+k]=D::weights[k]*rho[counter]/3.0;
                g[counter*D::NPOP+k]=D::weights[k]*phase[counter]/3.0;
                sum=sum+f[counter*D::NPOP+k];
                sum_phase=sum_phase+g[counter*D::NPOP+k];
            }

            f[counter*D::NPOP]=rho[counter]-sum;
            g[counter*D::NPOP]=phase[counter]-sum_phase;

        }

        int counter=iZ*NX*NY+(NY-1)*NX+NX-1;
        double sum=0.0;
        double sum_phase=0.0;

        for(int k=1;k<D::NPOP;k++)
        {
            f[counter*D::NPOP+k]=D::weights[k]*rho[counter];
            g[counter*D::NPOP+k]=D::weights[k]*phase[counter];
            sum=sum+f[counter*D::NPOP+k];
            sum_phase=sum_phase+g[counter*D::NPOP+k];
        }

        f[counter*D::NPOP]=rho[counter]-sum;
        g[counter*D::NPOP]=phase[counter]-sum_phase;

    }


}

template<typename D>
void Lattice<D>::collide_stream()
{
	//Collision and streaming for the bulk
	for (int counter=0;counter<NUMLOCAL;counter++)
	{
        int iZ=counter/(NX*NY);
        int iY=(counter%(NX*NY))/NX;
        int iX=(counter%(NX*NY))%NX;

        if ((iX==0)||(iX==NX-1)||(iY==0)||(iY==NY-1))
            continue;

        double gradx_temp=0.0,grady_temp=0.0,gradz_temp=0.0,laplace_temp=0.0;

        for (int k=0;k<D::NPOP;k++)
        {
            int iZ2=iZ+D::cz[k];
            int iY2=(iY+D::cy[k]+NY)%NY;
            int iX2=(iX+D::cx[k]+NX)%NX;
            double phase_temp;
            if (iZ2==-1)
                phase_temp=phase_bottom[iY2*NX+iX2];
            else if (iZ2==zend-zbegin+1)
                phase_temp=phase_top[iY2*NX+iX2];
            else
				phase_temp=phase[iZ2*NX*NY+iY2*NX+iX2];

            gradx_temp+=D::gradx_stencil[k]*phase_temp;
            grady_temp+=D::grady_stencil[k]*phase_temp;
            gradz_temp+=D::gradz_stencil[k]*phase_temp;
            laplace_temp+=D::laplace_stencil[k]*phase_temp;
        }

        //Force addition
        ux[counter]+=force_x/(2.0*rho[counter]);
        uy[counter]+=force_y/(2.0*rho[counter]);
        uz[counter]+=force_z/(2.0*rho[counter]);

        //std::cout<<"force_z="<<force_z;

        //Reconstruction of equilibrium populations
        double phase_temp=phase[counter];
        double dense_temp=rho[counter];
        double ux_temp=ux[counter];
        double uy_temp=uy[counter];
        double uz_temp=uz[counter]; //+force_z/(2.0*rho[counter]);

        double phase_square=phase_temp*phase_temp;
        double pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
        double chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);

        double sum=0.0;
        double sum_phase=0.0;
        double ux_temp_sq=ux_temp*ux_temp;
        double uy_temp_sq=uy_temp*uy_temp;
        double uz_temp_sq=uz_temp*uz_temp;
        double uxuy_temp=ux_temp*uy_temp;
        double uxuz_temp=ux_temp*uz_temp;
        double uyuz_temp=uy_temp*uz_temp;

        double gradx_temp_sq=gradx_temp*gradx_temp;
        double grady_temp_sq=grady_temp*grady_temp;
        double gradz_temp_sq=gradz_temp*gradz_temp;
        double gradxy_temp=gradx_temp*grady_temp;
        double gradxz_temp=gradx_temp*gradz_temp;
        double gradyz_temp=grady_temp*gradz_temp;

        double feq[D::NPOP], geq[D::NPOP];
        for (int k=1; k<D::NPOP; k++)
        {
            feq[k]=D::weights[k]*(pressure_bulk+dense_temp*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp)
						+1.5*dense_temp*(D::qxx[k]*ux_temp_sq+D::qyy[k]*uy_temp_sq+D::qzz[k]*uz_temp_sq
										 +2.0*(uxuy_temp*D::qxy[k]+uxuz_temp*D::qxz[k]+uyuz_temp*D::qyz[k])))
                +kconst*(D::wxx[k]*gradx_temp_sq+D::wyy[k]*grady_temp_sq+D::wzz[k]*gradz_temp_sq
                        +D::wxy[k]*gradxy_temp+D::wyz[k]*gradyz_temp+D::wzx[k]*gradxz_temp);
            geq[k]=D::weights[k]*(chemical_pot+phase_temp*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp)
						+1.5*phase_temp*(D::qxx[k]*ux_temp_sq+D::qyy[k]*uy_temp_sq+D::qzz[k]*uz_temp_sq
										 +2.0*(uxuy_temp*D::qxy[k]+uxuz_temp*D::qxz[k]+uyuz_temp*D::qyz[k])));

            sum+=feq[k];
            sum_phase+=geq[k];
        }

        feq[0]=dense_temp-sum;
        geq[0]=phase_temp-sum_phase;

        //Omega calculation
        double tau_temp=tau_gas+(phase_temp+1.0)/2.0*(tau_liq-tau_gas);
        double omega_rho=1.0/tau_temp;
        double omega_phi=1.0/tau_phi;

        //Force calculation
        double sum_force=0.0;
        double force_pop[D::NPOP];
        for (int k=1;k<D::NPOP;k++)
        {
            force_pop[k]=(1.0-0.5*omega_rho)*D::weights[k]*(force_x*((D::cx[k]-ux_temp)+3.0*D::cx[k]*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp))+
                force_y*((D::cy[k]-uy_temp)+3.0*D::cy[k]*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp))+
                force_z*((D::cz[k]-uz_temp)+3.0*D::cz[k]*(D::cx[k]*ux_temp+D::cy[k]*uy_temp+D::cz[k]*uz_temp)));

            sum_force+=force_pop[k];
        }
        force_pop[0]=-sum_force;

        //Collision procedure
        for (int k=0; k<D::NPOP; k++)
        {
            //cout<<"Increment="<<-omega*(lattice->f[counter*NPOP+k]-feq[k])<<"\n";
            f[counter*D::NPOP+k]+=-omega_rho*(f[counter*D::NPOP+k]-feq[k])+force_pop[k];
            g[counter*D::NPOP+k]+=-omega_phi*(g[counter*D::NPOP+k]-geq[k]);
        }

        //Streaming procedure
        for (int k=0; k<D::NPOP; k++)
        {
            int iZ2=(iZ+D::cz[k]+zend-zbegin+1)%(zend-zbegin+1);
            int iY2=(iY+D::cy[k]+NY)%NY;
            int iX2=(iX+D::cx[k]+NX)%NX;

            f2[(iZ2*NX*NY+iY2*NX+iX2)*D::NPOP+k]=f[counter*D::NPOP+k];
            g2[(iZ2*NX*NY+iY2*NX+iX2)*D::NPOP+k]=g[counter*D::NPOP+k];
        }

	}

	//Symmetric nodes plus update velocity which is changed
	for(int iZ=0;iZ<zend-zbegin+1;iZ++)
    {
        for(int iX=1;iX<NX-1;iX++)
        {
            int counter=iZ*NX*NY+iX;
            int counter2=iZ*NX*NY+NX+iX;
            int iY=0;

            ux[counter]+=force_x/(2*rho[counter]);
            uy[counter]+=force_y/(2*rho[counter]);
            uz[counter]+=force_z/(2*rho[counter]);

            for(int k=0;k<D::NPOP;k++)
            {
                f[counter*D::NPOP+symmetricy[k]]=f[counter2*D::NPOP+k];
                g[counter*D::NPOP+symmetricy[k]]=g[counter2*D::NPOP+k];
            }

            for (int k=0; k<D::NPOP; k++)
            {
                int iZ2=(iZ+D::cz[k]+zend-zbegin+1)%(zend-zbegin+1);
                int iY2=(iY+D::cy[k]+NY)%NY;
                int iX2=(iX+D::cx[k]+NX)%NX;
                f2[(iZ2*NX*NY+iY2*NX+iX2)*D::NPOP+k]=f[counter*D::NPOP+k];
                g2[(iZ2*NX*NY+iY2*NX+iX2)*D::NPOP+k]=g[counter*D::NPOP+k];
            }


        }

        for(int iY=1;iY<NY-1;iY++)
        {
            int counter=iZ*NX*NY+iY*NX;
            int counter2=iZ*NX*NY+iY*NX+1;
            int iX=0;

            ux[counter]+=force_x/(2*rho[counter]);
            uy[counter]+=force_y/(2*rho[counter]);
            uz[counter]+=force_z/(2*rho[counter]);

            for(int k=0;k<D::NPOP;k++)
            {
                f[counter*D::NPOP+symmetricx[k]]=f[counter2*D::NPOP+k];
                g[counter*D::NPOP+symmetricx[k]]=g[counter2*D::NPOP+k];
            }

            for (int k=0; k<D::NPOP; k++)
            {
                int iZ2=(iZ+D::cz[k]+zend-zbegin+1)%(zend-zbegin+1);
                int iY2=(iY+D::cy[k]+NY)%NY;
                int iX2=(iX+D::cx[k]+NX)%NX;
                f2[(iZ2*NX*NY+iY2*NX+iX2)*D::NPOP+k]=f[counter*D::NPOP+k];
                g2[(iZ2*NX*NY+iY2*NX+iX2)*D::NPOP+k]=g[counter*D::NPOP+k];
            }

        }

        int counter=iZ*NX*NY;
        int counter2=iZ*NX*NY+NX+1;

        ux[counter]+=force_x/(2*rho[counter]);
        uy[counter]+=force_y/(2*rho[counter]);
        uz[counter]+=force_z/(2*rho[counter]);


        for(int k=0;k<D::NPOP;k++)
        {
            f[counter*D::NPOP+symmetricxy[k]]=f[counter2*D::NPOP+k];
            g[counter*D::NPOP+symmetricxy[k]]=g[counter2*D::NPOP+k];
        }

        int iX=0;
        int iY=0;
        for (int k=0; k<D::NPOP; k++)
        {
            int iZ2=(iZ+D::cz[k]+zend-zbegin+1)%(zend-zbegin+1);
            int iY2=(iY+D::cy[k]+NY)%NY;
            int iX2=(iX+D::cx[k]+NX)%NX;
            f2[(iZ2*NX*NY+iY2*NX+iX2)*D::NPOP+k]=f[counter*D::NPOP+k];
            g2[(iZ2*NX*NY+iY2*NX+iX2)*D::NPOP+k]=g[counter*D::NPOP+k];
        }
    }


	//We perform wet collision as far it's easier to track
	for(int iZ=0;iZ<zend-zbegin+1;iZ++)
    {
        for(int iX=0;iX<NX-1;iX++)
        {
            int counter=iZ*NX*NY+(NY-2)*NX+iX;
            int iY=NY-2;
            for(int k=0;k<D::NPOP;k++)
            {
                int iZ2=(iZ+D::cz[k]+zend-zbegin+1)%(zend-zbegin+1);
                int iY2=(iY+D::cy[k]+NY)%NY;
                int iX2=(iX+D::cx[k]+NX)%NX;
                if (iY2==NY-1)
                {
                    f2[counter*D::NPOP+D::compliment[k]]=f[counter*D::NPOP+k];
                    g2[counter*D::NPOP+D::compliment[k]]=g[counter*D::NPOP+k];
                }
            }

        }

        for(int iY=0;iY<NY-1;iY++)
        {
            int counter=iZ*NX*NY+iY*NX+NX-2;
            int iX=NX-2;

            for(int k=0;k<D::NPOP;k++)
            {
                int iZ2=(iZ+D::cz[k]+zend-zbegin+1)%(zend-zbegin+1);
                int iY2=(iY+D::cy[k]+NY)%NY;
                int iX2=(iX+D::cx[k]+NX)%NX;
                if (iX2==NX-1)
                {
                    f2[counter*D::NPOP+D::compliment[k]]=f[counter*D::NPOP+k];
                    g2[counter*D::NPOP+D::compliment[k]]=g[counter*D::NPOP+k];
                }
            }

        }

        int counter=iZ*NX*NY+(NY-2)*NX+NX-2;
        int iY=NY-2;
        int iX=NX-2;
        for(int k=0;k<D::NPOP;k++)
        {
            int iZ2=(iZ+D::cz[k]+zend-zbegin+1)%(zend-zbegin+1);
            int iY2=(iY+D::cy[k]+NY)%NY;
            int iX2=(iX+D::cx[k]+NX)%NX;
            int counter2=iZ2*NX*NY+iY2*NX+iX2;
            if ((iY2==NY-1)||(iX2==NX-1))
            {
                f2[counter*D::NPOP+D::compliment[k]]=f[counter*D::NPOP+k];
                g2[counter*D::NPOP+D::compliment[k]]=g[counter*D::NPOP+k];
            }
        }
    }

    //Update symmetric after-collision nodes in order to have proper macroscopic parameters

	//Symmetric nodes
	for(int iZ=0;iZ<zend-zbegin+1;iZ++)
    {
        for(int iX=1;iX<NX-1;iX++)
        {
            int counter=iZ*NX*NY+iX;
            int counter2=iZ*NX*NY+NX+iX;
            int iY=0;

            for(int k=0;k<D::NPOP;k++)
            {
                f2[counter*D::NPOP+symmetricy[k]]=f2[counter2*D::NPOP+k];
                g2[counter*D::NPOP+symmetricy[k]]=g2[counter2*D::NPOP+k];
            }

         }

        for(int iY=1;iY<NY-1;iY++)
        {
            int counter=iZ*NX*NY+iY*NX;
            int counter2=iZ*NX*NY+iY*NX+1;
            int iX=0;
            for(int k=0;k<D::NPOP;k++)
            {
                f2[counter*D::NPOP+symmetricx[k]]=f2[counter2*D::NPOP+k];
                g2[counter*D::NPOP+symmetricx[k]]=g2[counter2*D::NPOP+k];
            }

        }

        int counter=iZ*NX*NY;
        int counter2=iZ*NX*NY+NX+1;

        for(int k=0;k<D::NPOP;k++)
        {
            f2[counter*D::NPOP+symmetricxy[k]]=f2[counter2*D::NPOP+k];
            g2[counter*D::NPOP+symmetricxy[k]]=g2[counter2*D::NPOP+k];
        }

    }

}

template<typename D>
void Lattice<D>::preparePopulations()
{
   	for (int iCount=0; iCount<NX*NY; iCount++)
	{
		for (int iPop=0; iPop<D::NPOP; iPop++)
		{
			layer_top_f[iCount*D::NPOP+iPop]=f[NX*NY*(zend-zbegin)*D::NPOP+iCount*D::NPOP+iPop];
			layer_bottom_f[iCount*D::NPOP+iPop]=f[iCount*D::NPOP+iPop];
			layer_top_g[iCount*D::NPOP+iPop]=g[NX*NY*(zend-zbegin)*D::NPOP+iCount*D::NPOP+iPop];
			layer_bottom_g[iCount*D::NPOP+iPop]=g[iCount*D::NPOP+iPop];
		}
	}
}

template<typename D> void Lattice<D>::updateWall()
{

    putDensity(0,NY-1,0,NX-1,NY-1,zend-zbegin,rho_wall);
    putDensity(NX-1,0,0,NX-1,NY-1,zend-zbegin,rho_wall);

    for(int iZ=0;iZ<zend-zbegin+1;iZ++)
    {
        int iY=NY-1;
        for(int iX=0;iX<NX-2;iX++)
        {

            double phase_temp=phase[iZ*NX*NY+(iY-1)*NX+iX]-phase_gradient;
            putPhase(iX,iY,iZ,iX,iY,iZ,phase_temp);
        }

        int iX=NX-1;
        for(int iY=0;iY<NY-2;iY++)
        {
            double phase_temp=phase[iZ*NX*NY+iY*NX+iX-1]-phase_gradient;
            putPhase(iX,iY,iZ,iX,iY,iZ,phase_temp);
        }
        //Corner node
        putPhase(NX-1,NY-1,iZ,NX-1,NY-1,iZ,phase[iZ*NX*NY+(NY-2)*NX+NX-2]-phase_gradient);
    }

}

template<typename D> void Lattice<D>::updateSymmetric()
{

    for(int iZ=0;iZ<zend-zbegin+1;iZ++)
    {
        int iY=0;
        for(int iX=1;iX<NX-2;iX++)
        {
            double phase_temp=phase[iZ*NX*NY+NX+iX];
            double rho_temp=rho[iZ*NX*NY+NX+iX];
            putPhase(iX,iY,iZ,iX,iY,iZ,phase_temp);
            putDensity(iX,iY,iZ,iX,iY,iZ,rho_temp);
        }

        int iX=0;
        for(int iY=1;iY<NY-2;iY++)
        {
            double phase_temp=phase[iZ*NX*NY+iY*NX+1];
            double rho_temp=rho[iZ*NX*NY+iY*NX+1];
            putPhase(iX,iY,iZ,iX,iY,iZ,phase_temp);
            putDensity(iX,iY,iZ,iX,iY,iZ,rho_temp);
        }
        //Corner node
        putPhase(0,0,iZ,0,0,iZ,phase[iZ*NX*NY+NX+1]);
        putDensity(0,0,iZ,0,0,iZ,rho[iZ*NX*NY+NX+1]);
    }

}

template<typename D>
void Lattice<D>::updateSymmetricPopulations()
{
    //Symmetric nodes
	for(int iZ=0;iZ<zend-zbegin+1;iZ++)
    {
        for(int iX=1;iX<NX-1;iX++)
        {
            int counter=iZ*NX*NY+iX;
            int counter2=iZ*NX*NY+NX+iX;
            int iY=0;

            for(int k=0;k<D::NPOP;k++)
            {
                f2[counter*D::NPOP+symmetricy[k]]=f2[counter2*D::NPOP+k];
                g2[counter*D::NPOP+symmetricy[k]]=g2[counter2*D::NPOP+k];
            }

         }

        for(int iY=1;iY<NY-1;iY++)
        {
            int counter=iZ*NX*NY+iY*NX;
            int counter2=iZ*NX*NY+iY*NX+1;
            int iX=0;
            for(int k=0;k<D::NPOP;k++)
            {
                f2[counter*D::NPOP+symmetricx[k]]=f2[counter2*D::NPOP+k];
                g2[counter*D::NPOP+symmetricx[k]]=g2[counter2*D::NPOP+k];
            }

        }

        int counter=iZ*NX*NY;
        int counter2=iZ*NX*NY+NX+1;

        for(int k=0;k<D::NPOP;k++)
        {
            f2[counter*D::NPOP+symmetricxy[k]]=f2[counter2*D::NPOP+k];
            g2[counter*D::NPOP+symmetricxy[k]]=g2[counter2*D::NPOP+k];
        }

    }

}


template<typename D>
void Lattice<D>::finishPropagation()
{
   	for (int iCount=0; iCount<NX*NY; iCount++)
	{
        int iY=iCount/NX;
        int iX=iCount%NX;
        if ((iX==NX-1)||(iY==NY-1))
            continue;

        for(int k=0;k<D::NPOP;k++)
        {
			int iY2=(iY+D::cy[k]+NY)%NY;
			int iX2=(iX+D::cx[k]+NX)%NX;
            //cout<<mpi_wrapper().get_rank()<<"\n";
            if (D::cz[k]<0)
            {
                f2[((zend-zbegin)*NX*NY+iY2*NX+iX2)*D::NPOP+k]=layer_top_f[(iY*NX+iX)*D::NPOP+k];
                g2[((zend-zbegin)*NX*NY+iY2*NX+iX2)*D::NPOP+k]=layer_top_g[(iY*NX+iX)*D::NPOP+k];
            }

            if (D::cz[k]>0)
            {
                    f2[(iY2*NX+iX2)*D::NPOP+k]=layer_bottom_f[(iY*NX+iX)*D::NPOP+k];
                    g2[(iY2*NX+iX2)*D::NPOP+k]=layer_bottom_g[(iY*NX+iX)*D::NPOP+k];
            }
        }

	}

	updateSymmetricPopulations();
}

template<typename D>
void Lattice<D>::exchangeMatrices()
{
		double* pointer;

		pointer=f;
		f=f2;
		f2=pointer;

		pointer=g;
		g=g2;
		g2=pointer;

}
