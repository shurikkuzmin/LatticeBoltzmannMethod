#ifndef DESCRIPTOR_H
#define DESCRIPTOR_H
#include <vector>
#include <cmath>
#include <iostream>
class Descriptor
{
	public:
	static const int d=3;
	static const int NPOP=19;

	static const char cx[NPOP];
    static const char cy[NPOP];
    static const char cz[NPOP];

    static const double qxx[NPOP];
    static const double qyy[NPOP];
    static const double qzz[NPOP];
    static const double qxy[NPOP];
    static const double qxz[NPOP];
    static const double qyz[NPOP];

    //Stencils
    static const double wxx[NPOP];
    static const double wyy[NPOP];
    static const double wzz[NPOP];
    static const double wxy[NPOP];
    static const double wyz[NPOP];
    static const double wzx[NPOP];
    static const double weights[NPOP];

    static const double laplace_stencil[NPOP];
    static const double gradx_stencil[NPOP];
    static const double grady_stencil[NPOP];
    static const double gradz_stencil[NPOP];
    static const int compliment[NPOP];


	public:
	Descriptor(){}
};

const char Descriptor::cx[19]={0,1,-1,0, 0,0, 0,1,-1, 1,-1,0, 0, 0, 0,1,-1, 1,-1};
const char Descriptor::cy[19]={0,0, 0,1,-1,0, 0,1, 1,-1,-1,1,-1, 1,-1,0, 0, 0, 0};
const char Descriptor::cz[19]={0,0, 0,0, 0,1,-1,0, 0, 0, 0,1, 1,-1,-1,1, 1,-1,-1};

const double Descriptor::qxx[19]={-1.0/3.0, 2.0/3.0, 2.0/3.0, -1.0/3.0, -1.0/3.0, -1.0/3.0,
    -1/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, -1.0/3.0, -1.0/3.0, -1.0/3.0,
    -1.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0};

const double Descriptor::qyy[19]={-1.0/3.0, -1.0/3.0, -1.0/3.0, 2.0/3.0, 2.0/3.0, -1.0/3.0,
    -1.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0,
    -1.0/3.0, -1.0/3.0, -1.0/3.0, -1.0/3.0};

const double Descriptor::qzz[19]={-1.0/3.0, -1.0/3.0, -1.0/3.0, -1.0/3.0, -1.0/3.0, 2.0/3.0,
    2.0/3.0, -1.0/3.0, -1.0/3.0, -1.0/3.0, -1.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0,
    2.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0};

const double Descriptor::qxy[19]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
const double Descriptor::qxz[19]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0};
const double Descriptor::qyz[19]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0};

const double Descriptor::wxx[19]={0.0,5.0/12.0,5.0/12.0,-1.0/3.0,-1.0/3.0,-1.0/3.0,-1.0/3.0,
		-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,
		-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0};
const double Descriptor::wyy[19]={0.0,-1.0/3.0,-1.0/3.0,5.0/12.0,5.0/12.0,-1.0/3.0,-1.0/3.0,
		-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,
		1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0};
const double Descriptor::wzz[19]={0.0,-1.0/3.0,-1.0/3.0,-1.0/3.0,-1.0/3.0,5.0/12.0,5.0/12.0,
		1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,
		-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0};
const double Descriptor::wxy[19]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		1.0/4.0,-1.0/4.0,-1.0/4.0,1.0/4.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
const double Descriptor::wyz[19]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,
		1.0/4.0,-1.0/4.0,-1.0/4.0,1.0/4.0,
		0.0,0.0,0.0,0.0};
const double Descriptor::wzx[19]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		1.0/4.0,-1.0/4.0,-1.0/4.0,1.0/4.0};
const double Descriptor::weights[19]={0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,
		1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,
		1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0,1.0/12.0};

const double Descriptor::gradx_stencil[19] = {0.0, 1.0/6.0, -1.0/6.0, 0.0, 0.0, 0.0, 0.0, 1.0/12.0, -1.0/12.0, 1.0/12.0,
                                -1.0/12.0, 0.0, 0.0, 0.0, 0.0, 1.0/12.0, -1.0/12.0, 1.0/12.0, -1.0/12.0};

const double Descriptor::grady_stencil[19] = {0.0, 0.0, 0.0, 1.0/6.0, -1.0/6.0, 0.0, 0.0, 1.0/12.0, 1.0/12.0, -1.0/12.0,
                                -1.0/12.0, 1.0/12.0, -1.0/12.0, 1.0/12.0, -1.0/12.0, 0.0, 0.0, 0.0, 0.0};

const double Descriptor::gradz_stencil[19] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0/6.0, -1.0/6.0, 0.0, 0.0, 0.0, 0.0, 1.0/12.0, 1.0/12.0,
                                -1.0/12.0, -1.0/12.0, 1.0/12.0, 1.0/12.0, -1.0/12.0, -1.0/12.0};

const double Descriptor::laplace_stencil[19]={-4.0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/6.0, 1.0/6.0, 1.0/6.0,
                                1.0/6.0, 1.0/6.0,1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0};

const int Descriptor::compliment[19]={0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15};

#endif

