#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "Initialization.h"
#include"MatrixOperation.h"
using namespace std;

void main()
{
	//Constant input and vector matrix initialization//
	int NodeNumber=4,Elements=5;
	int **Nodes;
	double *angle,*length,*EA;
	double **Kg;//K global
	double ***Kt;//K transformed for each element

	angle=(double*)malloc(Elements*sizeof(double));
	length=(double*)malloc(Elements*sizeof(double));
	EA=(double*)malloc(Elements*sizeof(double));
	
	Nodes=(int**)malloc(Elements*sizeof(int*));
	for(int i=0;i<Elements;i++)
	{
		*(Nodes+i)=(int*)malloc(2*sizeof(int));
	}

	Kg=(double**)malloc(2*NodeNumber*sizeof(double*));
	for(int i=0;i<2*NodeNumber;i++)
	{
		*(Kg+i) = (double*)malloc(2*NodeNumber*sizeof(double));
		for(int j=0;j<2*NodeNumber;j++)
		{
			Kg[i][j]=0.0;
		}
	}

	Kt = (double***)malloc(Elements*sizeof(double**));
	for(int i=0;i<Elements;i++)
	{
		*(Kt+i) = (double**)malloc(4*sizeof(double*));
		for(int j=0;j<4;j++)
		{
			*(*(Kt+i)+j) = (double*)malloc(4*sizeof(double));
			for(int k=0;k<4;k++)
			{
				Kt[i][j][k] = 0.0;
			}
		}
}
	for(int i=0;i<Elements;i++)
	{
		*(Kt+i) = (double**)malloc(4*sizeof(double*));
		for(int j=0;j<4;j++)
		{
			*(Kt[i]+j) = (double*)malloc(4*sizeof(double));
		}
	}

	InitializeArray(angle,Elements,0.0);
	InitializeArray(length,Elements,0.0);
	InitializeArray(EA,Elements,0.0);

	InputTruss(Elements,angle,length,EA,Nodes);
	
	
	//Assembly stiffness matrix//

    TransformedStiffnessMatrix(Kt,EA,length,angle,Elements);
	GlobalStiffnessMatrix(Kg,Kt,Nodes,Elements,NodeNumber);

	//Reduced Linear System//
	double *ReducedForce;
	double **Kr;
	ReducedForce=(double*)malloc(NodeNumber*sizeof(double));
	Kr=(double**)malloc(NodeNumber*sizeof(double*));
	for(int i=0;i<NodeNumber;i++)
	{
		*(Kr+i) = (double*)malloc(NodeNumber*sizeof(double));
	}
	InitializeArray(ReducedForce,NodeNumber,0.0);
	InitializeMatrix(Kr,NodeNumber,NodeNumber,0.0);

	for(int i=0;i<NodeNumber;i++)
	{
		for(int j=0;j<NodeNumber;j++)
		{
			Kr[i][j]=Kg[i+2][j+2];
		}
	}
	cout<<"Please set Fx at Node2"<<endl;
	cin>>ReducedForce[0];		
	cout<<"Please set Fy at Node2"<<endl;
	cin>>ReducedForce[1];
	cout<<"Please set Fx at Node3"<<endl;
	cin>>ReducedForce[2];		
	cout<<"Please set Fy at Node3"<<endl;
	cin>>ReducedForce[3];

	//Solving Linear system by Gaussian elimination OR Jacobi //
	double *ur;
	ur=(double*)malloc(NodeNumber*sizeof(double));
	InitializeArray(ur,NodeNumber,0.0);
	//Jacobi(Kr,ur,ReducedForce,NodeNumber);
	Gauss(Kr,ur,ReducedForce,NodeNumber);

	//Calculate displacement and force vector//
	double *u,*Force;
	Force=(double*)malloc(2*NodeNumber*sizeof(double));
	u=(double*)malloc(2*NodeNumber*sizeof(double));
	InitializeArray(Force,2*NodeNumber,0.0);
	InitializeArray(u,2*NodeNumber,0.0);
	Displacement(ur,u,NodeNumber);
	ForceVector(u,Force,Kg,NodeNumber);

	//Calculate strain,strain energy and stress//
	double *strain,*stress,*StrainEnergy;
	strain=(double*)malloc(Elements*sizeof(double));
	stress=(double*)malloc(Elements*sizeof(double));
	StrainEnergy=(double*)malloc(Elements*sizeof(double));
	StressStrain(strain,stress,StrainEnergy,angle,EA,u,Nodes,Elements,length);

	system("Pause");
	//

}