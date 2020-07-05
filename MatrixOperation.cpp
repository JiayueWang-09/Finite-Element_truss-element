#include <math.h>
#include<iostream>
using namespace std;
#include "Initialization.h"

void TransformedStiffnessMatrix(double ***Kt,double *EA,double *length,double *angle,int element)//Transformed elemental stiffness matrix
{
	for(int i=0;i<element;i++)
	{
		angle[i]=angle[i]*atan(1.0)*4/180;
		Kt[i][0][0]=Kt[i][2][2]=EA[i]*cos(angle[i])*cos(angle[i])/length[i];//c^2
		Kt[i][0][1]=Kt[i][1][0]=Kt[i][2][3]=Kt[i][3][2]=EA[i]*cos(angle[i])*sin(angle[i])/length[i];//cs
		Kt[i][1][1]=Kt[i][3][3]=EA[i]*sin(angle[i])*sin(angle[i])/length[i];//s^2
		Kt[i][0][2]=Kt[i][2][0]=-EA[i]*cos(angle[i])*cos(angle[i])/length[i];//-c^2
		Kt[i][0][3]=Kt[i][1][2]=Kt[i][2][1]=Kt[i][3][0]=-EA[i]*cos(angle[i])*sin(angle[i])/length[i];//-cs
		Kt[i][1][3]=Kt[i][3][1]=-EA[i]*sin(angle[i])*sin(angle[i])/length[i];//-s^2

	}
	for(int i=0;i<element;i++)
	{
		cout<<"The Transformed Stiffness Matrix for element "<<i+1<<endl;
		for(int j=0;j<4;j++)
		{
			for(int k=0;k<4;k++)
			{
				cout<<Kt[i][j][k]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;
		cout<<endl;
	}

}

void GlobalStiffnessMatrix(double **Kg,double ***Kt,int **Nodes,int Elements,int NodeNumber)
{
	for(int k=0;k<Elements;k++)
	{
		
		for(int i=0;i<4;i++)
		{
			for(int j=0;j<4;j++)
			{
				if(i<2&&j<2)
				{
					Kg[i+2*Nodes[k][1]][j+2*Nodes[k][1]]+=Kt[k][i][j];
				}
				else if(i>=2&&j<2)
				{
					Kg[i+2*(Nodes[k][2]-1)][j+2*Nodes[k][1]]+=Kt[k][i][j];
				}
				else if(i<2&&j>=2)
				{
					Kg[i+2*Nodes[k][1]][j+2*(Nodes[k][2]-1)]+=Kt[k][i][j];
				}
				else
				{
					Kg[i+2*(Nodes[k][2]-1)][j+2*(Nodes[k][2]-1)]+=Kt[k][i][j];
				}
			}
		}
	}

	cout<<"The Global Stiffness Matrix for truss element"<<endl;
	for(int i=0;i<2*NodeNumber;i++)
	{
		for(int j=0;j<2*NodeNumber;j++)
		{
			cout<<Kg[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
}


//void Jacobi(double **Kr,double *u,double *Force,int NodeNumber)
//{
//	double *u_new;
//	u_new=(double*)malloc(NodeNumber*sizeof(double));
//	double err,sum;
//	err=2.0;
//	while( err>10^-4)
//	{
//		for(int i=0;i<NodeNumber;i++)
//		{
//			for(int j=0;j<NodeNumber;j++)
//			{
//				sum=0.0;
//				if( i!=j)
//				{
//					sum+=Kr[i][j]*u[j];
//				}
//			}
//			u_new[i]=(Force[i]-sum)/Kr[i][i];
//			err+=(u_new[i]-u[i])*(u_new[i]-u[i]);
//		}
//		u=u_new;
//		err=sqrt(err);
//	}
//}

void Gauss(double **K,double *u,double *Force,int n)
{
	double c;
	for(int i=0;i<n-1;i++)
	{
		for(int k=i+1;k<n;k++)
		{
			c=-K[k][i]/K[i][i];
			for(int j=0;j<n;j++)
			{
				K[k][j]+=c*K[i][j];
			}
			Force[k]+=c*Force[i];
		}	
	}

	cout<<"The eliminated linear system is"<<endl;
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			cout<<K[i][j]<<" ";
		}
		cout<<Force[i]<<endl;
	}

	u[n-1]=Force[n-1]/K[n-1][n-1];

	for(int i=n-2;i>=0;i--)
	{
		double sum=0.0;
		for(int j=i+1;j<n;j++)
		{
			sum+=K[i][j]*u[j];
		}
		u[i]=(Force[i]-sum)/K[i][i];
	}

}

void Displacement(double *ur,double *u,int NodeNumber)
{
	for(int i=0;i<2*NodeNumber;i++)
	{
		if(1<i&&i<6)
		{
			u[i]=ur[i-2];
		}
		else
		{
			u[i]=0.0;
		}
	}

	cout<<"The displacement for the truss elements is"<<endl;
	for(int i=0;i<NodeNumber;i++)
	{
		cout<<"u"<<i+1<<"="<<u[2*i]<<"mm"<<endl;
		cout<<"v"<<i+1<<"="<<u[2*i+1]<<"mm"<<endl;
	}
}

void ForceVector(double *u,double *Force,double **K,int NodeNumber)
{
	for(int i=0;i<2*NodeNumber;i++)
	{
		for(int j=0;j<2*NodeNumber;j++)
		{
			Force[i]+=K[i][j]*u[j];
		}
		
	}

	cout<<"The Reaction Force for the truss elements is"<<endl;
	for(int i=0;i<NodeNumber;i++)
	{
		cout<<"Fx"<<i+1<<"="<<Force[2*i]<<"N"<<endl;
		cout<<"Fy"<<i+1<<"="<<Force[2*i+1]<<"N"<<endl;
	}
}

void StressStrain(double* strain,double *stress,double *StrainEnergy,double *angle,double *EA,double *u,int **Nodes,int Elements,double *length)
{
	for(int i=0;i<Elements;i++)
	{
		double E=70000;
		strain[i]=(-(u[2*Nodes[i][1]]*cos(angle[i])+u[2*Nodes[i][1]+1]*sin(angle[i]))+u[2*Nodes[i][2]]*cos(angle[i])+u[2*Nodes[i][2]+1]*sin(angle[i]))*0.001/length[i];
		stress[i]=E*strain[i];
		StrainEnergy[i]=0.5*length[i]*EA[i]*1000*strain[i]*strain[i];
		cout<<"The strain for element "<<i+1<<" is "<<strain[i]<<endl;
		cout<<"The stress for element "<<i+1<<" is "<<stress[i]<<"MPa"<<endl;
		cout<<"The strain energy for element "<<i+1<<" is "<<StrainEnergy[i]<<"N*m"<<endl;

	}
}