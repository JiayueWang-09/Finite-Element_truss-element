#include <iostream>
using namespace std;

void InitializeArray(double *Array,int m,double InitialValue)
{
	int i;
	for(i=0;i<m;i++)
	{
		Array[i] = InitialValue;
	}
}
void InitializeMatrix(double **Matrix,int m,int n,double InitialValue)
{
	int i,j;
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			Matrix[i][j] = InitialValue;
		}
	}
}
void InputTruss(int elements,double *angle,double *length,double *EA,int **nodes)
{
for (int i=0;i<elements;i++)
	{
		cout<<"Please set the rotation angle for element"<<i+1<<endl;
		cin>>angle[i];
		cout<<"Please set the length for element"<<i+1<<endl;
		cin>>length[i];
		cout<<"Please set EA for element"<<i+1<<endl;
		cin>>EA[i];	
		cout<<"Please set Node number for element"<<i+1<<endl;
		cin>>nodes[i][1];	
		cin>>nodes[i][2];	
	}
}