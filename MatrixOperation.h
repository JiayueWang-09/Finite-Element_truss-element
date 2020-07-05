#ifndef _MATRIXOPERATION_ 
#define _MATRIXOPERATION_ 

void TransformedStiffnessMatrix(double ***Kt,double *EA,double *length,double *angle,int element);

void GlobalStiffnessMatrix(double **Kg,double ***Kt,int **Nodes,int Elements,int NodeNumber);

//void Jacobi(double **Kr,double *u,double *Force,int NodeNumber);

void Gauss(double **K,double *u,double *Force,int n);

void Displacement(double *ur,double *u,int NodeNumber);

void ForceVector(double *u,double *Force,double **K,int NodeNumber);

void StressStrain(double *strain,double *stress,double *StrainEnergy,double *angle,double *EA,double *u,int **Nodes,int Elements,double *length);

#endif