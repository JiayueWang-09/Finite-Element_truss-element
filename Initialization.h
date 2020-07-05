#ifndef _INITIALIZATION_ 
#define _INITIALIZATION_ 

void InitializeMatrix(double **Matrix,int m,int n,double InitialValue);

void InitializeArray(double *Array,int m,double InitialValue);

void InputTruss(int elements,double *angle,double *length,double *EA,int **nodes);

#endif