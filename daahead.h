#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
void MergeSort(double* A, int n);
void Merge(double *A,int n,int m);
double mode(double *A,int n);
int LomutoPartition(double *A,int l,int r);
double QuickSelect(double *A,int l,int r,int k);
double median(double *A,int n);
double mean(double *A,int n);
double stddev(double *A,int n);
double hm(double *arr, int n);
void getCofactor(float **A, float **temp, int p, int q, int n);
float determinant(float **A, int n);
void adjoint(float **A,float **adj,int n);
void inverse(float **A, float **inv,int n);
void strassen(float *A,float *B,float *C, int m, int n);
void transpose(float **A,float **T,int n);
void two_to_one(float **A,float *b,int n);
void one_to_two(float *b,float **A,int n);
void leastsqrfit(float **A,float **B,int n,int m);
