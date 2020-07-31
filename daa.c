#include"daahead.h"

void MergeSort(double* A, int n)
{
	int m; 
	if(n<=1)
		return;
	m = n/2;
	MergeSort(A,m);
	MergeSort(A+m,n-m);
	Merge(A,n,m);
}

void Merge(double *A,int n,int m)
{
	int i,k,j = m;
	double *B = (double *)malloc(sizeof(double) * n);
	i=k=0;
	while(i<m && j<n)
	{
		if(A[i]<=A[j])
			B[k++] = A[i++];
		else
			B[k++] = A[j++];
	}
	if(j == n)
	{
		while(i<m && k<n)
			B[k++] = A[i++];
	}
	else
	{
		while(j<n && k<n)
			B[k++] = A[j++];
	}
	for(i=0;i<n;i++)
		A[i] = B[i];
	free(B);
}

double mode(double *A,int n)
{
	int i,mf,rl;
	double rv,mv;
	MergeSort(A,n);
	i = mf = 0;
	while(i<n)
	{
		rl = 1;
		rv = A[i];
		while(i+rl < n && A[i+rl] == rv)
			rl+=1;
		if(rl>mf)
		{
			mf = rl;
			mv = rv;
		}
		i+=rl;
	}
	return mv;
}

int LomutoPartition(double *A,int l,int r)
{
	double p=A[l],t;
	int s=l;
	for(int i=l+1;i<=r;i++)
	{
		if(A[i]<p)
		{
			s += 1;
			t = A[s];
			A[s] = A[i];
			A[i] = t;
		}
	}
	t = A[l];
	A[l] = A[s];
	A[s] = t;
	return s;
}

double QuickSelect(double *A,int l,int r,int k)
{
	int s = LomutoPartition(A,l,r);
	if(s == l+k-1)
		return A[s];
	if(s > l+k-1)
		return QuickSelect(A,l,s-1,k);
	return QuickSelect(A,s+1,r,k-s-1);
}

double median(double *A,int n)
{
	if(n%2 != 0)
		return QuickSelect(A,0,n-1,n/2 + 1);
	double B[n];
	for(int i=0;i<n;i++)
		B[i]=A[i];
	return (QuickSelect(A,0,n-1,n/2) + QuickSelect(B,0,n-1,n/2 + 1))/2;
}

double mean(double *A,int n)
{
	int sum = 0;
	for(int i=0;i<n;i++)
		sum+=A[i];
	return sum/n;
}

double stddev(double *A,int n)
{
	double sumofsqr = 0;
	for(int i=0;i<n;i++)
		sumofsqr+=(A[i]*A[i]);
	double m = mean(A,n);
	return sqrt(sumofsqr/(double)n - m*m);
}

double hm(double *arr, int n) 
{ 
    double sum = 0; 
    for (int i = 0; i < n; i++) 
        sum = sum + (double)1 / arr[i];   
    return (double)n/sum; 
} 

void getCofactor(float **A, float **temp, int p, int q, int n) 
{ 
    int i = 0, j = 0;  
    for (int row = 0; row < n; row++) 
    { 
        for (int col = 0; col < n; col++) 
        { 
            if (row != p && col != q) 
            { 
                temp[i][j++] = A[row][col]; 
                if (j == n - 1) 
                { 
                    j = 0; 
                    i++; 
                } 
            } 
        } 
    } 
} 
  
float determinant(float **A,int n)
{ 
    float D = 1,t,B[n][n]; 
    if (n == 1) 
        return A[0][0]; 
    for (int i=0; i<n; i++) 
    { 
        for (int j=0; j<n; j++)
        	B[i][j] = A[i][j];
    }
	int pivotrow,i,j,k,sgn=1;
    for(i=0;i<n-1;i++)
    {
    	pivotrow = i;
    	for(j=i+1;j<n;j++)
    	{
    		if(abs(B[j][i]) > abs(B[pivotrow][i]))
    			pivotrow = j;
    	}
    	if(pivotrow != i)
    	{
    		sgn = -sgn;
	    	for(k=i;k<n;k++)
	    	{
	    		t = B[i][k];
	    		B[i][k] = B[pivotrow][k];
	    		B[pivotrow][k] = t;
	    	}
	    }
    	for(j=i+1;j<n;j++)
    	{
    		t = B[j][i]/B[i][i];
    		for(k=i;k<n;k++)
    			B[j][k] -= B[i][k]*t;
    	}
    }
    for(i=0;i<n;i++) 
    	D *= B[i][i];
    D*=sgn;
    return D;
} 
  
void adjoint(float **A,float **adj,int n)
{ 
    if (n == 1) 
    { 
        adj[0][0] = 1; 
        return; 
    } 
    int sign = 1; 
    float **temp = (float **)malloc(sizeof(float *)*n); 
    for(int i=0;i<n;i++)
		temp[i] = (float *)malloc(sizeof(float)*n);
    for (int i=0; i<n; i++) 
    { 
        for (int j=0; j<n; j++) 
        { 
            getCofactor(A, temp, i, j, n); 
            sign = ((i+j)%2==0)? 1: -1;  
            adj[j][i] = (sign)*(determinant(temp,n-1));
        } 
    } 
    for(int i=0;i<n;i++)
		free(temp[i]);
    free(temp);
} 

void inverse(float **A, float **inv,int n) 
{ 
    float det = determinant(A,n),**adj = (float **)malloc(sizeof(float *)*n); 
    for(int i=0;i<n;i++)
		adj[i] = (float *)malloc(sizeof(float)*n); 
    adjoint(A,adj,n);  
    for (int i=0; i<n; i++) 
    {
        for (int j=0; j<n; j++)
            inv[i][j] = adj[i][j]/det;  
    }
    for(int i=0;i<n;i++)
		free(adj[i]);
    free(adj);
} 

void strassen(float *A,float *B,float *C, int m, int n)
{
    if(m==2)
    {
        float P=(*A+*(A+n+1))*(*B+*(B+n+1));  
        float Q=(*(A+n)+*(A+n+1))*(*B);   
        float R=(*A)*(*(B+1)-*(B+n+1));   
        float S=(*(A+n+1))*(*(B+n)-*B);   
        float T=(*A+*(A+1))*(*(B+n+1));   
        float U=(*(A+n)-*A)*(*B+*(B+1));  
        float V=(*(A+1)-*(A+n+1))*(*(B+n)+*(B+n+1));  

        *C=*C+P+S-T+V;  
        *(C+1)=*(C+1)+R+T;  
        *(C+n)=*(C+n)+Q+S;  
        *(C+n+1)=*(C+n+1)+P+R-Q+U;  
    }
    else
    {
        m=m/2;
        strassen(A,B,C,m,n);
        strassen(A,B+m,C+m,m,n);
        strassen(A+m,B+m*n,C,m,n);
        strassen(A+m,B+m*(n+1),C+m,m,n);
        strassen(A+m*n,B,C+m*n,m,n);
        strassen(A+m*n,B+m,C+m*(n+1),m,n);
        strassen(A+m*(n+1),B+m*n,C+m*n,m,n);
        strassen(A+m*(n+1),B+m*(n+1),C+m*(n+1),m,n);
    }
}

void transpose(float **A,float **T,int n)
{
     for(int i=0; i<n; ++i)
     {
        for(int j=0; j<n; ++j)
        {
            T[j][i] = A[i][j];
        }
     }
}

void two_to_one(float **A,float *b,int n)
{
	int k=0;
	for (int i=0; i<n; i++) 
    	{ 
        	for (int j=0; j<n; j++) 
    			b[k++] = A[i][j];
	}
}

void one_to_two(float *b,float **A,int n)
{
	int k=0;
	for (int i=0; i<n; i++) 
    	{ 
        	for (int j=0; j<n; j++) 
    			A[i][j] = b[k++];
	}	
}

void leastsqrfit(float **A,float **B,int n,int m)
{
	int i,j;
	float **At = (float**)malloc(sizeof(float*)*n); 
    for(int i=0;i<n;i++)
		At[i] = (float*)malloc(sizeof(float)*n);
	float *at = (float*)malloc(sizeof(float)*(n*n));
	float *a = (float*)malloc(sizeof(float)*(n*n));
	float *b = (float*)malloc(sizeof(float)*(n*n));
	float *t1 = (float*)malloc(sizeof(float)*(n*n)); 
	float *t2 = (float *)malloc(sizeof(float)*(n*n)); 
	float *t3 = (float *)malloc(sizeof(float)*(n*n)); 
	float **T1 = (float **)malloc(sizeof(float *)*n); 
    for(int i=0;i<n;i++)
		T1[i] = (float *)malloc(sizeof(float)*n);
	float **T3 = (float **)malloc(sizeof(float *)*n); 
    for(int i=0;i<n;i++)
		T3[i] = (float *)malloc(sizeof(float)*n);
	float *inv = (float *)malloc(sizeof(float)*(n*n));
	transpose(A,At,n);
	for(i=0;i<n*n;i++)
	{
        t1[i]=0;
		t2[i]=0;
		t3[i]=0;
	}
	two_to_one(At,at,n);
	two_to_one(A,a,n);
	strassen(at,a,t1,n,n);
	one_to_two(t1,T1,n);
	inverse(T1,T1,m);
	two_to_one(T1,inv,n);
	two_to_one(B,b,n);
	strassen(at,b,t2,n,n);
	strassen(inv,t2,t3,n,n);
	one_to_two(t3,T3,n);
	i=1;
	printf("y = %f",T3[0][0]);
	while(i<m)
	{
		printf(" + %f x%d",T3[i][0],i);
		i++;
	}
	printf("\n");
	for(int i=0;i<n;i++)
		free(At[i]);
    free(At);
	free(a);
	free(b);
	free(t1);
	free(t2);
	free(t3);
	for(int i=0;i<n;i++)
		free(T1[i]);
    free(T1);
	for(int i=0;i<n;i++)
		free(T3[i]);
    free(T3);
	free(inv);
}
