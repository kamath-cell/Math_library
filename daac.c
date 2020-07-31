#include"daahead.h"

int main()
{
	int choice,i,j,loop=1,n,m;
	struct timespec start,end;
	double time;
	printf("1.Mean\n2.Median\n3.Mode\n4.Standard Deviation\n5.Harmonic mean\n6.Least Square Fit\n\n");
	while(loop)
	{
		printf("Enter the choice\n");
		scanf("%d",&choice);
		switch(choice)
		{
			case 1:printf("Enter the no. of elements\n");
				scanf("%d",&n);
				double *A = (double *)malloc(sizeof(double)*n);
				printf("Enter the elements\n");
				for(i=0;i<n;i++)
					scanf("%lf",&A[i]);
			    clock_gettime(CLOCK_REALTIME,&start);
				printf("Mean = %lf\n\n",mean(A,n));
				clock_gettime(CLOCK_REALTIME,&end);
				time = end.tv_nsec - start.tv_nsec;
				printf("Time taken = %lf nsec\n\n",time); 
				break;
			case 2:printf("Enter the no. of elements\n");
				scanf("%d",&n);
				double *B = (double *)malloc(sizeof(double)*n);
				printf("Enter the elements\n");
				for(i=0;i<n;i++)
					scanf("%lf",&B[i]);
				clock_gettime(CLOCK_REALTIME,&start);
				printf("Median = %lf\n\n",median(B,n));
				clock_gettime(CLOCK_REALTIME,&end);
				time = end.tv_nsec - start.tv_nsec;
				printf("Time taken = %lf nsec\n\n",time); 
				break;
			case 3:printf("Enter the no. of elements\n");
				scanf("%d",&n);
				double *C = (double *)malloc(sizeof(double)*n);
				printf("Enter the elements\n");
				for(i=0;i<n;i++)
					scanf("%lf",&C[i]);
				clock_gettime(CLOCK_REALTIME,&start);
				printf("Mode = %lf\n\n",mode(C,n));
				clock_gettime(CLOCK_REALTIME,&end);
				time = end.tv_nsec - start.tv_nsec;
				printf("Time taken = %lf nsec\n\n",time);
				break;
			case 4:printf("Enter the no. of elements\n");
				scanf("%d",&n);
				double *D = (double *)malloc(sizeof(double)*n);
				printf("Enter the elements\n");
				for(i=0;i<n;i++)
					scanf("%lf",&D[i]);
				clock_gettime(CLOCK_REALTIME,&start);
				printf("Standard Deviation = %lf\n\n",stddev(D,n));
				clock_gettime(CLOCK_REALTIME,&end);
				time = end.tv_nsec - start.tv_nsec;
				printf("Time taken = %lf nsec\n\n",time);
				break;
			case 5:printf("Enter the no. of elements\n");
				scanf("%d",&n);
				double *E = (double *)malloc(sizeof(double)*n);
				printf("Enter the elements\n");
				for(i=0;i<n;i++)
					scanf("%lf",&E[i]);
				clock_gettime(CLOCK_REALTIME,&start);
				printf("Harmonic mean of numbers = %lf\n\n",hm(E,n));
				clock_gettime(CLOCK_REALTIME,&end);
				time = end.tv_nsec - start.tv_nsec;
				printf("Time taken = %lf nsec\n\n",time);
				break;
			case 6:printf("Enter the appropriate order of matrix A\n");
				scanf("%d",&n);
				float **a = (float **)malloc(sizeof(float *)*n); 
    			for(int i=0;i<n;i++)
					a[i] = (float *)malloc(sizeof(float)*n);
				float **b = (float **)malloc(sizeof(float *)*n); 
    			for(int i=0;i<n;i++)
					b[i] = (float *)malloc(sizeof(float)*n);
				printf("Enter the elements of A\n");
				for(i=0;i<n;i++)
				{
					for(j=0;j<n;j++)
						scanf("%f",&a[i][j]);
				}
				printf("Enter no. of columns of A\n");
				scanf("%d",&m);
				printf("Enter the elements of b\n");
				for(i=0;i<n;i++)
				{
					for(j=0;j<n;j++)
						scanf("%f",&b[i][j]);
				}
				printf("The least square solution for the given problem\n");
				clock_gettime(CLOCK_REALTIME,&start);
				leastsqrfit(a,b,n,m);
				clock_gettime(CLOCK_REALTIME,&end);
				time = end.tv_nsec - start.tv_nsec;
				printf("Time taken = %lf nsec\n\n",time);
				break;
			case 7:loop = 0;
		}
	}
	return 0;
}
