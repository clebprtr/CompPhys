#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complex.h"
#define N 200
#define nc 5 //nciclos, in this case 5
#define l 1 //lambda = .3, the example in the slides
#define PI 3.141593
#define T 1000

/*To-Do
- check that the array sizes are right, note that we're starting from 0 and going to N in the distance steps
- we can calculate the normalization factor by (at each step) adding the absolute square of the value and then we apply it by dividing all function values by the sqrt of the norma
- maybe I should be multiplying by complex conjugates? - don't think this necesarily applies
*/

int main()
{

	//Files
	FILE *fdata = NULL;
	fdata = fopen("schroe.txt", "w");
	FILE *ftest = NULL;
	ftest = fopen("schroetest.txt", "w");
	FILE *fpot = NULL;
	fpot = fopen("pot.txt", "w");
	FILE *fnorm = NULL;
	fnorm = fopen("norm.txt", "w");

	//initializing variables
	double k,s, norm;
	double V[N+1]; 

	fcomplex phi[N+1], y[N], a[N], B[N], X[N+1], A[N], b[N];
	fcomplex multb, one;

	//setup
	k = (2*PI*nc)/N; //initializing k0
	s = 1/(4*pow(k,2));
	multb = Complex(0, 4/s);
	one = Complex(1, 0);

	for(int j = 0; j <= N; j++) //setting up V
	{
		if(j >= (.4*N) && j <= (.6*N))
		{
			V[j] = l*pow(k,2);
			//V[j] = .02764;
		}
		else
		{
			V[j] = 0;
		}
	}

	//setting up phi
	phi[0] = Complex(0,0);
	phi[N] = Complex(0,0);
	norm = 0;
	for(int j=1; j < N; j++) //excluding the first and last elements of the array
	{
		double factor;
		factor = exp(-8*pow((4*j - N),2)/(pow(N,2)));
		phi[j].r = cos(k*j);
		phi[j].i = sin(k*j);
		phi[j] = RCmul(factor, phi[j]);
		norm = norm + ((Cmul(phi[j],Conjg(phi[j])).r));
	}
	
	//apply the norm to the inital function:
	for (int j = 0; j < N + 1; j++)
	{
		phi[j] = RCmul(1.0/sqrt(norm), phi[j]);
	}

	//initializing a and A
	a[N-1] = Complex(0,0);
	A[N-1] = Complex(-2 - V[N-1], 2/s);
	y[N-1] = Cpow(Cadd(A[N-1], a[N-1]),-1);
	for(int j = N-1; j >=1; j--)
	{
		A[j] = Complex(-2.0 - V[j], 2.0/s); //j+1 because we're using j as the index instead of j-1 for a		
		y[j] = Cdiv(one,Cadd(A[j], a[j]));
		a[j-1] = RCmul(-1, y[j]);
		
		fprintf(ftest, "\nj = %i\n", j);
		fprintf(ftest, "V[%i] = %lf\n", j, V[j]);
		fprintf(ftest, "A[%i] = %lf + i%lf\n", j, A[j].r, A[j].i);
		fprintf(ftest, "y[%i] = %lf + i%lf\n", j, y[j].r, y[j].i);
		fprintf(ftest, "a[%i] = %lf + i%lf\n", j, a[j].r, a[j].i);
	}
	
	//printing all initialized variables
		fprintf(ftest, "\nk = %lf\n", k); //k
		fprintf(ftest, "s = %lf\n", s); //s
		for(int j = 0; j < N; j++) //V
		{
	//		fprintf(ftest, "V[%i] = %lf\n", j, V[j]);
	//		fprintf(ftest, "phi[%i] = %lf + i%lf\n", j, phi[j].r, phi[j].i);
	//		fprintf(ftest, "y[%i] = %lf + i%lf\n", j, y[j].r, y[j].i);
		}
		for(int j=0; j < N; j++)
		{
		//	fprintf(ftest, "a[%i] = %lf + i%lf\n", j, a[j].r, a[j].i);
		}
	
	//time loop
	
	for(int t=0; t < T; t++)
	{
		//print data to file
		for(int j = 0; j <= N; j++)
		{
			fprintf(fdata, "%i\t%lf\n", j, Cmul(phi[j],Conjg(phi[j])).r);
			fprintf(fpot, "%i\t%lf\n", j, V[j]);
		}
		fprintf(fnorm, "%i\t%lf\n", t, norm);
		fprintf(fdata, "\n\n");
		fprintf(fpot, "\n\n");
		fprintf(ftest, "\ntime is %i\n", t);
		
		//Calculate B, B[0] is calculated with a[0]
		for(int j = 0; j <= N-1; j++)
		{ 
		//others did this with two separate loops, calculating b at each spot and then B B to deal with the different indices and fix B(N-1) at zero.
			b[j] = Cmul(multb,phi[j]);
		}
		
		B[N-1] = Complex(0,0);
		for(int j = N-2; j >=0 ; j--)
		{
			B[j] = Cmul(y[j+1],Csub(b[j+1],B[j+1]));
		}
		
		
		//calculate X
		X[0] = Complex(0,0);
		X[N] = Complex(0,0);
		for(int j = 0; j <= N-1; j++)
		{
			X[j+1] = Cadd(Cmul(a[j],X[j]),B[j]);
			//fprintf(ftest, "X[%i] = %lf + i%lf\n", j, X[j].r, X[j].i);
		}
		
		//calculate phi and norma
		norm = 0;
		for(int j = 0; j <= N-1; j++)
		{
			phi[j] = Csub(X[j],phi[j]);
			//fprintf(ftest, "phi[%i] = %lf + i%lf\n", j, phi[j].r, phi[j].i);
			norm = norm + Cmul(phi[j],Conjg(phi[j])).r;
		}
		
		for(int j = 29; j >=0; j--)
		{
			fprintf(ftest, "\nmultb= %lf + i%lf\n", multb.r, multb.i);
			//fprintf(ftest, "b = %lf + i%lf\n", b.r, b.i);
			fprintf(ftest, "y[%i] = %lf + i%lf\n", j, y[j].r, y[j].i);
			fprintf(ftest, "B[%i] = %lf + i%lf\n", j, B[j].r, B[j].i);
			fprintf(ftest, "X[%i] = %lf + i%lf\n", j, X[j].r, X[j].i);
			fprintf(ftest, "phi[%i] = %lf + i%lf\n", j, phi[j].r, phi[j].i);
		//printing all variables to test file
		}
	}
}


