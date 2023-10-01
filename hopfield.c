//the difference between this one and the previous attempt is that I'm going to use the re-calculated version of deltaH and will use a network of 0 and 1, not 1 and -1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gsl_rng.h"
#define MIN(i, j) (((i) < (j)) ? (i) : (j))
#define N 50
#define pat 1
#define T .05
#define time 20*N*N

//location -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm

//functions
void randomizeSpins(int s[][N]);
int Hopfield(int s[][N], int P[][N][N], double a[]);
double deltaH(int s[][N], int P[][N][N], int k, int l, double a[]);
double calcTheta(int k, int l, double w[][N]);
double calcA(int P[][N][N], int u);
void createPatron(FILE *p, int P[][N][N]);
double solapamiento(int s[][N], int P[][N][N], double a[], int n);

//external variable things
gsl_rng *tau;

int main()
{
	//variables
	int s[N][N];
	int P[pat][N][N];
	double a[pat];
	double sol;
	int changes = 0;
	
	//files
    	FILE *p1, *fresult, *fsol;
	p1 = fopen("yoshi.txt", "r");  
	fresult = fopen("hopfield.txt", "w");
    	fsol = fopen("solapamiento.txt", "w");
    	
	//set up network
	randomizeSpins(s);
	
	//set up patrons
	createPatron(p1, P);

	extern gsl_rng *tau;
    	int semilla = 24187353;
    	tau=gsl_rng_alloc(gsl_rng_taus);
    	gsl_rng_set(tau,semilla);

	//calculate a values - remain the same during the whole program, one value per patron
	for(int n = 0; n < pat; n++)
	{
		a[n] = calcA(P, n);
		sol = solapamiento(s, P, a, n);
		fprintf(fsol, "%lf\n", sol);
	}
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			fprintf(fresult, "%i", s[i][j]);
		}
		fprintf(fresult, "\n");
	}	
	//Run Hopfield
	for(int t = 0; t < time; t++)
	{
		//printf("\ntime = %i\n", t);
		changes = changes + Hopfield(s, P, a);

		
		//printing to data file - only after each paso Monte Carlo
		if(t % (N*N) == 0)
		{	
			for(int n = 0; n < pat; n++)
			{
				sol = solapamiento(s, P, a, n);
				fprintf(fsol, "%lf\n", sol);
			}
            	}
	}
	
	/*for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			fprintf(fresult, "%i", s[i][j]);
		}
		fprintf(fresult, "\n");
	} */
	//Hopfield ends
	printf("changes = %i\n", changes);
	return 0;
}

void randomizeSpins(int s[][N])
{
	extern gsl_rng *tau;
    	int semilla = 18237247;
    	tau=gsl_rng_alloc(gsl_rng_taus);
    	gsl_rng_set(tau,semilla);
    	
    	double x;
    	
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
		    x = gsl_rng_uniform(tau);
		    if(x > .5)
		    {
		        s[i][j] = 0;
		    }
		    else
		    {
		        s[i][j] = 1;
		    }
		}
	}
}

void createPatron(FILE *p1, int P[][N][N])
{
	int row = 0;
	while (!feof(p1))
	{
		if(ferror(p1))
		{
			printf("Error reading file.\n");
		}
		
		for (int i = 0; i < N; i++)
		{
			if(fscanf(p1, "%d", &P[0][row][i]) == EOF)
				break;
		}
		
		row++;
	}
}

double calcA(int P[][N][N], int u)
{
	double a, sum;
	
	a = (1/pow(N,2));
	sum = 0;
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			sum = sum + P[u][i][j];			
		}
	}
	return a*sum;
}

int Hopfield(int s[][N], int P[][N][N], double a[])
{
	int k, l;
	double dh,p,r;
	
	//choose random point
	k = gsl_rng_uniform_int(tau,N);
        l = gsl_rng_uniform_int(tau,N);
        
        dh = deltaH(s, P, k, l, a);
        p = MIN(1,exp(-1.*((double)dh/T)));
      
        r = gsl_rng_uniform(tau);
        
        //printf("k = %i, l = %i, s[k][l] = %i, P[0][k][l] = %i, dh = %lf, p = %lf, r = %lf\n\n", k, l, s[k][l], P[0][k][l], dh, p, r); 
        if(r < p)
        {
            s[k][l] = 1 - s[k][l];
            return 1;
        }
        return 0;
}

double deltaH(int s[][N], int P[][N][N], int k, int l, double a[])
{
	double dh, sum, theta;
	double w[N][N];
	
	dh = 1.0 - (2.0*s[k][l]);
	
	//printf("dhi = %lf\n", dh);
	sum = 0;
	
	//iterate through each point to find relation with point (k,l)
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			
			w[i][j] = 0.0;	
			
			if((k != i) || (l !=j))
			{
				for(int n = 0; n < pat; n++)
				{
					w[i][j] += ((1.0) / (N*N)) * (P[n][i][j] - a[n]) * (P[n][k][l] - a[n]);
				}
				sum = sum + w[i][j] * s[i][j];
			}
		}
	}
	
	theta = calcTheta(k, l, w); //w values are correct for the current k and l
	
	dh = dh*(theta - sum);
	return dh;
}

double calcTheta(int k, int l, double w[][N])
{
	double sum = 0;
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			if((k != i) || (l != j))
			{
				sum = sum + w[i][j];
			}
		}
	}
	sum = sum * .5;
	return sum;
}

double solapamiento(int s[][N], int P[][N][N], double a[], int n)
{
	double m, sum;
	m = (1/(N*N*a[n]*(1-a[n])));
	//printf("m = %lf\n", m);
	sum = 0;
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			sum += (P[n][i][j] - a[n])*(s[i][j] - a[n]);
		}
	}
	//printf("sum = %lf\n\n", sum);
	m = m * sum;
	return m;
}
