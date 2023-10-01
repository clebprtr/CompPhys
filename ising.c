#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gsl_rng.h"
#define MIN(i, j) (((i) < (j)) ? (i) : (j))
#define N 32
//location -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm
/*    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
*/

//external variable things
gsl_rng *tau;

//Functions
//double magnetizacion(int s[][N]);

int main()
{
    //generating random numbers
    extern gsl_rng *tau;
    int semilla = 18237247;
    tau=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau,semilla);

    //files
    FILE *fdata = NULL;
    fdata = fopen("ising.txt", "w");
    FILE *ftest = NULL;
    ftest = fopen("isingtest.txt", "w");

    //Variables
    double T = 3;
    double x,p,r,mag;
    int n, m, e, n_plus_1, n_minus_1, m_plus_1, m_minus_1, counter;
    int s[N][N];
    
    //set up spins    
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            x = gsl_rng_uniform(tau);
            if(x > .5)
            {
                s[i][j] = -1;
            }
            else
            {
                s[i][j] = 1;
            }
        }
    }
    //change Spins
    for (int i = 0; i < 30000; i++)
    {
        
        //choose random point
        n = gsl_rng_uniform_int(tau,N);
        m = gsl_rng_uniform_int(tau,N);

        //get e and p
        n_plus_1 = n + 1;
        n_minus_1 = n - 1;
        m_plus_1 = m + 1;
        m_minus_1 = m - 1;
        
        if(n_plus_1 == N){n_plus_1 = 0;}
        if(n_minus_1 == -1) {n_minus_1 = N;}
        if(m_plus_1 == N){m_plus_1 = 0;}
        if(m_minus_1 == -1) {m_minus_1 = N;}

        e = 2*(s[n][m])*(s[n_plus_1][m] + s[n_minus_1][m] + s[n][m_plus_1] + s[n][m_minus_1]);
        p = MIN(1,exp(-1.*((double)e/T)));
        r = gsl_rng_uniform(tau);
        
        if(r < p)
        {
            s[n][m] = -1* s[n][m];
        }

        //Print spins into document with each paso Monte Carlo
        //Note: for printing we want to have x (tab) y (tab) color (tab), for the color we want 1 or 2
        if(i % (N*N) == 0)
        {
            for(int i = 0; i < N; i++)
            {
                for(int j = 0; j < N; j++)
                {
                    fprintf(fdata, "%i\t%i\t%i\t", i, j, s[i][j]+2);
                }
            }   
            fprintf(fdata, "\n"); 
        }
        /*if(i % (N*N*100) == 0)
        {
            mag = mag + magnetizacion(s);
            fprintf(ftest, "mag at step %i is %lf\n", i, mag);
            counter = counter + 1;
        }*/
    }    //end for loop for changing spins
    //mag = mag / counter;
    //printf("mag at temp %lf = %lf\n", T, mag);    
    return 0;
}
/*
double magnetizacion(int s[][N])
{
    int sum;    
    double m;
    m = 1/pow(N,2);
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {    
            sum = sum + s[i][j];  
        }
    }
    m = m * abs(sum);
    return m;
} */

