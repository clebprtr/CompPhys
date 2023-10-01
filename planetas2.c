#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#define MS 1.9891e30 //Mass of the Sun - Kg
#define C 1.496e11 //Distance from Earth to Sun - m
#define H .05 //step size
#define G 6.67e-11 // Nm^2/Kg^2
#define K sqrt((G*MS)/(pow(C,3)))

//Functions
void convertPosition(double r[]);
void convertMass(double m[]);
void convertVelocity(double v[]);
void getAccelerations(double ax[], double ay[], double m[], double x[], double y[]);
void getPosition(double r[], double v[], double a[]);
void getW(double w[], double v[], double a[]);
void getVelocity(double v[], double w[], double a[]);
void checkPosition(double xi[], double x[], double yi[], double y[], double p[], int t);
void getEnergy(double e[], double x[], double y[], double m[]);

//Main
int main()
{
    //Files
    FILE *fdata = NULL;
    FILE *ftest = NULL;
    ftest = fopen("planets.txt", "w");
    fdata = fopen("planets_data.txt", "w");    
    //Variables
    double y[9],yi[9],vx[9],ax[9],ay[9],wx[9],wy[9],p[9];    

    double x[9] = {57.9, 108.2, 149.6, 227.9, 778.6, 1433.5, 2872.5, 4495.1, 0}; //distance from sun en 10^6 km
    double xi[9] = {57.9, 108.2, 149.6, 227.9, 778.6, 1433.5, 2872.5, 4495.1, 0}; //distance from sun en 10^6 km
    double m[9] = {.330, 4.87, 5.97, .642, 1899, 568, 86.8, 102, 1988500}; //mass (x10^24) kg
    double vy[9] = {47.9, 35.0, 29.8, 24.1, 13.1, 9.7, 6.8, 5.4, 0}; //orbital velocity km/s
    
   //Variable Conversions
    convertPosition(x);
    convertPosition(xi);
    convertMass(m);
    convertVelocity(vy);

    getAccelerations(ax,ay,m,x,y);

    for(int t=0; t<150; t++)
    {
        getPosition(x, vx, ax);
        getPosition(y, vy, ay);

        checkPosition(xi,x,yi,y,p,t);
        for(int i=0; i < 9; i++)
        {        
            double div;
            div = x[i] / xi[i]; //check if they're the same
            fprintf(ftest,"%i: div %lf at time %i\n", i, div, t);
            fprintf(ftest,"%i: %lf at time %i\n", i, p[i], t);
        }
        fprintf(ftest, "\n");

        for(int i=0;i<9;i++)
        {
         //   if(t % 10 == 0)
           // {            
                fprintf(fdata,"%lf\t%lf\t", x[i], y[i]);
            //}        
        }

        //    if(t % 10 == 0)
          //  {            
                fprintf(fdata,"\n");
           // }     

        getW(wx, vx, ax);
        getW(wy, vy, ay);

        getAccelerations(ax,ay,m,x,y);

        getVelocity(vx, wx, ax);
        getVelocity(vy, wy, ay);
    }
    for(int i = 0; i < 9; i++)
    {
        printf("%lf\n", p[i]);
    }
    return 0;
}

void getVelocity(double v[], double w[], double a[])
{
    int i;
    for(i=0; i<9; i++)
    {
        v[i] = w[i] + (H/2)*a[i];
    }
}

void getW(double w[], double v[], double a[])
{
    int i;
    for(i=0; i<9; i++)
    {
        w[i] = v[i] + (H/2)*a[i];
    }
}

void getPosition(double r[], double v[], double a[])
{
    int i;    
    for(i = 0; i<9; i++)
    {
        r[i] = r[i] + H*v[i] + (pow(H,2)/2)*a[i];
    }
}

void checkPosition(double xi[], double x[], double yi[], double y[], double p[], int t)
{
    int i;
    double div;
    for(i = 0; i < 9; i++)
    {
        if(p[i] == 0 && t>28) //if we don't already have a period and the time isn't zero
        {
            div = x[i] / xi[i]; //check if they're the same
            
            if(0.99 < div && div < 1.01)
            {
                p[i] = t * 1.; //convert t to double to put in array
            }  
            else
            {
                p[i] = 0;
            }
        }
    }
}
//idfk

void getEnergy(double e[], double x[], double y[], double m[])
{
    double r;    
    for(int i = 0; i < 9; i++)
    {
        r = sqrt((pow(x[i],2)) + (pow(y[i],2)));
        e[i] = (-1*G*m[i])/(2*r);
    }
}

void getAccelerations(double ax[], double ay[], double m[], double x[], double y[])
{

    int p,j;
    double sumx, sumy, mod, addingx, addingy;
  
    for(p=0; p<9; p++)
    {
        sumx = 0; //reset sum
        sumy = 0;        
        for(j=0; j<9; j++)
        {
            if(p==j)
            {
                sumx = sumx;
                sumy = sumy;
            }
            else
            {
               mod = 0;
               addingx = 0;
               addingy = 0;
               mod = sqrt((pow((x[p]-x[j]),2) + (pow((y[p]-y[j]),2))));
               addingx = m[j]*(x[p] - x[j])/(pow(mod,3));
               addingy = m[j]*(y[p] - y[j])/(pow(mod,3));
               sumx = sumx + addingx;
               sumy = sumy + addingy;
            }
        }
        ax[p] = -1.0*sumx;
        ay[p] = -1.0*sumy;
    }
}

void convertPosition(double r[])
{
    int i;
    for(i = 0; i < 9; i++)
    {
        r[i] = r[i]*pow(10,9);     
        r[i] = r[i]/C;
    }
}

void convertMass(double m[])
{
    int i;
    for(i = 0; i < 9; i++) 
    {
        m[i] = m[i]*pow(10,24); //convert to just kg
        m[i] = m[i]/MS; //convert to our variables
    }
}

void convertVelocity(double v[])
{
    int i;
    for(i=0; i<9; i++)
    {
        v[i] = v[i]*1000; //convert km/s to m/s        
        v[i] = v[i]/C;
        v[i] = v[i]/K;        
    }
}
