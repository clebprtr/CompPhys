#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 6.67e-11
#define MT 5.9736e24
#define ML 0.07349e24
#define dTL 3.844e8
#define w 2.6617e-6 //s^-1
#define RT 6.378160e6
#define RL 1.7374e6
#define del (G*MT)/pow(dTL,3)
#define u ML/MT
#define T 604800
#define h 60
#define pi 3.1415927

/* Notes:
- Earth will just be represented by printing 0,0 every time step
- arrays of k for each variable...
- .007825
- .005511
- .003237
- .001028 -.108
- .001180 -.109
- .000113 -.1085
- .000107 - .10848
- .000108 - .10849
- apollo mission orbit was at 100 - 122 km above

*/

//functions
double rdot(double pr);
double phidot(double r, double pphi);
double prdot(double r, double phi, double pr, double pphi, double t);
double pphidot(double r, double phi, double t);
double calcDist(double xR, double yR, double t);

int main()
{
	//files
	FILE *fdata = NULL;
	fdata = fopen("cohete.txt", "w");
	FILE *ftest = NULL;
	ftest = fopen("cohetetest.txt", "w");

	//variables
	int count;
	double r, phi, pr, pphi, xL, yL, xR, yR, zero, vi, theta, dist, distmin;
	double kr[4], kphi[4], kpr[4], kpphi[4];
	
	//give initial conditions of r, phi, pr, pphi - problem 2
	theta = 0;
	vi = 0; //1600m/s rescaled
	r = 1 - (RL/dTL) - (110000/dTL); //starts 110km to the right of the moon's surface
	phi = 0; //starts in line with moon and earth
	pr = 0; //if in orbit and changing direction (w.r.t. radius from earth, should be zero
	pphi = 0; //equal to the moon
	distmin = 10;
	
	/*for initial problem
	//give initial conditions of r, phi, pr, pphi
	theta = 0;
	vi = 1600/dTL; //1600m/s rescaled
	r = 1 - (110000/dTL); //starts 110km to the right of the moon
	phi = 0;
	pr = vi*cos(phi-theta); //m/s then rescaled
	pphi = r*vi*sin(phi-theta); //.0005 r/s then rescaled
	distmin = 10;
	*/	
	
	//time loop
	for(double t = 0; t < T; t = t + h)
	{
		//k1
		kr[0] = h*rdot(pr);
		kphi[0] = h*phidot(r, pphi);
		kpr[0] = h*prdot(r, phi, pr, pphi,t);
		kpphi[0] = h*pphidot(r, phi, t);
		
		//k2-k3
		for(int i = 1; i < 3; i++)
		{
			kr[i] = h*rdot(pr + kpr[i-1]/2);
			kphi[i] = h*phidot(r + kr[i-1]/2, pphi + kpphi[i-1]/2);
			kpr[i] = h*prdot(r + kr[i-1]/2, phi + kphi[i-1]/2, pr + kpr[i-1]/2, pphi + kpphi[i-1]/2,t + h/2);
			kpphi[i] = h*pphidot(r + kr[i-1]/2, phi + kphi[i-1]/2, t + h/2);
		}
		
		//k4
		kr[3] = h*rdot(pr + kpr[2]);
		kphi[3] = h*phidot(r + kr[2], pphi + kpphi[2]);
		kpr[3] = h*prdot(r + kr[2], phi + kphi[2], pr + kpr[2], pphi + kpphi[2],t + h);
		kpphi[3] = h*pphidot(r + kr[2], phi + kphi[2], t + h);
	
		//update values	
		r = r + (1.0/6.0)*(kr[0] + 2*kr[1] + 2*kr[2] + kr[3]);
		phi = phi + (1.0/6.0)*(kphi[0] + 2*kphi[1] + 2*kphi[2] + kphi[3]);
		pr = pr + (1.0/6.0)*(kpr[0] + 2*kpr[1] + 2*kpr[2] + kpr[3]);
		pphi = pphi + (1.0/6.0)*(kpphi[0] + 2*kpphi[1] + 2*kpphi[2] + kpphi[3]);
		
		//Luna
		xL = cos(w*t);
		yL = sin(w*t);
		
		//get x and y for cohete
		xR = r*cos(phi);
		yR = r*sin(phi);
		
		//calculate distance between moon surface (point closest to earth) and rocket
		dist = calcDist(xR, yR, t);
		
		//test printing
		fprintf(ftest, "\nt = %lf\n", t);
		/*for(int i = 0; i < 4; i++)
		{
			fprintf(ftest, "kr[%i] = %lf\n", i, kr[i]);	
			fprintf(ftest, "kphi[%i] = %lf\n", i, kphi[i]);	
		}
		fprintf(ftest, "r = %lf\n", r);
		fprintf(ftest, "phi = %lf\n", phi);
		fprintf(ftest, "pr = %lf\n", pr);
		fprintf(ftest, "pphi = %lf\n", pphi);*/
		fprintf(ftest, "dist = %lf\n", dist);
		
		//data printing
		if(count % 10 == 0) //should give 1008 frames
		{
			fprintf(fdata, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", xR, yR, xL, yL, zero, zero);
		}
		count++;
		
		//find smallest distance
		if(dist < distmin)
		{
			distmin = dist;
		}
	}
	
	printf("\ndistmin = %lf\n", distmin);
	printf("distmin(m) = %lf\n", distmin*dTL);
	return 0;
}

double rdot(double pr)
{
	return pr;
}

double phidot(double r, double pphi)
{
	return pphi/(pow(r,2));
}

double prdot(double r, double phi, double pr, double pphi, double t)
{
	double prdot, rprime;
	rprime = sqrt(1 + pow(r,2) - 2*r*cos(phi - w*t));
	prdot = pow(pphi,2)/(pow(r,3)) - del*((1/pow(r, 2)) + (u/pow(rprime,3))*(r - cos(phi - w*t)));
	return prdot;
}

double pphidot(double r, double phi, double t)
{
	double pphidot, rprime;
	rprime = sqrt(1 + pow(r,2) - 2*r*cos(phi - w*t));
	pphidot = ((-1.0*del*u*r)/(pow(rprime,3)))*sin(phi - w*t);
	return pphidot;
}

double calcDist(double xR, double yR, double t)
{
	double xLs, yLs, dist;
	xLs = (1-RL/dTL)*cos(w*t);
	yLs = (1-RL/dTL)*sin(w*t);
	dist = sqrt(pow((xR - xLs),2) + pow((yR - yLs),2));
	return dist;
}
