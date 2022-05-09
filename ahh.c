#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
int i,j;
double x,t;
#define bplus 101
#define meshout 102
#define meshsize 100
#define boundary 30
#define gamma 1.4
#define tstep 0.0005
#define xstep 0.01
#define CFL 0.5

double **genq(double *rho, double *mom,double *velocity,double *E, double *energy)
{
    double *pressure = malloc(meshout*sizeof(double));
    double **Q = malloc(3*sizeof(double));
    for(i=0;i<=2;i++)
    {
        Q[i] = malloc(meshout*sizeof(double));
    }
    for(i=1;i<=meshsize;i++)
    {
        pressure[i] = ((gamma-1)*rho[i]*energy[i]);
        Q[0][i] = rho[i];
        Q[1][i] = mom[i];
        Q[2][i] = E[i];//e not epsilon
    }
    for(j=0;j<=2;j++)
    {
        Q[j][0]=Q[j][1];
        Q[j][bplus]=Q[j][meshsize];
    }
    return Q;
}

double **genf(double *rho, double *mom,double *velocity,double *E, double *energy)
{
    double *pressure = malloc(meshout*sizeof(double));
    double **F = malloc(3*sizeof(double));
    for(i=0;i<=2;i++)
    {
        F[i] = malloc(meshout*sizeof(double));
    }
    for(i=1;i<=meshsize;i++)
    {
        pressure[i] = ((gamma-1)*rho[i]*energy[i]); 
        velocity[i]=mom[i]/rho[i];
        F[0][i]=mom[i];
        F[1][i] = ((rho[i]*velocity[i]*velocity[i])+pressure[i]);
        F[2][i] = (velocity[i]*E[i])+(velocity[i]+pressure[i]);
    }

    for(j=0;j<=2;j++)
    {
        F[j][0]=F[j][1];
        F[j][bplus]=F[j][meshsize];
    }
    return F;
}

int main()
{
    FILE *ptr;
    ptr = fopen("q1a.txt","w");
    //initialise state arrays
    double *Energy = malloc(meshout*sizeof(double));
    double *energy =  malloc(meshout*sizeof(double));
    double *pressure =  malloc(meshout*sizeof(double));
    double *rho =  malloc(meshout*sizeof(double));
    double *mom =  malloc(meshout*sizeof(double));
    double *velocity =  malloc(meshout*sizeof(double)); 

    //initialise state vectors 
    //load initial conditions into state arrays 
    for(i=1;i<=meshsize;i++)
    {
        if(i<=boundary)
        {
            pressure[i] = 1;
            rho[i] = 1;
            velocity[i] = 0.75;
        }
        else if(i>boundary && i<=meshsize)
        {
            pressure[i] = 0.1;
            rho[i] = 0.125;
            velocity[i] = 0;
        }
        Energy[i]=((pressure[i])/(gamma-1)*rho[i])+(0.5*(rho[i]*(velocity[i]*velocity[i])));
        energy[i]=(Energy[i]/rho[i])-(0.5*(velocity[i]*velocity[i]));
        mom[i]= rho[i]*velocity[i];
    }
    double **Q = genq(rho,mom, velocity,Energy,energy);
    double **qnext = genq(rho,mom, velocity,Energy,energy);
    
    double **F = genq(rho, mom, velocity,Energy, energy);
    
    //combine state arrays into vectors
    
    //start iterative timeloop
    for(t=0;t<=0.2;t=t+tstep)
    {
        for(i=1;i<=meshsize;i++)
        {
            for(j=0;j<=2;j++)
            {
                int kp = i+1;
                int km = i-1;
                double fminus = (double) 0.5*(F[j][km]+F[j][i])- (((xstep)/(2*tstep))*(Q[j][i]-Q[j][km]));
                double fplus = (double) 0.5*(F[j][i]+F[j][kp])+ (((xstep)/(2*tstep))*(Q[j][i]-Q[j][kp]));
                qnext[j][i] = Q[j][i] - (((tstep)/(xstep))*(fplus-fminus));
                //qnext[j][i] = (0.5*(Q[j][i+1]+Q[j][i-1]))-((tstep/(2*xstep))*(F[j][i+1]-F[j][i-1]));
                
            }
        }
        for(i=1;i<=meshsize;i++)
        {
            rho[i]=qnext[0][i];
            mom[i]=qnext[1][i];
            velocity[i] = mom[i]/rho[i];
            Energy[i] = qnext[2][i];
            energy[i]=(Energy[i]/rho[i])-(0.5*(velocity[i]*velocity[i]));
            
        }
        Q = genq(rho, mom, velocity, Energy, energy);
        F = genf(rho,mom, velocity, Energy, energy);
        //fill new updated state vectors 
        
                
        //new state calculated above -> load into new vectors 
    }
    x=0;
    for(i=0;i<=bplus;i++)
    {
        x=x+xstep;
        rho[i] = Q[0][i];
        velocity[i] = Q[1][i]/rho[i];
        Energy[i] = Q[2][i];
        energy[i]=(Energy[i]/rho[i])-(0.5*(velocity[i]*velocity[i]));
        pressure[i] = ((gamma-1)*rho[i]*energy[i]); 
        fprintf(ptr,"%f,%f,%f,%f,%f\n", rho[i], velocity[i], energy[i], pressure[i],x);
    }
    return 0;
}
