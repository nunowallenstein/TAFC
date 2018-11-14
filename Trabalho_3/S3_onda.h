//======================================
//Guard
#ifndef __S3_ONDA_H_INCLUDED
#define __S3_ONDA_H_INCLUDED

//============================

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>

//Criação de um eixo uniformemente distribuido segundo dois limites
using namespace std;
void create_grid(double *&y,int Npoints,double x_inf, double x_sup)
{
  y=new double[Npoints];
  //Discretização do espaço
  for(int i=0;i<Npoints;i++)
    {
      y[i]=x_inf +i*(x_sup-x_inf)/(Npoints-1);
    }
}





//Inicialização de Psi
void psi0( double *x,double *&gauss,int Npointsx, double sigma, double x_0)//criação de psi0
{
  gauss=new double[Npointsx];
  for(int i=0;i<Npointsx;i++)
    {
    
      gauss[i]=exp(-pow((x[i]-x_0),2)/(2*pow(sigma,2)));

    }

}


//Diferenças finitas de 2ª ordem multiplicada por uma constante C, tem duas opções:
//-"none" derivar sem condições fronteira
//-"periodicas" derivar com condições fronteira periódicas
double *diffs_2nd_order(double *f,int Npoints,double dx,double C,string option)
{
  double *diff=new double[Npoints];
  
  if(option=="none")
    {
      diff[0]=C*(-3.0*f[0]+4.0*f[1]-1.0*f[2])/(2.0*dx);  
      for(int i=1;i<Npoints-1;i++)
	diff[i]=C*(1.0*f[i+1]-1.0*f[i-1])/(2.0*dx);
    
      diff[Npoints-1]=C*(3.0*f[Npoints-1]-4.0*f[Npoints-2]+1.0*f[Npoints-3])/(2.0*dx);
    }

  if(option=="periodicas")
    {
      diff[0]=C*(f[1]-f[Npoints-1])/(2.0*dx);     
      for(int i=1;i<Npoints-1;i++)
	diff[i]=C*(f[i+1]-f[i-1])/(2.0*dx);
    
      diff[Npoints-1]=C*(f[0]-f[Npoints-2])/(2.0*dx);
      
    }

  return diff;
  delete diff;
}

//Diferenças finitas de 2ª ordem para o caso de condições fronteira reflectoras
double *diffs_reflect(double *f, double *g,int Npoints, double dx, double C,double Kreflect_0,double Kreflect_1)
//para o RHS de phi, f tem que ser pi, para o RHS de pi f tem que ser phi
{
  double *diff;
  diff=new double [Npoints];
  diff[0]=C*(1.0/2)*(1+Kreflect_0)*(-3.0*(f[0]+g[0])+4.0*(f[1]+g[1])-1.0*(f[2]+g[2]))/(2.0*dx);
  for(int i=1;i<Npoints-1;i++)
    {
      diff[i]=C*(f[i+1]-f[i-1])/(2.0*dx);
    } 
  diff[Npoints-1]=C*(1.0/2)*(1+Kreflect_1)*(3.0*(f[Npoints-1]-g[Npoints-1])-4.0*(f[Npoints-2]-g[Npoints-2])+1.0*(f[Npoints-3]-g[Npoints-3]))/(2.0*dx);
  return diff;
  delete diff;
}
  


//soma de dois arrays
double *array_sum(double *f, double *g,int Npoints)
{
  double *sum=new double[Npoints];
  for (int i=0;i<Npoints;i++)
    {
      sum[i]=f[i]+g[i];

    }
  return sum;
  delete sum;
}


//multiplicação de um array por uma constante
double *array_multiply(double *f,double C, int Npoints)
{
  double *g=new double[Npoints];
  for(int i=0;i<Npoints;i++)
    {
      g[i]=C*f[i];
    }
  return g;
  delete g;
}

//Runge kutta

void runge_kutta(double *phi,double *pi,double *psi,int N_x, double dt, double dx, string Bconditions,double kR0,double kR1)
{
  double **k1=new double *[3];
  double **k2=new double *[3];
  double **k3=new double *[3];
  double **k4=new double *[3];

 
  for(int i=0;i<3;i++)
    {
      k1[i]=new double [N_x];
      k2[i]=new double [N_x];
      k3[i]=new double [N_x];
      k4[i]=new double [N_x];
    }

  
  if(Bconditions!="reflectoras")
    {
      //phi
      k1[0]=diffs_2nd_order(pi,N_x,dx,dt,Bconditions);
    
      k2[0]=diffs_2nd_order(array_sum(pi,array_multiply(k1[0],0.5,N_x),N_x),N_x,dx,dt,Bconditions);
     
      k3[0]=diffs_2nd_order(array_sum(pi,array_multiply(k2[0],0.5,N_x),N_x),N_x,dx,dt,Bconditions);
      
      k4[0]=diffs_2nd_order(array_sum(pi,k3[0],N_x),N_x,dx,dt,Bconditions);
    
      
      //pi
      k1[1]=diffs_2nd_order(phi,N_x,dx,dt,Bconditions);

      
      k2[1]=diffs_2nd_order(array_sum(phi,array_multiply(k1[1],0.5,N_x),N_x),N_x,dx,dt,Bconditions);
   
      k3[1]=diffs_2nd_order(array_sum(phi,array_multiply(k2[1],0.5,N_x),N_x),N_x,dx,dt,Bconditions);
   
      k4[1]=diffs_2nd_order(array_sum(phi,k3[1],N_x),N_x,dx,dt,Bconditions);
     
    }


    if(Bconditions=="reflectoras")
    {
      //phi
   
      k1[0]=diffs_reflect(pi,phi,N_x,dx,dt,kR0,kR1);
      k2[0]=diffs_reflect(array_sum(pi,array_multiply(k1[0],0.5,N_x),N_x),array_sum(phi,array_multiply(k1[0],0.5,N_x),N_x),N_x,dx,dt,kR0,kR1);
      k3[0]=diffs_reflect(array_sum(pi,array_multiply(k2[0],0.5,N_x),N_x),array_sum(phi,array_multiply(k2[0],0.5,N_x),N_x),N_x,dx,dt,kR0,kR1);
      k4[0]=diffs_reflect(array_sum(pi,k3[0],N_x),array_sum(phi,k3[0],N_x),N_x,dx,dt,kR0,kR1);

      //pi      
      k1[1]=diffs_reflect(phi,pi,N_x,dx,dt,-kR0,-kR1);   
      k2[1]=diffs_reflect(array_sum(phi,array_multiply(k1[1],0.5,N_x),N_x),array_sum(pi,array_multiply(k1[1],0.5,N_x),N_x),N_x,dx,dt,-kR0,-kR1);
      k3[1]=diffs_reflect(array_sum(phi,array_multiply(k2[1],0.5,N_x),N_x),array_sum(pi,array_multiply(k2[1],0.5,N_x),N_x),N_x,dx,dt,-kR0,-kR1);
      k4[1]=diffs_reflect(array_sum(phi,k3[1],N_x),array_sum(pi,k3[1],N_x),N_x,dx,dt,-kR0,-kR1);

      
    }

  //psi
  for(int i=0;i<N_x;i++)
    {
      k1[2][i]=pi[i]*dt;
      k2[2][i]=(pi[i]+k1[2][i]*0.5)*dt;
      k3[2][i]=(pi[i]+k2[2][i]*0.5)*dt;
      k4[2][i]=(pi[i]+k3[2][i])*dt;
      
    }
    
  
									
  for(int i=0;i<N_x;i++)
    {
      
      phi[i]=phi[i]+(1.0/6)*(k1[0][i]+2*k2[0][i]+2*k3[0][i]+k4[0][i]);
      pi[i]=pi[i]+(1.0/6)*(k1[1][i]+2*k2[1][i]+2*k3[1][i]+k4[1][i]);
      psi[i]=psi[i]+(1.0/6)*(k1[2][i]+2*k2[2][i]+2*k3[2][i]+k4[2][i]);
    
    }


  for(int i=0;i<3;i++)
    {
      delete k1[i];
      delete k2[i];
      delete k3[i];
      delete k4[i];
    }
  delete k1;
  delete k2;
  delete k3;
  delete k4;

			    
}


//Escrita para ficheiros
void write_file(ofstream& file,double *x,int N)
{
  for(int i=0 ;i<N;i++)
    file<<x[i]<<" ";

  file<<endl;   
}
#endif //__S3_ONDA_H_INCLUDED__
