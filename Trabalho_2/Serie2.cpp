#include <iostream>
#include <stdlib.h>
#include "gnuplot-iostream.h"

using namespace std;

static const int N=1;

static const int N_iter=10;

static const double dt =1;


int main(int argc, char *argv[])
{
  Gnuplot gp;
  if (argc!=3)
    {
      cout << "This program works with two arguments, execute it by the following line" << endl << "[directory]/thisfile [length] [N_sheets]" << endl;
      return 0;
    }
  //Aquisição dos argumentos, comprimento da caixa e o número de folhas 
  double length=atol(argv[1]);    
  int N_sheets=atoi(argv[2]);
   
  //Posições iniciais
  double x_i=1;
  double v_i=0;
  gp << "f(x)=x**2";
  gp << "plot";

  
    
  /*
    cout << "teste" << endl;
    cout	<< "Nsheet=" << N_sheets<< endl;
    cout	<<"length="<<length<<endl ;
  */
  double *x_eq;
  double *x;
  double *v;
  double *v_eq;
  


  
  x_eq=new double[N_sheets];
  v_eq=new double[N_sheets]; 
  x=new double[N_sheets];
  v=new double[N_sheets];

  //Inicialização
  
  for (int i=0;i<N_sheets;i++)
    {
      x_eq[i]=(length/(2*N_sheets))+ (i * length/N_sheets);
      x[i]=x_i;
      v[i]=v_i;
      v_eq[i]=0;
  
    }


  for(int i=0;i<N_iter;i++)
    {
      for(int j=0;j<N_sheets;j++)
	{
	  double dx,dv;
	  dx=x[j]-x_eq[j];	   
	  dv=v[j]-v_eq[j];

	  //Força=-dx
	  
	  dv=dv-dx*dt;

	  //incremento nas velocidades

	  v[j]=dv+v_eq[j];
	  
	  //incremento nas posições
	  dx=dx+dv*dt;

	  x[j]=x_eq[j]+dx;
	    

	}

    }

  
  delete[] x_eq;
  delete[] v_eq;
  delete[] x;
  delete[] v;

  //delete[] x_i;
  // delete[] v_i;
  return 0;   
}
