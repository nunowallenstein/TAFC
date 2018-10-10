#include <iostream>
#include <stdlib.h>
using namespace std;

static const int N=1;

static const int N_iter=10;

static const double dt =1;


int main(int argc, char *argv[])
{
  if (argc!=3)
    {
      cout << "This program works with two arguments, execute it by the following line" << endl << "./thisfile [length] [N_sheets]" << endl;
      return 0;
    }
  //Aquisição dos argumentos, comprimento da caixa e o número de folhas 
  double length=atol(argv[1]);    
  int N_sheets=atoi(argv[2]);
   
  //Posições iniciais
  double x_i=1;
  double v_i=0;

  
    
  /*
  cout << "teste" << endl;
  cout	<< "Nsheet=" << N_sheets<< endl;
  cout	<<"length="<<length<<endl ;
  */
  double *x_eq;
  double *x;
  double *v;

  
  double *Dlt_x;
  double *Dlt_v;

  
  x_eq=new double[N_sheets];
  x=new double[N_sheets];
  v=new double[N_sheets];
  Dlt_x=new double[N_sheets];
  Dlt_v=new double[N_sheets];
  for (int i=0;i<N_sheets;i++)
    {
      x_eq[i]=(length/(2*N_sheets))+ (i * length/N_sheets);
      cout << x_eq[i] << endl;
      x[i]=x_i;
      v[i]=v_i;
    }


  for(int i=0;i<N_iter;i++)
    {

    }

  
  delete[] x_eq;
  //delete[] x_i;
  // delete[] v_i;
  return 0;   
}
