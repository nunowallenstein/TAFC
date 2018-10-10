#include <iostream>
#include <stdlib.h>
using namespace std;

static const int N=1;

static const int N_iter=10;

static const double dt =0.01;

int main(int argc, char *argv[])
{
  if (argc!=3)
    {
      cout << "This program works with two arguments, execute it by the following line" << endl << "./thisfile [length] [N_sheets]" << endl;
      return 0;
    }
  double length=atol(argv[1]);
    
  int N_sheets=atoi(argv[2]);
   
    
  /*
  cout << "teste" << endl;
  cout	<< "Nsheet=" << N_sheets<< endl;
  cout	<<"length="<<length<<endl ;
  */
  double x_eq[N_sheets];
  //depois pomos o método de euler para as velocidades e posições e tá a andar
    

  return 0;   
}
