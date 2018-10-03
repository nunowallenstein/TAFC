
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>
#include <sys/time.h>



// Number of particles
static const int N = 16384;

//static const int N = 10;

// Number of iterations
static const int N_iter = 100;

// Timestep
static const double dt = 0.001;

// Softening factor
static const double soft = 1e-40;

// Timing routine based on gettimeofday
double get_wtime() {
  struct timeval tp;
  gettimeofday( &tp, NULL );
  return 1.0e-6 * tp.tv_usec + tp.tv_sec;
}


void init_rand01( double*  buffer, const int size ) {
  const double r_rand_max = 1.0/RAND_MAX;  
  for( int i = 0; i < size; i ++)
    buffer[i] = rand() * r_rand_max;
}



/*

int free2ddouble(double ***array) {
    // free the memory - the first element of the array is at the start 
    free(&((*array)[0][0]));

     free the pointers into the memory 
    free(*array);

    return 0;
}
*/



double **alloc_2d_init(int rows, int cols) {
    double *data = (double *)malloc(rows*cols*sizeof(double));
    double **array= (double **)malloc(rows*sizeof(double*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}

double kinetic_energy_local(double  **param, int N_particles, int world_rank,int n_proc)
{
  int i,my_first,my_last;
  double energy_local=0;
 
  my_first =(world_rank)*N_particles/n_proc;
  my_last=(world_rank +1)*N_particles/n_proc;
  for(i=my_first;i<my_last;i++)
  energy_local=energy_local +param[3][i]*param[3][i]+param[4][i]*param[4][i]+param[5][i]*param[5][i];

  return 0.5*energy_local;
}
 
    

int main(int argc, char *argv[])
{

  int world_rank,world_size,err;
  double **param;
  double energy,energy_local=0;
  
   
  param=alloc_2d_init(6,N);
 err= MPI_Init(&argc,&argv);
 err= MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
 err= MPI_Comm_size(MPI_COMM_WORLD,&world_size);


  if(world_rank==0)
    {
      printf("Initializing %d particles...\n", N );
      srand(0);
       for(int i=0;i<6;i++)
       init_rand01(param[i],N);	
      printf("Initialization complete.\n");
    }
 MPI_Bcast(&(param[0][0]), 6*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
 energy_local=kinetic_energy_local(param,N,world_rank,world_size);
  
 err=MPI_Allreduce(&energy_local,&energy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
 printf("energy %f\n",energy);

 free(param[0]);
 free(param);
  MPI_Finalize();

  return 0;
}
