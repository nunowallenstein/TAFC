
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
/*
int malloc2ddouble(double ***array, int n, int m) {

     allocate the n*m contiguous items 
    double *p = (double *)malloc(n*m*sizeof(double));
    if (!p) return -1;

     allocate the row pointers into the memory 
    (*array) = (double **)malloc(n*sizeof(double*));
    if (!(*array)) {
       free(p);
       return -1;
    }
}
*/
int main(int argc, char *argv[])
{

  int world_rank,world_size;
  double **param;
  double y=2;
 

   
  param=alloc_2d_init(6,N);
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);


  if(world_rank==0)
    {
      printf("Initializing %d particles...\n", N );
      srand(0);
       for(int i=0;i<6;i++)
       init_rand01(param[i],N);	
      printf("Initialization complete.\n");
    }


    
  MPI_Bcast(&(param[0][0]), 6*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if(world_rank==2)
    {
  printf("3-Processo %d\n",world_rank);
  printf("%f\n",param[0][2]);
  printf("%f\n",param[2][10000]);
  printf("%f\n",param[4][500]);
 printf("---------------------");
    }

     /* free(param[0]);
  free(param);
     */
   MPI_Finalize();

  return 0;
}
