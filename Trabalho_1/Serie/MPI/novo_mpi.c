//Este programa está organizado da seguinte forma, substituímos a struct das partículas por um array 2d de doubles,param, em que os o primeiro indice do array refere-se à propriedade a estudar e o segundo a partícula em questão. Deste modo ,relativamente à posição do primeiro indice do array, os indices 0 a 2 referem-se a x,y,z.Os indices 3 a 5 referem-se a vx,vy,vz.
//Deste modo:
// param[0][i]=x da partícula i
// param[1][i]=y da partícula i
// param[2][i]=z da partícula i
// param[3][i]=vx da partícula i
// param[4][i]=vy da partícula i
// param[5][i]=vz da partícula i
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>
#include <sys/time.h>



// Number of particles
//static const int N = 16384;

static const int N = 20;



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

//-------------------------------------------------------------------------------

void advance_particles(double **param,const int N_particles,const double dt,int world_rank, int world_size)
{
  double *x,*y,*z,*vx,*vy,*vz,Fx,Fy,Fz;
  int my_first,my_last;
  x=param[0];
  y=param[1];
  z=param[2];
  vx=param[3];
  vy=param[4];
  vz=param[5];



  my_first =(world_rank)*N/world_size;
  my_last=(world_rank +1)*N/world_size;

 if(world_rank==0)
    {
      for(int i=0;i<N;i++)
	printf("%f\n",x[i]);
      printf("-------------");
    }
  

  for(int i =0;i<N;i++)
    {
      Fx=Fy=Fz=0;
      for(int j =my_first;j<my_last;j++)
	{
	  double dx=x[j]-x[i];
	  double dy=y[j]-y[i];
	  double dz=z[j]-z[i];
	  double dr2 = dx*dx + dy*dy + dz*dz + soft;
	  double r_dr3 = 1.0 / (dr2 * sqrt(dr2));
	  Fx+=dx*r_dr3;
	  Fy+=dy*r_dr3;
	  Fz+=dz*r_dr3;
	  
	}
        MPI_Allreduce(&Fx,&Fx,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
       MPI_Allreduce(&Fy,&Fy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
       MPI_Allreduce(&Fz,&Fz,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

 
      
      vx[i] += Fx * dt;
      vy[i] += Fy * dt;
      vz[i] += Fz * dt;
      // printf("teste");
    }

  
  for( int i =my_first ; i < my_last; i++)
    {
      x[i] += vx[i] * dt;
      y[i] += vy[i] * dt;
      z[i] += vz[i] * dt;
    }

 
  
 // MPI_Allreduce(&x,&param[0],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  // MPI_Allreduce(&y,&y,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  // MPI_Allreduce(&z,&z,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	
if(world_rank==0)
    {for(int i=0;i<N;i++)
	printf("%f\n",x[i]);
    }
	
}

//-----------------------------------------------------------------------------------------------------

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
  double energy;
  
   
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
 
  energy=kinetic_energy_local(param,N,world_rank,world_size);
 
  
  err=MPI_Allreduce(&energy,&energy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
 
  advance_particles(param,N,dt,world_rank,world_size);


  free(param[0]);
  free(param);
  MPI_Finalize();

  return 0;
}
