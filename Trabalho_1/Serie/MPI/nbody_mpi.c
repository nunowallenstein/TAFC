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
static const int N = 16384;

//static const int N = 100;



//static const int N = 10;

 // Number of iterations 
static const int N_iter = 100;

 // Timestep 
 static const double dt = 0.001; 

 // Softening factor 
 static const double soft = 1e-40; 

/* // Timing routine based on gettimeofday */
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



double **alloc_2d_init(int rows, int cols) {
  double *data = (double *)malloc(rows*cols*sizeof(double));
  double **array= (double **)malloc(rows*sizeof(double*));
  for (int i=0; i<rows; i++)
    array[i] = &(data[cols*i]);

  return array;
}

double kinetic_energy(double  **param, int N_particles, int world_rank,int n_proc)
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
  double energy,energy_total,start,stop;
  
   
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
  start=MPI_Wtime();

  for (int i=0;i<N_iter;i++)
    {
      energy=kinetic_energy(param,N,world_rank,world_size);
      MPI_Allreduce(&energy,&energy_total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      if(world_rank==0)
      	printf("i = %3d, kin = %g\n", i,energy_total );
      // advance_particles(param,N,dt,world_rank,world_size);

      double Fx,Fy,Fz,Fx_tot,Fy_tot,Fz_tot;
      int my_first,my_last;
      int num_elem_proc;
      int indexA;
      /* double *xaux,*yaux,*zaux;
      double *vxaux,*vyaux,*vzaux;
      xaux=param[0];
      yaux=param[1];
      zaux=param[2];
      vxaux=param[3];
      vyaux=param[4];
      vzaux=param[5];
      */
      int N_particles=N;
      my_first =(world_rank)*N_particles/world_size;
      my_last=(world_rank +1)*N_particles/world_size;
      num_elem_proc=N_particles/world_size;
      indexA=world_rank*num_elem_proc; 

      for(int i =0;i<num_elem_proc;i++)
	{
	  int iA=i+indexA;
	      	Fx=Fy=Fz=0;
	  for(int j =0;j<N_particles;j++)
	    {
	      double dx=param[0][j]-param[0][iA];
	      double dy=param[1][j]-param[1][iA];
	      double dz=param[2][j]-param[2][iA];
	      double dr2 = dx*dx + dy*dy + dz*dz + soft;
	      double r_dr3 = 1.0 / (dr2 * sqrt(dr2));
	      Fx+=dx*r_dr3;
	      Fy+=dy*r_dr3;
	      Fz+=dz*r_dr3;
	  
	    }

	  param[3][iA]=param[3][iA]+Fx*dt;
	  param[4][iA]=param[4][iA]+Fy*dt;
	  param[5][iA]=param[5][iA]+Fz*dt;
	 
	}
   
  
    
      for( int i =0; i < num_elem_proc;i++)
	{
	  int iA=i+indexA;
	  param[0][iA]=param[0][iA]+param[3][iA]*dt;
	  param[1][iA]=param[1][iA]+param[4][iA]*dt;
	  param[2][iA]=param[2][iA]+param[5][iA]*dt;
	  
	}
      
      	for (int i=0;i<world_size;i++)
	{
	  int a=i*num_elem_proc;
	MPI_Bcast(&(param[0][a]), num_elem_proc, MPI_DOUBLE, i, MPI_COMM_WORLD);
        MPI_Bcast(&(param[1][a]), num_elem_proc, MPI_DOUBLE, i, MPI_COMM_WORLD);
	MPI_Bcast(&(param[2][a]), num_elem_proc, MPI_DOUBLE, i, MPI_COMM_WORLD);
	MPI_Bcast(&(param[3][a]), num_elem_proc, MPI_DOUBLE, i, MPI_COMM_WORLD);
        MPI_Bcast(&(param[4][a]), num_elem_proc, MPI_DOUBLE, i, MPI_COMM_WORLD);
	MPI_Bcast(&(param[5][a]), num_elem_proc, MPI_DOUBLE, i, MPI_COMM_WORLD);
      
	}
	
    }
	  
  stop=MPI_Wtime();
  energy=kinetic_energy(param, N,world_rank,world_size);
  MPI_Allreduce(&energy,&energy_total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if(world_rank==0)
    {
      printf("---\n");
      printf("Total iterations = %d\n", N_iter );
      printf("Final kinetic energy: %g\n",energy_total );
      printf("Total execution time: %g s\n", stop - start);
    }
  free(param[0]);
  free(param);
  MPI_Finalize();

  return 0;
}
