#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>
#define NUM_THREADS 2



// Number of particles
static const int N = 16384;

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

typedef struct {
  double *x, *y, *z;
  double *vx, *vy, *vz;
} t_particles;

/**
 * Initializes buffer values in array buffer with random values between 0.0 and 1.0,
 * inclusive
 * @param buffer [description]
 * @param size   [description]
 */
void init_rand01( double* const buffer, const int size ) {
  const double r_rand_max = 1.0/RAND_MAX;
  for( int i = 0; i < size; i ++)
    buffer[i] = rand() * r_rand_max;
}

void init_particles(t_particles* const particles, const int n_part) {

  particles -> x  = (double *) malloc( n_part * sizeof(double) );
  particles -> y  = (double *) malloc( n_part * sizeof(double) );
  particles -> z  = (double *) malloc( n_part * sizeof(double) );
  particles -> vx = (double *) malloc( n_part * sizeof(double) );
  particles -> vy = (double *) malloc( n_part * sizeof(double) );
  particles -> vz = (double *) malloc( n_part * sizeof(double) );

  printf("Initializing %d particles...\n", n_part );

  srand(0);
  init_rand01( particles ->  x, n_part );
  init_rand01( particles ->  y, n_part );
  init_rand01( particles ->  z, n_part );
  init_rand01( particles -> vx, n_part );
  init_rand01( particles -> vy, n_part );
  init_rand01( particles -> vz, n_part );

  printf("Initialization complete.\n");
}

void cleanup_particles(t_particles* const particles){
  free( particles -> x);
  free( particles -> y);
  free( particles -> z);
  free( particles -> vx);
  free( particles -> vy);
  free( particles -> vz);
}

// Version 0
void advance_particles(t_particles* const particles, const int n_part, const double dt) {
  // omp_set_num_threads(NUM_THREADS);
  // Restrict pointers
  double* restrict x  = particles -> x;
  double* restrict y  = particles -> y;
  double* restrict z  = particles -> z;
  double* restrict vx = particles -> vx;
  double* restrict vy = particles -> vy;
  double* restrict vz = particles -> vz;

  // Advance velocities

  int i;
  double Fx, Fy, Fz;
 

  //#pragma omp for reduction(+:vx[i],vy[i],vz[i],Fx,Fy,Fz) schedule (auto) 
#pragma omp parallel
      {
#pragma omp for reduction(+:vx[:n_part],vy[:n_part],vz[:n_part],Fx,Fy,Fz) schedule (auto)	
  for( int i = 0; i < n_part; i++)
    {
      Fx=Fy=Fz=0;	

      //#pragma omp for reduction(+:Fx,Fy,Fz) schedule (auto)	
	for( int j = 0; j < n_part; j++)
	  {
	    const double dx = x[j] - x[i];
	    const double dy = y[j] - y[i];
	    const double dz = z[j] - z[i];

	    // The softening (soft) value avoids a division by zero
	    // if two particles happen to have the same position
	    const double dr2 = dx*dx + dy*dy + dz*dz + soft;
	    const double r_dr3 = 1.0 / (dr2 * sqrt(dr2));

	    Fx += dx * r_dr3;
	    Fy += dy * r_dr3;
	    Fz += dz * r_dr3;
	  }
         
      // Since m = 1, F = a
      vx[i] += Fx * dt;
      vy[i] += Fy * dt;
      vz[i] += Fz * dt;
    };
  }	

  // Advance positions
#pragma parallel for reduction(+:x[:n_part],y[:n_part],z[:n_part]) schedule (auto)
  for( int i = 0; i < n_part; i++){
    x[i] += vx[i] * dt;
    y[i] += vy[i] * dt;
    z[i] += vz[i] * dt;
  };
}

double kinetic_energy( t_particles* const particles, const int n_part ){
  omp_set_num_threads(NUM_THREADS);
  double energy ;
#pragma parallel for reduction(+:energy) schedule (auto)
  for(int i = 0; i < n_part; i++) {
    energy += particles->vx[i] * particles->vx[i] +
      particles->vy[i] * particles->vy[i] +
      particles->vz[i] * particles->vz[i];
  }

  return 0.5*energy;
}

int main (int argc, const char * argv[]) {
  omp_set_num_threads(NUM_THREADS);
  omp_set_nested(1);
  t_particles particles;

  // Initialize particle positions and velocities
  init_particles( &particles, N );

  // Time execution
  const double start = omp_get_wtime();

  // Advance particles
  for( int i = 0; i < N_iter; i ++ )
    {
      printf("i = %3d, kin = %g\n", i, kinetic_energy( &particles, N ));
      advance_particles( &particles, N, dt );
      
    }

  const double stop = omp_get_wtime();

  // Print execution time
  printf("---\n");
  printf("Total iterations = %d\n", N_iter );
  printf("Final kinetic energy: %g\n", kinetic_energy(&particles, N));
  printf("Total execution time: %g s\n", stop - start);

  cleanup_particles(&particles);
}
