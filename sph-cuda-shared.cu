/****************************************************************************
 *
 * sph.c -- Smoothed Particle Hydrodynamics
 *
 * https://github.com/cerrno/mueller-sph
 *
 * Copyright (C) 2016 Lucas V. Schuermann
 * Copyright (C) 2022 Moreno Marzolla
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 ****************************************************************************/
#include "hpc.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* "Particle-Based Fluid Simulation for Interactive Applications" by
   Müller et al. solver parameters */

const float Gx = 0.0, Gy = -10.0;   // external (gravitational) forces
const float REST_DENS = 300;    // rest density
const float GAS_CONST = 2000;   // const for equation of state
const float H = 16;             // kernel radius
const float EPS = 16;           // equal to H
const float MASS = 2.5;         // assume all particles have the same mass
const float VISC = 200;         // viscosity constant
const float DT = 0.0007;        // integration timestep
const float BOUND_DAMPING = -0.5;

// rendering projection parameters
// (the following ought to be "const float", but then the compiler
// would give an error because VIEW_WIDTH and VIEW_HEIGHT are
// initialized with non-literal expressions)

const int MAX_PARTICLES = 20000;
// Larger window size to accommodate more particles
#define WINDOW_WIDTH 3000
#define WINDOW_HEIGHT 2000

const int DAM_PARTICLES = 500;

const float VIEW_WIDTH = 1.5 * WINDOW_WIDTH;
const float VIEW_HEIGHT = 1.5 * WINDOW_HEIGHT;

#define PRINT_AVERANGE 10

#define BLKDIM 1024

/* Particle data structure; stores position, velocity, and force for
   integration stores density (rho) and pressure values for SPH.

   You may choose a different layout of the particles[] data structure
   to suit your needs. */
typedef struct {
    float x, y;         // position
    float vx, vy;       // velocity
    float fx, fy;       // force
    float rho, p;       // density, pressure
} particle_t;


float *pos_x, *pos_y;
float *vx, *vy;
float *fx, *fy;
float *rho, *p; 

int n_particles = 0;    // number of currently active particles

#define SHARED_MEM_PER_BLOCK 49152

/**
 * Return a random value in [a, b]
 */
float randab(float a, float b)
{
    return a + (b-a)*rand() / (float)(RAND_MAX);
}

/**
 * Set initial position of particle `*p` to (x, y); initialize all
 * other attributes to default values (zeros).
 */
void init_particle( int index , float x, float y )
{
    pos_x[index] = x;
    pos_y[index] = y;
    vx[index] = vy[index] = 0.0;
    fx[index] = fy[index] = 0.0;
    rho[index] = 0.0;
    p[index] = 0.0;
}

/**
 * Return nonzero iff (x, y) is within the frame
 */
int is_in_domain( float x, float y )
{
    return ((x < VIEW_WIDTH - EPS) &&
            (x > EPS) &&
            (y < VIEW_HEIGHT - EPS) &&
            (y > EPS));
}

/**
 * Initialize the SPH model with `n` particles. The caller is
 * responsible for allocating the `particles[]` array of size
 * `MAX_PARTICLES`.
 *
 * DO NOT parallelize this function, since it calls rand() which is
 * not thread-safe.
 *
 * For MPI and OpenMP: only the master must initialize the domain;
 *
 * For CUDA: the CPU must initialize the domain.
 */
void init_sph( int n )
{
    n_particles = 0;
    //printf("Initializing with %d particles\n", n);

    for (float y = EPS; y < VIEW_HEIGHT - EPS; y += H) {
        for (float x = EPS; x <= VIEW_WIDTH * 0.8f; x += H) {
            if (n_particles < n) {
                float jitter = rand() / (float)RAND_MAX;
                init_particle(n_particles, x+jitter, y);
                n_particles++;
            } else {
                return;
            }
        }
    }
    assert(n_particles == n);
}

/**
 ** You may parallelize the following four functions
 **/

__global__ void compute_density_pressure( float* d_rho, float* d_pos_x, float * d_pos_y, float * d_p, int n_particles)
{
    const int index_particle = threadIdx.x + blockIdx.x * blockDim.x;
    const int lindex = threadIdx.x;
    const int FLOAT_PER_SHARED_MEM = SHARED_MEM_PER_BLOCK / sizeof(float);
    
    
    __shared__ float sh_pos_x[FLOAT_PER_SHARED_MEM/2];
    __shared__ float sh_pos_y[FLOAT_PER_SHARED_MEM/2];

    const float HSQ = H * H;    // radius^2 for optimization

    /* Smoothing kernels defined in Muller and their gradients adapted
    to 2D per "SPH Based Shallow Water Simulation" by Solenthaler
    et al. */
    const float POLY6 = 4.0 / (M_PI * pow(H, 8));

    d_rho[index_particle] = 0.0;
    
    // per ogni particella memorizzi 2 float
    // numero di volte di cui devi fare una copia nella shared per tutto il kernel
    const int repetitions = (n_particles * 2 + FLOAT_PER_SHARED_MEM - 1) / FLOAT_PER_SHARED_MEM;
    const int max_particles_to_copy = FLOAT_PER_SHARED_MEM / 2;

    printf("[idx: %4d] [x: %f] [y: %f] [rho: %f] [p: %f]\n",
            index_particle,
            d_pos_x[index_particle],
            d_pos_y[index_particle],
            d_rho[index_particle],
            d_p[index_particle]
            );
    
    for (int r = 0; r < repetitions;  r++) {
        int end_copy = max_particles_to_copy;

        if (r == repetitions - 1) {
            end_copy = n_particles - max_particles_to_copy * r;
        }

        int copy_shift = 0;
        while (copy_shift * BLKDIM + lindex < end_copy) {
            sh_pos_x[copy_shift * BLKDIM + lindex] = d_pos_x[r * max_particles_to_copy + copy_shift * BLKDIM + lindex];
            sh_pos_y[copy_shift * BLKDIM + lindex] = d_pos_y[r * max_particles_to_copy + copy_shift * BLKDIM + lindex];
            copy_shift++;
        }

        __syncthreads();

        if (index_particle < n_particles) {
            for (int j = 0; j < end_copy; j++) {

                const float dx = sh_pos_x[j] - d_pos_x[index_particle];
                const float dy = sh_pos_y[j] - d_pos_y[index_particle];
                const float d2 = dx*dx + dy*dy;

                if (d2 < HSQ) {
                    d_rho[index_particle] += MASS * POLY6 * pow(HSQ - d2, 3.0);
                }
            }

            d_p[index_particle] = GAS_CONST * (d_rho[index_particle] - REST_DENS);
        }

        __syncthreads();
    }
}

__global__ void compute_forces( float* d_rho, float* d_pos_x, float * d_pos_y, float * d_p, float* d_vx, float* d_vy, float* d_fx, float* d_fy, int n_particles )
{
    const int index_particle = threadIdx.x + blockIdx.x * blockDim.x;
    /* Smoothing kernels defined in Muller and their gradients adapted
       to 2D per "SPH Based Shallow Water Simulation" by Solenthaler
       et al. */
    const int lindex = threadIdx.x;
    const int FLOAT_PER_SHARED_MEM = SHARED_MEM_PER_BLOCK / sizeof(float);

    const float SPIKY_GRAD = -10.0 / (M_PI * pow(H, 5));
    const float VISC_LAP = 40.0 / (M_PI * pow(H, 5));
    const float EPS = 1e-6;

    float fpress_x = 0.0, fpress_y = 0.0;
    float fvisc_x = 0.0, fvisc_y = 0.0;

    __shared__ float sh_pos_x[FLOAT_PER_SHARED_MEM/2];
    __shared__ float sh_pos_y[FLOAT_PER_SHARED_MEM/2];

    const int repetitions = (n_particles * 2 + FLOAT_PER_SHARED_MEM - 1) / FLOAT_PER_SHARED_MEM;
    const int max_particles_to_copy = FLOAT_PER_SHARED_MEM / 2;

    printf("[idx: %4d] [vx: %f] [vy: %f] [x: %f] [y: %f] [fx: %f] [fy: %f] [rho: %f] [p: %f]\n",
                index_particle,
                d_vx[index_particle],
                d_vy[index_particle],
                d_pos_x[index_particle],
                d_pos_y[index_particle],
                d_fx[index_particle],
                d_fy[index_particle],
                d_rho[index_particle],
                d_p[index_particle]
                );
    
    for (int r = 0; r < repetitions;  r++) {
        int end_copy = max_particles_to_copy;

        if (r == repetitions - 1) {
            end_copy = n_particles - max_particles_to_copy * r;
        }

        int copy_shift = 0;
        while (copy_shift * BLKDIM + lindex < end_copy) {
            sh_pos_x[copy_shift * BLKDIM + lindex] = d_pos_x[r * max_particles_to_copy + copy_shift * BLKDIM + lindex];
            sh_pos_y[copy_shift * BLKDIM + lindex] = d_pos_y[r * max_particles_to_copy + copy_shift * BLKDIM + lindex];
            copy_shift++;
        }

        __syncthreads();

        if (index_particle < n_particles)  {
            for (int j=0; j< end_copy; j++) {

                if (index_particle == r * max_particles_to_copy + j)
                    continue;

                const float dx = sh_pos_x[j] - d_pos_x[index_particle];
                const float dy = sh_pos_y[j] - d_pos_y[index_particle];
                const float dist = hypotf(dx, dy) + EPS; // avoids division by zero later on

                if (dist < H) {
                    const float norm_dx = dx / dist;
                    const float norm_dy = dy / dist;
                    // compute pressure force contribution
                    fpress_x += -norm_dx * MASS * (d_p[index_particle] + d_p[j]) / (2 * d_rho[j]) * SPIKY_GRAD * pow(H - dist, 3);
                    fpress_y += -norm_dy * MASS * (d_p[index_particle] + d_p[j]) / (2 * d_rho[j]) * SPIKY_GRAD * pow(H - dist, 3);
                    // compute viscosity force contribution
                    fvisc_x += VISC * MASS * (d_vx[j] - d_vx[index_particle]) / d_rho[j] * VISC_LAP * (H - dist);
                    fvisc_y += VISC * MASS * (d_vy[j] - d_vy[index_particle]) / d_rho[j] * VISC_LAP * (H - dist);
                }
            }
        }

        __syncthreads();


    }
    const float fgrav_x = Gx * MASS / d_rho[index_particle];
    const float fgrav_y = Gy * MASS / d_rho[index_particle];
    d_fx[index_particle] = fpress_x + fvisc_x + fgrav_x;
    d_fy[index_particle] = fpress_y + fvisc_y + fgrav_y;
    
}

__global__ void integrate( float* d_rho, float* d_x, float * d_y, float* d_vx, float* d_vy, float* d_fx, float* d_fy, int n_particles )
{
    const int index_particle = threadIdx.x + blockIdx.x * blockDim.x;
    if (index_particle < n_particles) {
        printf("[idx: %4d] [vx: %f] [vy: %f] [x: %f] [y: %f] [fx: %f] [fy: %f] [rho: %f]\n",
                index_particle,
                d_vx[index_particle],
                d_vy[index_particle],
                d_x[index_particle],
                d_y[index_particle],
                d_fx[index_particle],
                d_fy[index_particle],
                d_rho[index_particle]
                );
        // forward Euler integration
        d_vx[index_particle] += DT * d_fx[index_particle] / d_rho[index_particle];
        d_vy[index_particle] += DT * d_fy[index_particle] / d_rho[index_particle];
        d_x[index_particle] += DT * d_vx[index_particle];
        d_y[index_particle] += DT * d_vy[index_particle];

        // enforce boundary conditions
        if (d_x[index_particle] - EPS < 0.0) {
            d_vx[index_particle] *= BOUND_DAMPING;
            d_x[index_particle] = EPS;
        }
        if (d_x[index_particle] + EPS > VIEW_WIDTH) {
            d_vx[index_particle] *= BOUND_DAMPING;
            d_x[index_particle] = VIEW_WIDTH - EPS;
        }
        if (d_y[index_particle] - EPS < 0.0) {
            d_vy[index_particle] *= BOUND_DAMPING;
            d_y[index_particle] = EPS;
        }
        if (d_y[index_particle] + EPS > VIEW_HEIGHT) {
            d_vy[index_particle] *= BOUND_DAMPING;
            d_y[index_particle] = VIEW_HEIGHT - EPS;
        }
    }
}

__global__ void reduction(float* d_vx, float* d_vy, int n, float * d_sums) {
    const int index = threadIdx.x + blockIdx.x * blockDim.x;

    printf("[idx: %4d] [vx: %f] [vy: %f]\n", index, d_vx[index], d_vy[index]);

    /* reduction of averange velocity */
    if (index < n) {
        __shared__ float temp[BLKDIM];
        const int lindex = threadIdx.x;
        const int bindex = blockIdx.x;
        int bsize = blockDim.x / 2;
        temp[lindex] = hypot(d_vx[index], d_vy[index]) / n;

        __syncthreads();
        while ( bsize > 0 ) {
            if ( lindex < bsize ) {
                temp[lindex] += temp[lindex + bsize];
            }
            bsize = bsize / 2;
            __syncthreads();
        }
        if ( 0 == lindex ) {
            d_sums[bindex] = temp[0];
        }
    }
}

#define MAX_BLOCK (MAX_PARTICLES + BLKDIM - 1)/BLKDIM

int main(int argc, char **argv)
{
    srand(1234);

    int n = DAM_PARTICLES;
    int nsteps = 50;

    if (argc > 3) {
        fprintf(stderr, "Usage: %s [nparticles [nsteps]]\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (argc > 1) {
        n = atoi(argv[1]);
    }

    if (argc > 2) {
        nsteps = atoi(argv[2]);
    }

    if (n > MAX_PARTICLES) {
        fprintf(stderr, "FATAL: the maximum number of particles is %d\n", MAX_PARTICLES);
        return EXIT_FAILURE;
    }

    pos_x = (float *) malloc(n * sizeof(float)); assert( pos_x != NULL );
    pos_y = (float *) malloc(n * sizeof(float)); assert( pos_y != NULL );
    vx = (float *) malloc(n * sizeof(float)); assert( vx != NULL );
    vy = (float *) malloc(n * sizeof(float)); assert( vy != NULL );
    fx = (float *) malloc(n * sizeof(float)); assert( fx != NULL );
    fy = (float *) malloc(n * sizeof(float)); assert( fy != NULL );
    rho = (float *) malloc(n * sizeof(float)); assert( rho != NULL );
    p = (float *) malloc(n * sizeof(float)); assert( p != NULL );
    
    float *d_pos_x, *d_pos_y;
    float *d_vx, *d_vy;
    float *d_fx, *d_fy;
    float *d_rho, *d_p; 


    float h_sums[MAX_BLOCK];
    float *d_sums;

    int block_num = (n + BLKDIM - 1)/BLKDIM;

    init_sph(n);

    cudaMalloc((void **) &d_pos_x, sizeof(float) * n);
    cudaMemcpy(d_pos_x, pos_x, sizeof(float) * n, cudaMemcpyHostToDevice);

    cudaMalloc((void **) &d_pos_y, sizeof(float) * n);
    cudaMemcpy(d_pos_y, pos_y, sizeof(float) * n, cudaMemcpyHostToDevice);

    cudaMalloc((void **) &d_vx, sizeof(float) * n);
    cudaMemcpy(d_vx, vx, sizeof(float) * n, cudaMemcpyHostToDevice);

    cudaMalloc((void **) &d_vy, sizeof(float) * n);
    cudaMemcpy(d_vy, vy, sizeof(float) * n, cudaMemcpyHostToDevice);

    cudaMalloc((void **) &d_fx, sizeof(float) * n);
    cudaMemcpy(d_fx, fx, sizeof(float) * n, cudaMemcpyHostToDevice);

    cudaMalloc((void **) &d_fy, sizeof(float) * n);
    cudaMemcpy(d_fy, fy, sizeof(float) * n, cudaMemcpyHostToDevice);

    cudaMalloc((void **) &d_rho, sizeof(float) * n);
    cudaMemcpy(d_rho, rho, sizeof(float) * n, cudaMemcpyHostToDevice);

    cudaMalloc((void **) &d_p, sizeof(float) * n);
    cudaMemcpy(d_p, p, sizeof(float) * n, cudaMemcpyHostToDevice);
    
    cudaMalloc((void **) &d_sums, block_num * sizeof(float));

    double loop_start = hpc_gettime();
    
    for (int s=0; s<nsteps; s++) {
        double start = hpc_gettime();

        printf("Pressione e Densita:\n");

        compute_density_pressure<<<block_num, BLKDIM>>>(d_rho, d_pos_x, d_pos_y, d_p, n);
        
        cudaDeviceSynchronize();

        printf("Forze:\n");

        compute_forces<<<block_num, BLKDIM>>>(d_rho, d_pos_x, d_pos_y, d_p, d_vx, d_vy, d_fx, d_fy, n);

        cudaDeviceSynchronize();

        printf("Integrazione");

        integrate<<<block_num, BLKDIM>>>(d_rho, d_pos_x, d_pos_y, d_vx, d_vy, d_fx, d_fy, n);

        cudaDeviceSynchronize();

        printf("Reduction");

        reduction<<<block_num, BLKDIM>>>(d_vx, d_vy, n, d_sums);
        /* the average velocities MUST be computed at each step, even
        if it is not shown (to ensure constant workload per
        iteration) */
        cudaMemcpy(h_sums, d_sums, block_num * sizeof(float), cudaMemcpyDeviceToHost);
        
        float avg = 0.0;
        
        for (int i = 0; i < block_num; i++)
            avg += h_sums[i];
        
        double end = hpc_gettime() - start;

        if (s % PRINT_AVERANGE == 0){
            printf("step %5d, avgV=%f, took: %fs\n", s, avg, end);
            //printf("%f;", avg);
            //for (int i = 0; i < MAX_BLOCK; i++)
            //    printf("%f ", h_sums[i]);
            //printf("\n");
        }
    }

    double loop_end = hpc_gettime() - loop_start;
    printf("took: %fs\n", loop_end);

    cudaFree(d_rho);
    cudaFree(d_pos_x);
    cudaFree(d_pos_y);
    cudaFree(d_p);
    cudaFree(d_vx);
    cudaFree(d_vy);
    cudaFree(d_fx);
    cudaFree(d_fy);

    cudaFree(d_sums);


    free(rho);
    free(pos_x);
    free(pos_y);
    free(p);
    free(vx);
    free(vy);
    free(fx);
    free(fy);
    return EXIT_SUCCESS;
}
