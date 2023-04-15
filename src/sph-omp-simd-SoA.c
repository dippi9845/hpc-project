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
#include <omp.h>

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


#define DYNAMIC_SIZE 15

const int DAM_PARTICLES = 500;

const float VIEW_WIDTH = 1.5 * WINDOW_WIDTH;
const float VIEW_HEIGHT = 1.5 * WINDOW_HEIGHT;

typedef float v4f __attribute__ ((vector_size (16)));
#define VLEN (sizeof(v4f) /sizeof(float))

int particles_num;

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

//particle_t *particles;

float *pos_x, *pos_y;
float *vx, *vy;
float *fx, *fy;
float *rho, *p; 

int n_particles = 0;    // number of currently active particles

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
void init_particle( int index, float x, float y )
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

void compute_density_pressure( void )
{
    const float HSQ = H * H;    // radius^2 for optimization

    /* Smoothing kernels defined in Muller and their gradients adapted
       to 2D per "SPH Based Shallow Water Simulation" by Solenthaler
       et al. */
    const float POLY6 = 4.0 / (M_PI * pow(H, 8));
    v4f acc_rho = {0.f, 0.f, 0.f, 0.f};

    #pragma omp parallel for schedule(dynamic, DYNAMIC_SIZE) default(none) shared(n_particles, pos_x, pos_y, vx, vy, fx, fy, rho, p) firstprivate(acc_rho)
    for (int i=0; i<n_particles; i++) {
        //particle_t *pi = &particles[i];
        rho[i] = 0.0;
        int near = 0;
        for (int j=0; j<n_particles; j++) {
            //const particle_t *pj = &particles[j];

            const float dx = pos_x[j] - pos_x[i];
            const float dy = pos_y[j] - pos_y[i];
            const float d2 = dx*dx + dy*dy;

            if (d2 < HSQ) {
                acc_rho[near] = MASS * POLY6 * pow(HSQ - d2, 3.0);
                near++;
                if (near == VLEN) {
                    rho[i] += acc_rho[0] + acc_rho[1] + acc_rho[2] + acc_rho[3];
                    near = 0; 
                }
            }
        }
        
        /* handle remaining */
        for (int index = 0; index < near; index++) {
            rho[i] += acc_rho[index];
        }
        /* end of simd computation */
        p[i] = GAS_CONST * (rho[i] - REST_DENS);
    }
}

void compute_forces( void )
{
    /* Smoothing kernels defined in Muller and their gradients adapted
       to 2D per "SPH Based Shallow Water Simulation" by Solenthaler
       et al. */
    const float SPIKY_GRAD = -10.0 / (M_PI * pow(H, 5));
    const float VISC_LAP = 40.0 / (M_PI * pow(H, 5));
    const float EPS = 1e-6;

    v4f acc_press_x = {0.f, 0.f, 0.f, 0.f};
    v4f acc_press_y = {0.f, 0.f, 0.f, 0.f};
    v4f acc_visc_x = {0.f, 0.f, 0.f, 0.f};
    v4f acc_visc_y = {0.f, 0.f, 0.f, 0.f};

    #pragma omp parallel for schedule(dynamic, DYNAMIC_SIZE) default(none) shared(n_particles, pos_x, pos_y, vx, vy, fx, fy, rho, p) firstprivate(acc_press_x, acc_press_y, acc_visc_x, acc_visc_y)
    for (int i=0; i<n_particles; i++) {
        //particle_t *pi = &particles[i];
        float fpress_x = 0.0, fpress_y = 0.0;
        float fvisc_x = 0.0, fvisc_y = 0.0;
        int near = 0;

        for (int j=0; j<n_particles; j++) {
            //const particle_t *pj = &particles[j];

            if (i == j)
                continue;

            const float dx = pos_x[j] - pos_x[i];
            const float dy = pos_y[j] - pos_y[i];
            const float dist = hypotf(dx, dy) + EPS; // avoids division by zero later on

            if (dist < H) {
                const float norm_dx = dx / dist;
                const float norm_dy = dy / dist;
                
                // compute pressure force contribution
                acc_press_x[near] = -norm_dx * MASS * (p[i] + p[j]) / (2 * rho[j]) * SPIKY_GRAD * pow(H - dist, 3);
                acc_press_y[near] = -norm_dy * MASS * (p[i] + p[j]) / (2 * rho[j]) * SPIKY_GRAD * pow(H - dist, 3);
                // compute viscosity force contribution
                acc_visc_x[near] = VISC * MASS * (vx[j] - vx[i]) / rho[j] * VISC_LAP * (H - dist);
                acc_visc_y[near] = VISC * MASS * (vy[j] - vy[i]) / rho[j] * VISC_LAP * (H - dist);
                
                near++;
                if (near == VLEN) {
                    near = 0;
                    fpress_x += acc_press_x[0] + acc_press_x[1] + acc_press_x[2] + acc_press_x[3];
                    fpress_y += acc_press_y[0] + acc_press_y[1] + acc_press_y[2] + acc_press_y[3];
                    fvisc_x += acc_visc_x[0] + acc_visc_x[1] + acc_visc_x[2] + acc_visc_x[3];
                    fvisc_y += acc_visc_y[0] + acc_visc_y[1] + acc_visc_y[2] + acc_visc_y[3];
                }
            }
        }

        const float fgrav_x = Gx * MASS / rho[i];
        const float fgrav_y = Gy * MASS / rho[i];


        for (int index = 0; index < near; index++) {
            fpress_x += acc_press_x[index];
            fpress_y += acc_press_y[index];
            fvisc_x += acc_visc_x[index];
            fvisc_y += acc_visc_y[index];
        }
        /*-------------------------------------*/

        fx[i] = fpress_x + fvisc_x + fgrav_x;
        fy[i] = fpress_y + fvisc_y + fgrav_y;
    }
}

void integrate( void )
{
    #pragma omp parallel for schedule(dynamic, DYNAMIC_SIZE) default(none) shared(n_particles, pos_x, pos_y, vx, vy, fx, fy, rho, p)
    for (int i=0; i<n_particles; i++) {
        // forward Euler integration
        vx[i] += DT * fx[i] / rho[i];
        vy[i] += DT * fy[i] / rho[i];
        pos_x[i] += DT * vx[i];
        pos_y[i] += DT * vy[i];

        // enforce boundary conditions
        if (pos_x[i] - EPS < 0.0) {
            vx[i] *= BOUND_DAMPING;
            pos_x[i] = EPS;
        }
        if (pos_x[i] + EPS > VIEW_WIDTH) {
            vx[i] *= BOUND_DAMPING;
            pos_x[i] = VIEW_WIDTH - EPS;
        }
        if (pos_y[i] - EPS < 0.0) {
            vy[i] *= BOUND_DAMPING;
            pos_y[i] = EPS;
        }
        if (pos_y[i] + EPS > VIEW_HEIGHT) {
            vy[i] *= BOUND_DAMPING;
            pos_y[i] = VIEW_HEIGHT - EPS;
        }
    }
}


float avg_velocities( void )
{
    double result = 0.0;
    #pragma omp parallel for reduction(+:result)
    for (int i=0; i<n_particles; i++) {
        /* the hypot(x,y) function is equivalent to sqrt(x*x +
           y*y); */
        result += hypot(vx[i], vy[i]) / n_particles;
    }
    return result;
}


void update( void ) {

    compute_density_pressure();
    compute_forces();
    integrate();
}

int main(int argc, char **argv)
{
    srand(1234);

    //particles = (particle_t*)malloc(MAX_PARTICLES * sizeof(*particles));
    
    pos_x = (float *) malloc(MAX_PARTICLES * sizeof(float)); assert( pos_x != NULL );
    pos_y = (float *) malloc(MAX_PARTICLES * sizeof(float)); assert( pos_y != NULL );
    vx = (float *) malloc(MAX_PARTICLES * sizeof(float)); assert( vx != NULL );
    vy = (float *) malloc(MAX_PARTICLES * sizeof(float)); assert( vy != NULL );
    fx = (float *) malloc(MAX_PARTICLES * sizeof(float)); assert( fx != NULL );
    fy = (float *) malloc(MAX_PARTICLES * sizeof(float)); assert( fy != NULL );
    rho = (float *) malloc(MAX_PARTICLES * sizeof(float)); assert( rho != NULL );
    p = (float *) malloc(MAX_PARTICLES * sizeof(float)); assert( p != NULL );

    //assert( particles != NULL );

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

    particles_num = n;

    init_sph(n);
    double st = hpc_gettime();
    for (int s=0; s<nsteps; s++) {
        update();
        
        /* the average velocities MUST be computed at each step, even
    if it is not shown (to ensure constant workload per
    iteration) */
        const float avg = avg_velocities();

        if (s % 10 == 0) {
            printf("step %5d, avgV=%f\n", s, avg);
        }
    }
    printf("%f\n", hpc_gettime() - st);

    //free(particles);
    free(pos_x);
    free(pos_y);
    free(vx);
    free(vy);
    free(fx);
    free(fy);
    free(rho);
    free(p);
    return EXIT_SUCCESS;
}
