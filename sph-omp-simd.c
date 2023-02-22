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


const int DAM_PARTICLES = 500;

const float VIEW_WIDTH = 1.5 * WINDOW_WIDTH;
const float VIEW_HEIGHT = 1.5 * WINDOW_HEIGHT;

float *near_rho;
float *near_press_x;
float *near_press_y;
float *near_visc_x;
float *near_visc_y;

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

particle_t *particles;
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
void init_particle( particle_t *p, float x, float y )
{
    p->x = x;
    p->y = y;
    p->vx = p->vy = 0.0;
    p->fx = p->fy = 0.0;
    p->rho = 0.0;
    p->p = 0.0;
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
    printf("Initializing with %d particles\n", n);

    for (float y = EPS; y < VIEW_HEIGHT - EPS; y += H) {
        for (float x = EPS; x <= VIEW_WIDTH * 0.8f; x += H) {
            if (n_particles < n) {
                float jitter = rand() / (float)RAND_MAX;
                init_particle(particles + n_particles, x+jitter, y);
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

void compute_density_pressure( size_t start, size_t end, size_t step, size_t my_id )
{
    const float HSQ = H * H;    // radius^2 for optimization

    /* Smoothing kernels defined in Muller and their gradients adapted
       to 2D per "SPH Based Shallow Water Simulation" by Solenthaler
       et al. */
    const float POLY6 = 4.0 / (M_PI * pow(H, 8));
    
    for (int i = start; i < end; i += step) {
        particle_t *pi = &particles[i];
        pi->rho = 0.0;
        int near = 0;
        for (int j=0; j<n_particles; j++) {
            const particle_t *pj = &particles[j];

            const float dx = pj->x - pi->x;
            const float dy = pj->y - pi->y;
            const float d2 = dx*dx + dy*dy;

            if (d2 < HSQ) {
                near_rho[near + particles_num * my_id] = MASS * POLY6 * pow(HSQ - d2, 3.0);
                near++;
            }
        }

        /* evaluate rho 
            pi->rho += MASS * POLY6 * pow(HSQ - d2, 3.0);
        */
        int index = 0;
        
        if (near > VLEN) {
            //printf("density near : %d\n", near);
            v4f acc = {0.0, 0.0, 0.0, 0.0};
            v4f *vv = (v4f*)near_rho;

            for (; index < near - VLEN + 1; index+= VLEN) {
                acc += *vv;
                vv++;
            }
            
            pi->rho = acc[0] + acc[1] + acc[2] + acc[3];
            
        }
        
        for (; index < near; index++) {
            pi->rho += near_rho[index + particles_num * my_id];
        }
        
        /* end of simd computation */
        pi->p = GAS_CONST * (pi->rho - REST_DENS);
    }
}

void compute_forces( size_t start, size_t end, size_t step, size_t my_id )
{
    /* Smoothing kernels defined in Muller and their gradients adapted
       to 2D per "SPH Based Shallow Water Simulation" by Solenthaler
       et al. */
    const float SPIKY_GRAD = -10.0 / (M_PI * pow(H, 5));
    const float VISC_LAP = 40.0 / (M_PI * pow(H, 5));
    const float EPS = 1e-6;

    for (int i = start; i < end; i += step) {
        particle_t *pi = &particles[i];
        float fpress_x = 0.0, fpress_y = 0.0;
        float fvisc_x = 0.0, fvisc_y = 0.0;
        int near = 0;

        for (int j=0; j<n_particles; j++) {
            const particle_t *pj = &particles[j];

            if (pi == pj)
                continue;

            const float dx = pj->x - pi->x;
            const float dy = pj->y - pi->y;
            const float dist = hypotf(dx, dy) + EPS; // avoids division by zero later on

            if (dist < H) {
                const float norm_dx = dx / dist;
                const float norm_dy = dy / dist;
                // compute pressure force contribution
                near_press_x[near + particles_num * my_id] = -norm_dx * MASS * (pi->p + pj->p) / (2 * pj->rho) * SPIKY_GRAD * pow(H - dist, 3);
                near_press_y[near + particles_num * my_id] = -norm_dy * MASS * (pi->p + pj->p) / (2 * pj->rho) * SPIKY_GRAD * pow(H - dist, 3);
                // compute viscosity force contribution
                near_visc_x[near + particles_num * my_id] = VISC * MASS * (pj->vx - pi->vx) / pj->rho * VISC_LAP * (H - dist);
                near_visc_y[near + particles_num * my_id] = VISC * MASS * (pj->vy - pi->vy) / pj->rho * VISC_LAP * (H - dist);
                near++;
            }
        }
        /*
            fpress_x += -norm_dx * MASS * (pi->p + pj->p) / (2 * pj->rho) * SPIKY_GRAD * pow(H - dist, 3);
            fpress_y += -norm_dy * MASS * (pi->p + pj->p) / (2 * pj->rho) * SPIKY_GRAD * pow(H - dist, 3);
            fvisc_x += VISC * MASS * (pj->vx - pi->vx) / pj->rho * VISC_LAP * (H - dist);
            fvisc_y += VISC * MASS * (pj->vy - pi->vy) / pj->rho * VISC_LAP * (H - dist);
        */
        const float fgrav_x = Gx * MASS / pi->rho;
        const float fgrav_y = Gy * MASS / pi->rho;

        //int index = 0; // index for all simd operations
        
        /* if is possible to use simd ops */
        /*
        if (near > VLEN) {
            //printf("press near : %d\n", near);
            v4f pres_x = {0.0, 0.0, 0.0, 0.0};
            v4f pres_y = {0.0, 0.0, 0.0, 0.0};
            v4f visc_x = {0.0, 0.0, 0.0, 0.0};
            v4f visc_y = {0.0, 0.0, 0.0, 0.0};

            v4f *vv_press_x = (v4f*)near_press_x;
            v4f *vv_press_y = (v4f*)near_press_y;
            v4f *vv_visc_x = (v4f*)near_visc_x;
            v4f *vv_visc_y = (v4f*)near_visc_y;

            // TODO: problemi di cache provare anche 4 loop separati
            for (; index < near - VLEN + 1; index+= VLEN) {
                pres_x += *vv_press_x;
                pres_y += *vv_press_y;
                visc_x += *vv_visc_x;
                visc_y += *vv_visc_y;

                vv_press_x++;
                vv_press_y++;
                vv_visc_x++;
                vv_visc_y++;

            }
            
            
            fpress_x = pres_x[0] + pres_x[1] + pres_x[2] + pres_x[3];
            fpress_y = pres_y[0] + pres_y[1] + pres_y[2] + pres_y[3];
            fvisc_x = visc_x[0] + visc_x[1] + visc_x[2] + visc_x[3];
            fvisc_y = visc_y[0] + visc_y[1] + visc_y[2] + visc_y[3];

        }
        */
           
        /* remaining of everything TODO: cache pure qua */
        for (int index = 0; index < near; index++) {
            fpress_x += near_press_x[index + particles_num * my_id];
            fpress_y += near_press_y[index + particles_num * my_id];
            fvisc_x += near_visc_x[index + particles_num * my_id];
            fvisc_y += near_visc_y[index + particles_num * my_id];
        }
        /*-------------------------------------*/

        pi->fx = fpress_x + fvisc_x + fgrav_x;
        pi->fy = fpress_y + fvisc_y + fgrav_y;
    }
}

void integrate( size_t start, size_t end, size_t step )
{
    for (size_t i = start; i < end; i += step) {
        particle_t *p = &particles[i];
        // forward Euler integration
        p->vx += DT * p->fx / p->rho;
        p->vy += DT * p->fy / p->rho;
        p->x += DT * p->vx;
        p->y += DT * p->vy;

        // enforce boundary conditions
        if (p->x - EPS < 0.0) {
            p->vx *= BOUND_DAMPING;
            p->x = EPS;
        }
        if (p->x + EPS > VIEW_WIDTH) {
            p->vx *= BOUND_DAMPING;
            p->x = VIEW_WIDTH - EPS;
        }
        if (p->y - EPS < 0.0) {
            p->vy *= BOUND_DAMPING;
            p->y = EPS;
        }
        if (p->y + EPS > VIEW_HEIGHT) {
            p->vy *= BOUND_DAMPING;
            p->y = VIEW_HEIGHT - EPS;
        }
    }
}


float avg_velocities( size_t start, size_t end, size_t step )
{
    double result = 0.0;
    for (size_t i = start; i < end; i += step) {
        /* the hypot(x,y) function is equivalent to sqrt(x*x +
           y*y); */
        result += hypot(particles[i].vx, particles[i].vy) / n_particles;
    }
    return result;
}


float update( void ) {
    double avg = 0.f;

    #pragma omp parallel default(none) shared(n_particles, particles) reduction(+:avg)
    {
        const size_t my_id = omp_get_thread_num();
        const size_t num_threads = omp_get_num_threads();
        const size_t my_start = (n_particles*my_id)/num_threads;
        const size_t my_end = (n_particles*(my_id+1))/num_threads;
        const size_t my_step = 1;

        compute_density_pressure(my_start, my_end, my_step, my_id);
        
        #pragma omp barrier
        compute_forces(my_start, my_end, my_step, my_id);

        #pragma omp barrier
        integrate(my_start, my_end, my_step);
        /* the average velocities MUST be computed at each step, even
        if it is not shown (to ensure constant workload per
        iteration) */

        avg = avg_velocities(my_start, my_end, my_step);
    }

    return avg;
}

int main(int argc, char **argv)
{
    srand(1234);

    particles = (particle_t*)malloc(MAX_PARTICLES * sizeof(*particles));
    assert( particles != NULL );

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

    const int max_th = omp_get_max_threads();

    near_rho = malloc(max_th * n * sizeof(float)); assert(near_rho != NULL);
    
    near_press_x = malloc(max_th * n * sizeof(float)); assert(near_press_x != NULL);

    near_press_y = malloc(max_th * n * sizeof(float)); assert(near_press_y != NULL);

    near_visc_x = malloc(max_th * n * sizeof(float)); assert(near_visc_x != NULL);

    near_visc_y = malloc(max_th * n * sizeof(float)); assert(near_visc_y != NULL);

    init_sph(n);
    double st = hpc_gettime();
    for (int s=0; s<nsteps; s++) {
        
        const float avg = update();

        if (s % 10 == 0)
            printf("step %5d, avgV=%f\n", s, avg);
    }
    printf("time elapsed: %f\n", hpc_gettime() - st);

    free(particles);
    return EXIT_SUCCESS;
}