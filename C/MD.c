/*
 *  Simple molecular dynamics code.
 *  $Id: MD-c.c,v 1.2 2002/01/31 16:43:14 spb Exp spb $
 *
 * This program implements:
 *     long range inverse square forces between particles. F = G * m1*m2 / r**2
 *     viscosity term     F = -u V
 * If 2 particles approach closer than Size we flip the direction of the
 * interaction force to approximate a collision.
 *
 * Coordinates are relative to a large central mass and the entire system is moving relative to the
 * viscous media.
 * If 2 particles approach closer than Size we flip the direction of the
 * interaction force to approximate a collision.
 *
 * This program was developed as part of a code optimisation course
 * and is therefore deliberately inefficient.
 *
 */
#include <math.h>
#include "coord.h"

void visc_force(int N,double *f, double *visc, double *vel);
void add_norm(int N,double *r, double *delta);
double force(double W, double delta, double r);
void wind_force(int N,double *f, double *visc, double vel);





void evolve(int count,double dt){
int  step;
int i,j,k,l;
int collided;

/*
 * Loop over timesteps.
 */
      for(step = 1;step<=count;step++){
        printf("timestep %d\n",step);
        printf("collisions %d\n",collisions);

/* set the viscosity term in the force calculation */
        for(j=0;j<Ndim;j++){
          visc_force(Nbody,f[j],visc,vel[j]);
        }
/* add the wind term in the force calculation */
        for(j=0;j<Ndim;j++){
          wind_force(Nbody,f[j],visc,wind[j]);
        }
/* calculate distance from central mass */
        for(k=0;k<Nbody;k++){
          r[k] = 0.0;
        }
        for(i=0;i<Ndim;i++){
	  add_norm(Nbody,r,pos[i]);
        }
        for(k=0;k<Nbody;k++){
          r[k] = sqrt(r[k]);
        }
       /* calculate central force */
        for(i=0;i<Nbody;i++){
	  for(l=0;l<Ndim;l++){
                f[l][i] = f[l][i] -
                   force(G*mass[i]*M_central,pos[l][i],r[i]);
	  }
	}
/* calculate pairwise separation of particles */
        k = 0;
        for(i=0;i<Nbody;i++){
          for(j=i+1;j<Nbody;j++){
            for(l=0;l<Ndim;l++){
              delta_pos[l][k] = pos[l][i] - pos[l][j];
            }
            k = k + 1;
          }
        }

/* calculate norm of seperation vector */
        for(k=0;k<Npair;k++){
          delta_r[k] = 0.0;
        }
        for(i=0;i<Ndim;i++){
	  add_norm(Npair,delta_r,delta_pos[i]);
        }
        for(k=0;k<Npair;k++){
          delta_r[k] = sqrt(delta_r[k]);
        }

/*
 * add pairwise forces.
 */
        k = 0;
        for(i=0;i<Nbody;i++){
          for(j=i+1;j<Nbody;j++){
            collided=0;

/*  flip force if close in */
              if( delta_r[k] >= Size ){
                for(l=0;l<Ndim;l++){
                f[l][i] = f[l][i] -
                   force(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
                f[l][j] = f[l][j] +
                   force(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
              }}else{
                  for(l=0;l<Ndim;l++){
                f[l][i] = f[l][i] +
                   force(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
                f[l][j] = f[l][j] -
                   force(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
                  }
		collided=1;
              }

	    if( collided == 1 ){
	      collisions++;
	    }
            k = k + 1;
          }
        }

/* update positions */
        for(i=0;i<Nbody;i++){
          for(j=0;j<Ndim;j++){
            pos[j][i] = pos[j][i] + dt * vel[j][i];
          }
        }

/* update velocities */
        for(i=0;i<Nbody;i++){
          for(j=0;j<Ndim;j++){
            vel[j][i] = vel[j][i] + dt * (f[j][i]/mass[i]);
          }
        }


      }

}




