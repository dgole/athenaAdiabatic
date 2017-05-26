#include "copyright.h"
/*============================================================================*/
/*! \file blast.c
 *  \brief Problem generator for spherical blast wave problem.
 *
 * PURPOSE: Problem generator for spherical blast wave problem.
 *
 * REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for
 *   multidimensional MHD", ApJ, 530, 508 (2000), and references therein.     */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  Prim1DS W;
  Cons1DS U1d;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  int vdir, bdep;
  Real pressure,bpressure,drat,prat,rad,pa,da,h,x1,x2,x3,Q_AD,v0;
  Real b0x, b0y, b0z, Bx=0.0, By=0.0, Bz=0.0, x1p, x2p, x3p, roty, rotz;
  double theta;
  da  = par_getd("problem","damb");
  pa  = par_getd("problem","pamb");
  h  = par_getd("problem","h");
#ifdef MHD
  b0x = par_getd("problem","b0x");
  b0y = par_getd("problem","b0y");
  b0z = par_getd("problem","b0z");
  v0 = par_getd("problem","v0");
  vdir = par_geti("problem","vdir");
  bdep = par_geti("problem","bdep");
  if (bdep==-1) {
	roty = (3.1415926535/180.0)*par_getd("problem","roty");
	rotz = (3.1415926535/180.0)*par_getd("problem","rotz");
  }
		
#endif
	
  if (vdir==0) {
    W.Vx = v0;
    W.Vy = 0.0;
    W.Vz = 0.0;
  }
  if (vdir==1) {
    W.Vx = 0.0;
    W.Vy = v0;
    W.Vz = 0.0;
  }
  if (vdir==2) {
    W.Vx = 0.0;
    W.Vy = 0.0;
    W.Vz = v0;
  }

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
		cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	
		x1p=x1*cos(rotz)-x2*sin(rotz);
		x2p=x1*sin(rotz)+x2*cos(rotz);
		x3p=x3;

        if (bdep==-1) {
          Bx=b0x*exp(-(x1p*x1p)/(2.0*h*h));  
		  W.By=b0y*exp(-(x1p*x1p)/(2.0*h*h));
		  W.Bz=b0z*exp(-(x1p*x1p)/(2.0*h*h));		
		}
        if (bdep==0) {
          Bx=b0x*exp(-(x1p*x1p)/(2.0*h*h));  
		  W.By=b0y*exp(-(x1p*x1p)/(2.0*h*h));
		  W.Bz=b0z*exp(-(x1p*x1p)/(2.0*h*h));		
		}
        if (bdep==1) {
          Bx=b0x*exp(-(x2p*x2p)/(2.0*h*h));  
		  W.By=b0y*exp(-(x2p*x2p)/(2.0*h*h));
		  W.Bz=b0z*exp(-(x2p*x2p)/(2.0*h*h));		
		}
        if (bdep==2) {
          Bx=b0x*exp(-(x3p*x3p)/(2.0*h*h));  
		  W.By=b0y*exp(-(x3p*x3p)/(2.0*h*h));
		  W.Bz=b0z*exp(-(x3p*x3p)/(2.0*h*h));		
		}

 
	
#ifndef ISOTHERMAL
		bpressure=0.5*(Bx*Bx+W.By*W.By+W.Bz*W.Bz);
        W.P = pa - bpressure/0.666666666667;
#endif
        W.d = da;

	
		U1d = Prim1D_to_Cons1D(&(W),&Bx);
	
		pGrid->U[k][j][i].d  = U1d.d;
		pGrid->U[k][j][i].M1 = U1d.Mx;
		pGrid->U[k][j][i].M2 = U1d.My;
		pGrid->U[k][j][i].M3 = U1d.Mz;
#ifndef ISOTHERMAL
		pGrid->U[k][j][i].E  = U1d.E;
#endif
#ifdef MHD
		pGrid->B1i[k][j][i] = Bx;
		pGrid->B2i[k][j][i] = U1d.By;
		pGrid->B3i[k][j][i] = U1d.Bz;
		pGrid->U[k][j][i].B1c = Bx;
		pGrid->U[k][j][i].B2c = U1d.By;
		pGrid->U[k][j][i].B3c = U1d.Bz;
		if (i == ie && ie > is) pGrid->B1i[k][j][i+1] = Bx;
		if (j == je && je > js) pGrid->B2i[k][j+1][i] = U1d.By;
		if (k == ke && ke > ks) pGrid->B3i[k+1][j][i] = U1d.Bz;
#endif /* MHD */
      }
    }
  }
#ifdef RESISTIVITY 
  eta_Ohm = 0.0;
  Q_AD    = par_getd_def("problem","Q_AD",0.0);
  Q_Hall  = 0.0;
  d_ind   = 0.0;
#endif


}



/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

#ifdef RESISTIVITY
void get_eta_user(GridS *pG, int i, int j, int k,
                             Real *eta_O, Real *eta_H, Real *eta_A)
{

  *eta_O = 0.0;
  *eta_H = 0.0;
  *eta_A = 1.0;

  return;
}
#endif


void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
}

void Userwork_before_loop(MeshS *pM)
{
}
