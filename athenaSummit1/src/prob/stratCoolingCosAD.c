#include "copyright.h"
/*============================================================================*/
/*! \file strat.c
 *  \brief Problem generator for stratified 3D shearing sheet.
 *
 * PURPOSE:  Problem generator for stratified 3D shearing sheet.  Based on the 
 *   initial conditions described in "Three-dimensional Magnetohydrodynamic
 *   Simulations of Vertically Stratified Accretion Disks" by Stone, Hawley,
 *   Gammie & Balbus.
 *
 * Several different field configurations and perturbations are possible:
 *  ifield = 1 - Bz=B0 sin(x1) field with zero-net-flux [default]
 *  ifield = 2 - uniform Bz
 *  ifield = 3 - B=(0,B0cos(kx*x1),B0sin(kx*x1))= zero-net flux w helicity
 *  ifield = 4 - uniform By, but only for |z|<2
 *  ifield = 5 - By with constant \beta versus z
 *  ifield = 6 - flux tube at disk center from Hirose et al.
 *  ifield = 7 - zero field everywhere
 *  ifield = 8 - uniform Bz with sin(x1) field added on top
 *
 * - ipert = 1 - random perturbations to P and V [default, used by HGB]
 *
 * Code must be configured using --enable-shearing-box
 *
 * REFERENCE: Stone, J., Hawley, J. & Balbus, S. A., ApJ 463, 656-673 (1996)
 *            Hawley, J. F. & Balbus, S. A., ApJ 400, 595-609 (1992)	      */
/*============================================================================*/

#include <float.h>
#include <math.h>

#include <time.h>

#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()           - random number generator from NR
 * UnstratifiedDisk() - tidal potential in 3D shearing box
 * VertGrav()         - potential for vertical component of gravity
 * expr_*()         - computes new output variables
 * hst_*            - adds new history variables
 * strat_ix3        - vertical outflow boundary for bottom of grid
 * strat_ox3        - vertical outflow boundary for top of grid
 * output_1d()      - dumps horizontally averaged quantities to a text file
 *============================================================================*/

static double ran2(long int *idum);
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3);
static Real VertGrav(const Real x1, const Real x2, const Real x3);
static void strat_ix3(GridS *pG);
static void strat_ox3(GridS *pG);
static Real expr_dV2(const GridS *pG, const int i, const int j, const int k);
static Real expr_beta(const GridS *pG, const int i, const int j, const int k);
static Real expr_ME(const GridS *pG, const int i, const int j, const int k);
static Real expr_KE(const GridS *pG, const int i, const int j, const int k);
static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j,const int k);
static Real hst_rho_dVy2(const GridS *pG, const int i, const int j, const int k);
static Real mass_cons(MeshS *pM);
#ifdef ADIABATIC
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k);
#endif
#ifdef MHD
static Real hst_Bx(const GridS *pG, const int i, const int j, const int k);
static Real hst_By(const GridS *pG, const int i, const int j, const int k);
static Real hst_Bz(const GridS *pG, const int i, const int j, const int k);
static Real hst_BxBy(const GridS *pG, const int i, const int j, const int k);
#endif /* MHD */
static void combine_densities(MeshS *pM);
static void output_1d(MeshS *pM, OutputS *pOut);

/* top and bottom of root Domain, shared with outputs, etc. */
static Real ztop, zbtm;

/* Apply a density floor - useful for large |z| regions */
static Real D_FLOOR = 1.e-5;

/* Apply a pressure floor - neccessary if adiabatic and applying a density floor */
#ifdef ADIABATIC
static Real P_FLOOR = 5.e-6;
#endif /* ADIABATIC */

/* Apply a temperature floor */
#ifdef ADIABATIC
static Real T_FLOOR = 0.0;
#endif /* ADIABATIC */

static Real R_AU = 100.0;
static Real Sig_FUV=0.001;
static Real ionfrac_FUV=1.e-5;
static Real Am0=1.0;

/* Total mass in grid */
static Real mass;

/* Flag to determine whether or not to employ outflow boundaries
   in z */
static int zbc_out = 1;

#ifdef RESISTIVITY
static void density_profile(MeshS *pM);
static void ionization_rate(MeshS *pM);
static Real get_Am_FUV(const GridS *pG, const int i, const int j, const int k, Real Am0);
static Real expr_Lambda(const GridS *pG, const int i, const int j, const int k);
static Real expr_chi(const GridS *pG, const int i, const int j, const int k);
static Real expr_Am(const GridS *pG, const int i, const int j, const int k);
#endif

#ifdef RESISTIVITY
/* Horizontally averaged gas density and column densities */
static Real *rho_avg;
static Real *Sigma_top, *Sigma_bot;
static Real *ion_rate;    /* ionization rate   */
/* Smoothing length to connect low Am region to high Am region */
static Real delz = 0.05;
Real zit,zib; /* z location for top and bottom ionization front */
/* lookup table of magnetic diffusivities */
static int  nzeta, nrho;
static Real *rho_tab,*zeta_tab,dlogrho1,dlogzeta1;
static Real **etaO_tab,**Q_H_tab,**Q_A_tab;
static Real myrho,myB,myOmega,mySigma=1.0;
#endif
static double *densInit1;
static double *tempInit1;
static double *densInit2;
static double *tempInit2;
static double *densInit3;
static double *tempInit3;
static double *gdens1D;
static double ***ldens,***gdens;
static int kbase1,kbase2;
static int nxgrid,nygrid,nzgrid;
static int gridnumx1,gridnumy1,gridnumz1;
static int count1 = 0;
static Real lastRecalcTime = 0.0;
static Real dtRecalc = 0.05;



/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;

  FILE *fp;
  Real xFP[160],dFP[160],vxFP[160],vyFP[160];
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int ixs,jxs,kxs,i,j,k,ipert,ifield,icool;
  long int iseed = -1; /* Initialize on the first call to ran2 */
  Real x1,x2,x3,xmin,xmax,Lx,Ly,Lz,divb;
  Real den=1.0, pres=5.0e-7, rd, rp, rvx, rvy, rvz, T0n, T0nm1;
  Real beta,B0,kx,amp,mass_cell;
  int nwx,nwy,nwz;  /* input number of waves per Lx,Ly,Lz [default=1] */
  double rval;
  static int frst=1;  /* flag so new history variables enrolled only once */

  if (pGrid->Nx[1] == 1){
    ath_error("[problem]: Strat only works on a 2D or 3D grid\n");
  }

#ifdef STATIC_MESH_REFINEMENT
  ath_error("[problem]: Strat does not yet work with SMR\n");
#endif

/* Reset d_MIN to be 0.1 of D_FLOOR */
  D_FLOOR = par_getd("problem","dfloor");
#ifdef adiabatic  
  P_FLOOR = par_getd("problem","pfloor");
  T_FLOOR = par_getd("problem","tfloor");
#endif //adiabatic
  d_MIN = 0.1*D_FLOOR;

/* Read problem parameters.  Note Omega set to 10^{-3} by default */
#ifdef ISOTHERMAL
  pres=den*Iso_csound2;
#endif
  Omega_0 = par_getd_def("problem","omega",1.0e-3);
  qshear  = par_getd_def("problem","qshear",1.5);
  amp = par_getd("problem","amp");
  beta = par_getd("problem","beta");
  pres = par_getd("problem","pres");
  B0 = sqrt((double)(2.0*pres/beta));
  ifield = par_geti_def("problem","ifield", 1);
  ipert = par_geti_def("problem","ipert", 1);
  R_AU = par_getd("problem","R_AU");
  Am0 = par_getd("problem","Am0");
  Sig_FUV = par_getd("problem","Sig_FUV");
  ionfrac_FUV = par_getd("problem","ionfrac_FUV");

#ifdef ADIABATIC
  icool = par_geti_def("problem","icool", 0);
  if (icool==2) {
	CoolingFunc = CosCoolingFunc;
	tauCool= par_getd("problem","tauCool");
	zq= par_getd("problem","zq");
	Tmid= par_getd("problem","Tmid");
	Tatm= par_getd("problem","Tatm");
	delta= par_getd("problem","delta");
	den=pres/Tmid;
  }
  else if (icool==1) {
	CoolingFunc = FlatCoolingFunc;
	T0 = par_getd("problem","T0");
	tauCool = par_getd("problem","tauCool");
  }
  else {
	CoolingFunc = NoCoolingFunc;
  }
#endif /*ADIABATIC*/

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  kxs = pGrid->Disp[2];
  iseed = -1 - (ixs + pDomain->Nx[0]*(jxs + pDomain->Nx[1]*kxs));

/* Initialize boxsize */
  ztop = pDomain->RootMaxX[2];
  zbtm = pDomain->RootMinX[2];
  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  Lz = pDomain->RootMaxX[2] - pDomain->RootMinX[2];

/* initialize wavenumbers, given input number of waves per L */
  nwx = par_geti_def("problem","nwx",1);
  kx = (2.0*PI/Lx)*((double)nwx);








// NEW STUFF TO GET INITIAL DENSITY AND PRESSURE PROFILES RIGHT WITH TEMPERATURE PROFILE
  densInit1 = (double*)calloc_1d_array(pDomain->NGrid[2]*pGrid->Nx[2],sizeof(double));
  tempInit1 = (double*)calloc_1d_array(pDomain->NGrid[2]*pGrid->Nx[2],sizeof(double));
  densInit2 = (double*)calloc_1d_array(pDomain->NGrid[2]*pGrid->Nx[2],sizeof(double));
  tempInit2 = (double*)calloc_1d_array(pDomain->NGrid[2]*pGrid->Nx[2],sizeof(double));
  densInit3 = (double*)calloc_1d_array(pDomain->NGrid[2]*pGrid->Nx[2],sizeof(double));
  tempInit3 = (double*)calloc_1d_array(pDomain->NGrid[2]*pGrid->Nx[2],sizeof(double));
  // algorithm needs starting points, this is at z=0 (although there is no zone there)
  densInit1[0] = den;
  tempInit1[0] = Tmid;
  // calculate entire T0 profile, regardless of actual mpi domain we are currently in
  for (k=1; k<=pDomain->NGrid[2]*pGrid->Nx[2]/2; k++) {
		x3 = k*pGrid->dx3;
    if (x3<zq) {tempInit1[k]=Tatm+(Tmid-Tatm)*pow(cos(PI*x3/(2.0*zq)),2.0*delta);}
		else {tempInit1[k]=Tatm;}
  }
  // calculate entire density profile, regardless of actual mpi domain we are currently in
  for (k=1; k<=pDomain->NGrid[2]*pGrid->Nx[2]/2; k++) {
		x3 = k*pGrid->dx3;
		densInit1[k] = densInit1[k-1]*(2.0-((x3-pGrid->dx3)*pGrid->dx3/tempInit1[k-1])-(tempInit1[k]/tempInit1[k-1]));
  }
  // interpolate half a cell over and re assign
  for (k=0; k<(pDomain->NGrid[2]*pGrid->Nx[2]/2); k++) {
		densInit1[k]=(densInit1[k]+densInit1[k+1])/2.0;		
  }
  // mirror array for -z parts of the box
  for (k=0; k<(pDomain->NGrid[2]*pGrid->Nx[2]/2); k++) {
		int ktemp=(pDomain->NGrid[2]*pGrid->Nx[2]/2)-k-1;
 		tempInit2[ktemp]=tempInit1[k];
		densInit2[ktemp]=densInit1[k];
  }
  // combine arrays 
  for (k=0; k<(pDomain->NGrid[2]*pGrid->Nx[2]); k++) {
		if (k<(pDomain->NGrid[2]*pGrid->Nx[2]/2)) {	
	  	tempInit3[k]=tempInit2[k];
	 		densInit3[k]=MAX(densInit2[k],D_FLOOR);
		}
		else {
	  	tempInit3[k]=tempInit1[k-(pDomain->NGrid[2]*pGrid->Nx[2]/2)];
	 		densInit3[k]=MAX(densInit1[k-(pDomain->NGrid[2]*pGrid->Nx[2]/2)],D_FLOOR);
		}
  }
  // printer loop
  for (k=0; k<(pDomain->NGrid[2]*pGrid->Nx[2]); k++) {
		if (myID_Comm_world==0){
	  	printf("%s %i %i %i %i %i %0.9G \n", "a", myID_Comm_world, ks, ke, pGrid->Disp[2], k, densInit3[k]);    
  	}
  }


  for (k=ks; k<=ke; k++) {
 	 for (j=js; j<=je; j++) {
 	   for (i=is; i<=ie; i++) {
  	   cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	  
/* Initialize perturbations
 *  ipert = 1 - random perturbations to P and V [default, used by HGB]
 */
      if (ipert == 1) {
        rval = amp*(ran2(&iseed) - 0.5);
		rd =densInit3[k+pGrid->Disp[2]-nghost]; 
		rd*=(1.0+2.0*rval); 
        if (rd < D_FLOOR) {rd = D_FLOOR;}
		#ifdef ADIABATIC
		  rp = rd*tempInit3[k+pGrid->Disp[2]-nghost];  
		  if (rp < P_FLOOR) {rp = P_FLOOR;}
		#endif
		if (i==is && j==js){
		  printf("%s %i %i %i %i %i %i %0.9G %0.9G \n", "b", myID_Comm_world, ks, ke, pGrid->Disp[2], k, k+pGrid->Disp[2]-nghost, x3, rd);
		}

/* To conform to HGB, the perturbations to V/Cs are (1/5)amp/sqrt(Gamma)  */
        rval = amp*(ran2(&iseed) - 0.5);
        rvx = 0.4*rval*sqrt(pres/den);

        rval = amp*(ran2(&iseed) - 0.5);
        rvy = 0.4*rval*sqrt(pres/den);

        rval = amp*(ran2(&iseed) - 0.5);
        rvz = 0.4*rval*sqrt(pres/den);
      }

/* Initialize d, M, and P.  For 3D shearing box M1=Vx, M2=Vy, M3=Vz
 * With FARGO do not initialize the background shear */ 

      pGrid->U[k][j][i].d  = rd;
      pGrid->U[k][j][i].M1 = rd*rvx;
      pGrid->U[k][j][i].M2 = rd*rvy;
#ifndef FARGO
      pGrid->U[k][j][i].M2 -= rd*(qshear*Omega_0*x1);
#endif
      pGrid->U[k][j][i].M3 = rd*rvz;
#ifdef ADIABATIC
      pGrid->U[k][j][i].E = rp/Gamma_1
        + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2) 
             + SQR(pGrid->U[k][j][i].M3))/rd;
#endif

/* Initialize magnetic field.  For 3D shearing box B1=Bx, B2=By, B3=Bz
 *  ifield = 1 - Bz=B0 sin(x1) field with zero-net-flux [default]
 *  ifield = 2 - uniform Bz
 *  ifield = 3 - B=(0,B0cos(kx*x1),B0sin(kx*x1))= zero-net flux w helicity
 *  ifield = 4 - uniform By, but only for |z|<2
 *  ifield = 5 - By with constant \beta versus z
 *  ifield = 6 - flux tube at disk center from Hirose et al.
 *  ifield = 7 - zero field everywhere
 *  ifield = 8 - uniform Bz with sin(x1) field added on top
 */

#ifdef MHD
      pGrid->B1i[k][j][i] = 0.0;
      pGrid->B2i[k][j][i] = 0.0;
      pGrid->B3i[k][j][i] = 0.0;
      if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
      if (j==je) pGrid->B2i[k][je+1][i] = 0.0;
      if (k==ke) pGrid->B3i[ke+1][j][i] = 0.0;

      if (ifield == 1) {
        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = 0.0;
        pGrid->B3i[k][j][i] = B0*(sin((double)kx*x1));
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = 0.0;
        if (k==ke) pGrid->B3i[ke+1][j][i] = B0*(sin((double)kx*x1));
      }
      if (ifield == 2) {
        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = 0.0;
        pGrid->B3i[k][j][i] = B0;
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = 0.0;
        if (k==ke) pGrid->B3i[ke+1][j][i] = B0;
      }
      if (ifield == 3) {
        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = B0*(cos((double)kx*x1));
        pGrid->B3i[k][j][i] = B0*(sin((double)kx*x1));
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = B0*(cos((double)kx*x1));
        if (k==ke) pGrid->B3i[ke+1][j][i] = B0*(sin((double)kx*x1));
      }
      if (ifield == 4 && fabs(x3) < 2.0) {
        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = B0;
        pGrid->B3i[k][j][i] = 0.0;
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = B0;
        if (k==ke) pGrid->B3i[ke+1][j][i] = 0.0;
      }
      if (ifield == 5) {
        /* net toroidal field with constant \beta with height */
        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = sqrt(den*exp(-x3*x3)*SQR(Omega_0)/beta);
        pGrid->B3i[k][j][i] = 0.0;
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = sqrt(den*exp(-x3*x3)*SQR(Omega_0)/beta);
        if (k==ke) pGrid->B3i[ke+1][j][i] = 0.0;
      }
      if (ifield == 7) {
        /* zero field everywhere */
        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = 0.0;
        pGrid->B3i[k][j][i] = 0.0;
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = 0.0;
        if (k==ke) pGrid->B3i[ke+1][j][i] = 0.0;
      }
      if (ifield == 8) {
        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = 0.0;
        pGrid->B3i[k][j][i] = B0*(1.0+0.25*sin((double)kx*x1));
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = 0.0;
        if (k==ke) pGrid->B3i[ke+1][j][i] = B0*(1.0+0.25*sin((double)kx*x1));
      }
#endif /* MHD */
    }
  }}

#ifdef MHD

  if (ifield == 6) {
  /* flux tube of Hirose et al. We put this down here to break away from the
     for loops used above. */
    Real xc=0.0,zc=0.0,Bpratio=0.25,Bp0,rad0,rad,Ay0,x1h,x3h;
    Real Ay[pGrid->Nx[2]+2*nghost][pGrid->Nx[1]+2*nghost][pGrid->Nx[0]+2*nghost];
    for (k=ks; k<=ke+2; k++) {
      for (j=js; j<=je+2; j++) {
        for (i=is; i<=ie+2; i++) {
          cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
          x1h = x1-0.5*pGrid->dx1;
          x3h = x3-0.5*pGrid->dx3;
          Bp0 = B0*Bpratio;
  /* We are assuming that the scale height is H = 1 as usual */
          rad0 = Lx/4.0;
          Ay0 = Bp0*rad0/PI;
          rad = sqrt(SQR(x1h-xc)+SQR(x3h-zc));
  /* Calculate vector potential */
          if (rad < rad0) {
            Ay[k][j][i]=-Ay0*(1.0+cos(rad/rad0*PI));
          } else {
            Ay[k][j][i]=0.0;
          }
        }
      }
    }
  /* In this case, we calculate face fields from vect. potential */
    for (k=ks; k<=ke+1; k++) {
      for (j=js; j<=je+1; j++) {
        for (i=is; i<=ie+1; i++) {
          pGrid->B1i[k][j][i] = -(Ay[k+1][j][i]-Ay[k][j][i])/pGrid->dx3;
          pGrid->B3i[k][j][i] =  (Ay[k][j][i+1]-Ay[k][j][i])/pGrid->dx1;
        }
      }
    }
  /* Sync poloidal centered fields to face fields for the next step */
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je+1; j++) {
        for (i=is; i<=ie; i++) {
          pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i]+pGrid->B1i[k][j][i+1]);
          pGrid->U[k][j][i].B3c = 0.5*(pGrid->B3i[k][j][i]+pGrid->B3i[k+1][j][i]);
        }
      }
    }
  /* If mag. of total poloidal field (defined at cell center) is zero, then By is zero */
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je+1; j++) {
        for (i=is; i<=ie; i++) {
          if (SQR(pGrid->U[k][j][i].B1c)+SQR(pGrid->U[k][j][i].B3c) != 0.0) {
            pGrid->B2i[k][j][i] = sqrt(SQR(B0)-(SQR(pGrid->U[k][j][i].B1c)+SQR(pGrid->U[k][j][i].B3c)));
          } else {
            pGrid->B2i[k][j][i] = 0.0;
          }
        }
      }
  /* Finally, calculate cell centered By field from face fields */
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i]+pGrid->B2i[k][j+1][i]);
#ifdef ADIABATIC
          pGrid->U[k][j][i].E += 0.5*(SQR(pGrid->U[k][j][i].B1c)
               + SQR(pGrid->U[k][j][i].B2c) + SQR(pGrid->U[k][j][i].B3c));
#endif
        }
      }
    }
  } else {
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i]+pGrid->B1i[k][j][i+1]);
          pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i]+pGrid->B2i[k][j+1][i]);
          pGrid->U[k][j][i].B3c = 0.5*(pGrid->B3i[k][j][i]+pGrid->B3i[k+1][j][i]);
#ifdef ADIABATIC
          pGrid->U[k][j][i].E += 0.5*(SQR(pGrid->U[k][j][i].B1c)
               + SQR(pGrid->U[k][j][i].B2c) + SQR(pGrid->U[k][j][i].B3c));
#endif
        }
      }
    }
  } /* End if ifield = 6 */

/* Finally, let's check that the field is divergenceless */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        divb = (pGrid->B1i[k][j][i+1]-pGrid->B1i[k][j][i])/pGrid->dx1 +
               (pGrid->B2i[k][j+1][i]-pGrid->B2i[k][j][i])/pGrid->dx2 +
               (pGrid->B3i[k+1][j][i]-pGrid->B3i[k][j][i])/pGrid->dx3;
        if (fabs(divb) >= 1.e-12) {
          ath_error("[problem]: Nonzero divergence of initial magnetic field\n");
        }
      }
    }
  }

#endif /* MHD */

/* With viscosity and/or resistivity, read eta_Ohm and nu_V */
#ifdef RESISTIVITY
  eta_Ohm = par_getd_def("problem","eta_O",0.0);
  Q_Hall  = par_getd_def("problem","Q_H",0.0);
//  Q_AD    = par_getd_def("problem","Q_A",0.0);
  d_ind   = par_getd_def("problem","d_ind",0.0);
#endif
#ifdef VISCOSITY
  nu_iso = par_getd_def("problem","nu_iso",0.0);
  nu_aniso = par_getd_def("problem","nu_aniso",0.0);
#endif

/* enroll gravitational potential function */

  StaticGravPot = VertGrav;
  ShearingBoxPot = UnstratifiedDisk;

/* Enroll vertically stratified outflow boundaries */
  if (zbc_out == 1) {
    bvals_mhd_fun(pDomain, left_x3, strat_ix3);
    bvals_mhd_fun(pDomain, right_x3, strat_ox3);
  }

/* enroll new history variables */

  if (frst == 1) {
    dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
    dump_history_enroll(hst_rho_dVy2, "<rho dVy^2>");
#ifdef ADIABATIC
    dump_history_enroll(hst_E_total, "<E + rho Phi>");
#endif
#ifdef MHD
    dump_history_enroll(hst_Bx, "<Bx>");
    dump_history_enroll(hst_By, "<By>");
    dump_history_enroll(hst_Bz, "<Bz>");
    dump_history_enroll(hst_BxBy, "<-Bx By>");
#endif /* MHD */
    frst = 0;

  }

/*  mass_cell = 0.0;
  for (k=0; k<=pGrid->Nx[2]-1; k++) {
	  for (j=0; j<=pGrid->Nx[1]-1; j++) {
			for (i=0; i<=pGrid->Nx[0]-1; i++) {
		  	mass_cell += pGrid->U[k][j][i].d;
			}
	  }
  }
	mass_cell *= pGrid->dx1*pGrid->dx2*pGrid->dx3;
  mass = mass_cell;
	printf("%s %0.9G \n", "total mass: ", mass);
*/
/*
#ifdef MPI_PARALLEL
    ierr = MPI_Reduce(&(mtot), &(my_mtot), 1,
           MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
#else
    my_mtot = mtot;
#endif
*/
//    if(myID_Comm_world==0){  /* I'm the parent */
//      my_mtot /= dVol;
		
  mass_cell = 0.0;
  for (k=0; k<=pDomain->Nx[2]-1; k++) {
    x3 = pDomain->MinX[2] + ((Real)k + 0.5)*pGrid->dx3;
    mass_cell += den*exp(-x3*x3)*pGrid->dx3;
  }
  mass = mass_cell/Lz;
	printf("%s %0.9G \n", "total mass in initialization from old routine: ", mass);

  mass_cell = 0.0;
  for (k=0; k<=pDomain->Nx[2]-1; k++) {
		mass_cell += densInit3[k]*pGrid->dx3;
  }
  mass = mass_cell/Lz;
	printf("%s %0.9G \n", "total mass in initialization from new routine: ", mass);

  return;
}

/*==============================================================================
 * PUBLIC PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_before_loop    - problem specific work BEFORE main loop
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

/*
 * 'problem_read_restart' must enroll gravity on restarts
 */

void problem_read_restart(MeshS *pM, FILE *fp)
{

  GridS *pGrid;
  DomainS *pDomain;
  
	pDomain = (DomainS*)&(pM->Domain[0][0]);
  pGrid = pM->Domain[0][0].Grid;

	Real Lz;
	Lz = pDomain->RootMaxX[2] - pDomain->RootMinX[2];
  int icool;
	int i, j, k;
	Real x1, x2, x3;
	Real mass_cell;
	Real pres, den; 





/* Reset d_MIN to be 0.1 of D_FLOOR */
  D_FLOOR = par_getd("problem","dfloor");
#ifdef adiabatic
  P_FLOOR = par_getd("problem","pfloor");
  T_FLOOR = par_getd("problem","tfloor");
#endif
  d_MIN = 0.1*D_FLOOR;

  Omega_0 = par_getd_def("problem","omega",1.0e-3);
  qshear  = par_getd_def("problem","qshear",1.5);
  R_AU = par_getd("problem","R_AU");
  Am0 = par_getd("problem","Am0");
  Sig_FUV = par_getd("problem","Sig_FUV");
  ionfrac_FUV = par_getd("problem","ionfrac_FUV");


#ifdef ADIABATIC
  icool = par_geti_def("problem","icool", 0);
  if (icool==2) {
	CoolingFunc = CosCoolingFunc;
	tauCool= par_getd("problem","tauCool");
	zq= par_getd("problem","zq");
	Tmid= par_getd("problem","Tmid");
	Tatm= par_getd("problem","Tatm");
	delta= par_getd("problem","delta");
  }
  else if (icool==1) {
	CoolingFunc = FlatCoolingFunc;
	T0 = par_getd("problem","T0");
	tauCool = par_getd("problem","tauCool");
  }
  else {
	CoolingFunc = NoCoolingFunc;
  }
#endif /*ADIABATIC*/

  pres = par_getd("problem","pres");
	den=pres/Tmid;

/* With viscosity and/or resistivity, read eta_Ohm and nu_V */
#ifdef RESISTIVITY
  eta_Ohm = par_getd_def("problem","eta_O",0.0);
  Q_Hall  = par_getd_def("problem","Q_H",0.0);
//  Q_AD    = par_getd_def("problem","Q_A",0.0);
  d_ind   = par_getd_def("problem","d_ind",0.0);
#endif
#ifdef VISCOSITY
  nu_iso = par_getd_def("problem","nu_iso",0.0);
  nu_aniso = par_getd_def("problem","nu_aniso",0.0);
#endif

/* enroll gravitational potential function */

  StaticGravPot = VertGrav;
  ShearingBoxPot = UnstratifiedDisk;

/* enroll new history variables */

  dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
  dump_history_enroll(hst_rho_dVy2, "<rho dVy^2>");
#ifdef ADIABATIC
  dump_history_enroll(hst_E_total, "<E + rho Phi>");
#endif
#ifdef MHD
  dump_history_enroll(hst_Bx, "<Bx>");
  dump_history_enroll(hst_By, "<By>");
  dump_history_enroll(hst_Bz, "<Bz>");
  dump_history_enroll(hst_BxBy, "<-Bx By>");
#endif /* MHD */

/* If using outflow boundaries, have to enroll them here too */

  if (zbc_out == 1) {

    int nl, nd;
    ztop = pM->RootMaxX[2];
    zbtm = pM->RootMinX[2];

    for (nl=0; nl<(pM->NLevels); nl++){
      for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
        if (pM->Domain[nl][nd].Grid != NULL){

          bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x3, strat_ix3);
          bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x3, strat_ox3);

        }
      }
    }

  }

// NEW STUFF TO GET INITIAL DENSITY AND PRESSURE PROFILES RIGHT WITH TEMPERATURE PROFILE
  densInit1 = (double*)calloc_1d_array(pDomain->NGrid[2]*pGrid->Nx[2],sizeof(double));
  tempInit1 = (double*)calloc_1d_array(pDomain->NGrid[2]*pGrid->Nx[2],sizeof(double));
  densInit2 = (double*)calloc_1d_array(pDomain->NGrid[2]*pGrid->Nx[2],sizeof(double));
  tempInit2 = (double*)calloc_1d_array(pDomain->NGrid[2]*pGrid->Nx[2],sizeof(double));
  densInit3 = (double*)calloc_1d_array(pDomain->NGrid[2]*pGrid->Nx[2],sizeof(double));
  tempInit3 = (double*)calloc_1d_array(pDomain->NGrid[2]*pGrid->Nx[2],sizeof(double));
  // algorithm needs starting points, this is at z=0 (although there is no zone there)
  densInit1[0] = den;
  tempInit1[0] = Tmid;
  // calculate entire T0 profile, regardless of actual mpi domain we are currently in
  for (k=1; k<=pDomain->NGrid[2]*pGrid->Nx[2]/2; k++) {
	x3 = k*pGrid->dx3;
    if (x3<zq) {tempInit1[k]=Tatm+(Tmid-Tatm)*pow(cos(PI*x3/(2.0*zq)),2.0*delta);}
	else {tempInit1[k]=Tatm;}
  }
  // calculate entire density profile, regardless of actual mpi domain we are currently in
  for (k=1; k<=pDomain->NGrid[2]*pGrid->Nx[2]/2; k++) {
	x3 = k*pGrid->dx3;
	densInit1[k] = densInit1[k-1]*(2.0-((x3-pGrid->dx3)*pGrid->dx3/tempInit1[k-1])-(tempInit1[k]/tempInit1[k-1]));
  }
  // interpolate half a cell over and re assign
  for (k=0; k<(pDomain->NGrid[2]*pGrid->Nx[2]/2); k++) {
	densInit1[k]=(densInit1[k]+densInit1[k+1])/2.0;		
  }
  // mirror array for -z parts of the box
  for (k=0; k<(pDomain->NGrid[2]*pGrid->Nx[2]/2); k++) {
	int ktemp=(pDomain->NGrid[2]*pGrid->Nx[2]/2)-k-1;
 	tempInit2[ktemp]=tempInit1[k];
	densInit2[ktemp]=densInit1[k];
  }
  // combine arrays 
  for (k=0; k<(pDomain->NGrid[2]*pGrid->Nx[2]); k++) {
	if (k<(pDomain->NGrid[2]*pGrid->Nx[2]/2)) {	
	  tempInit3[k]=tempInit2[k];
	  densInit3[k]=MAX(densInit2[k],D_FLOOR);
	}
	else {
	  tempInit3[k]=tempInit1[k-(pDomain->NGrid[2]*pGrid->Nx[2]/2)];
	  densInit3[k]=MAX(densInit1[k-(pDomain->NGrid[2]*pGrid->Nx[2]/2)],D_FLOOR);
	}
  }
  // printer loop
  for (k=0; k<(pDomain->NGrid[2]*pGrid->Nx[2]); k++) {
	if (myID_Comm_world==0){
	  //printf("%i %i %i %i %i %0.9G \n", myID_Comm_world, ks, ke, pGrid->Disp[2], k, densInit3[k]);    
  	}
  }

  mass_cell = 0.0;
  for (k=0; k<=pDomain->Nx[2]-1; k++) {
		mass_cell += densInit3[k]*pGrid->dx3;
  }
  mass = mass_cell/Lz;
	printf("%s %0.9G \n", "total mass in restart from new routine: ", mass);


  return;
}

/* Get_user_expression computes dVy */
ConsFun_t get_usr_expr(const char *expr)
{
  if(strcmp(expr,"dVy")==0) return expr_dV2;
  else if(strcmp(expr,"beta")==0) return expr_beta;
  else if(strcmp(expr,"ME")==0) return expr_ME;
  else if(strcmp(expr,"KE")==0) return expr_KE;
  else if(strcmp(expr,"BxBy")==0) return hst_BxBy;
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name)
{
  if(strcmp(name,"1d")==0) return output_1d;
  return NULL;
}

#ifdef RESISTIVITY
void get_eta_user(GridS *pG, int i, int j, int k,
                            Real *eta_O, Real *eta_H, Real *eta_A)
{
  int k0,rj,zi;
  Real x3;
  Real zeta,rho,Bmag,Bsq,eta_A0;
  Real Q_H,Q_A,Am,beta;

  Bsq  = SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c);
  Am = get_Am_FUV(pG,i,j,k,Am0);
  *eta_A = Bsq / (pG->U[k][j][i].d*Omega_0*Am);  
  *eta_O=0.0;
  *eta_H=0.0;
//  printf("%i %i %i %.9G %.9G %.9G %.9G %.9G %.9G %.9G %.9G %s \n", i, j, k, pG->U[k][j][i].d, Bsq, *eta_A, Bsq / (pG->U[k][j][i].d*Omega_0*Am), Am, Am0, eta_A0, Q_A , "eol");
  return;
}


static Real get_Am_FUV(const GridS *pG, const int i, const int j, const int k, const Real Am0)
{
  Real Am, Am_FUV, Am_base, x, y, z;
  Real ndelz = 8.;
  static int kbase1,kbase2;
  int ki, kg;
  int idisp,jdisp,i_index,j_index;

  for (ki=pG->ks;ki<=pG->ke;ki++){
    cc_pos(pG,i,j,ki,&x,&y,&z);
		if (i==0 && j==0 && k==0) {printf("%s %i %i %i %0.9G %0.9G %0.9G \n", "get_Am_FUV a ", myID_Comm_world, k, ki, z, zib, zit);}
    if (z > -zib-delz && z <= -zib-delz+pG->dx3) kbase1 = ki;
    if (z < zit+delz && z >= zit+delz-pG->dx3) kbase2 = ki;
  }

  cc_pos(pG,i,j,k,&x,&y,&z);
	if (i==0 && j==0 && k==0) {printf("%s %i %i %0.9G %0.9G %0.9G %i %i\n", "get_Am_FUV b ", myID_Comm_world, k, z, zib, zit, kbase1, kbase2);}

	//myrho = 
	//Am_FUV = 2.36e16*(ionfrav_FUV/1.e-5)*myrho*pow(R_AU,1.5);
	Am_FUV = 3.3e7*(ionfrac_FUV/1.e-5)*pow(R_AU,-5./4.);

	kg = k+pG->Disp[2]-nghost;
	
	if ((z>=zbtm) && (z<=ztop)){
		if(z<=-(zib+delz)){
			Am = Am_FUV*rho_avg[kg];
		} else if (z>-zib-delz && z<-zib+ndelz*delz){
			Am_base = Am_FUV*rho_avg[kbase1];
			Am = (1.-erf((z+zib*0.9)/delz))*Am_base*0.5+Am0;
		} else if (z>zit-ndelz*delz & z<=zit-ndelz*delz){
			Am = Am0;
		} else if (z>zit-ndelz*delz && z<zit+delz){
			Am_base = Am_FUV*rho_avg[kbase2];
			Am = (erf((z-zit*0.9)/delz)+1.)*Am_base*0.5+Am0;
		} else if (z>=zit+delz){
			Am = Am_FUV*rho_avg[kg];
		}
	} else if (z<zbtm){
		Am = Am_FUV*rho_avg[0];
	} else if (z > ztop){
		Am = Am_FUV*rho_avg[pG->ke+pG->Disp[2]-nghost];
	}

/* Here is our equation for Am in the FUV active layer */
  //Am_FUV = 3.3e7*(ionfrac_FUV/1.e-5)*pow(R_AU,-5./4.);

  // 4 new lines
/*  idisp = pG->Disp[0];
  jdisp = pG->Disp[1];
  i_index = i+gridnumx1*2*nghost+idisp;
  j_index = j+gridnumy1*2*nghost+jdisp;

  if (z <= -(zib+delz)) {
    Am = Am_FUV*pG->U[k][j][i].d;
  } else if (z > -zib-delz && z < -zib+ndelz*delz) {
    Am_base = Am_FUV*gdens[kbase1][j_index][i_index];
    Am = (1.-erf((z+zib*0.9)/delz))*Am_base*0.5+Am0;
  } else if (z >= -zib+ndelz*delz && z <= zit-ndelz*delz) {
    Am = Am0;
  } else if (z > zit-ndelz*delz && z < zit+delz) {
    Am_base = Am_FUV*gdens[kbase2][j_index][i_index];
    Am = (erf((z-zit*0.9)/delz)+1.)*Am_base*0.5+Am0;
  } else if (z >= zit+delz) {
    Am = Am_FUV*pG->U[k][j][i].d;
  }
*/
  //if (i==0 && j==0) {printf("%s %i %i %0.9G %0.9G %0.9G %0.9G %0.9G %0.9G %i %i %0.9G \n", "get_Am_FUV", myID_Comm_world, k, z, Am, zib, zit, delz, ndelz, kbase1, kbase2, ionfrac_FUV);}
	//Am = 1000.0;
  return Am;
}

#endif

/* Enforce mass conservation */
static Real mass_cons(MeshS *pM)
{
  GridS *pG;
  DomainS *pD;
  int i,j,k,ierr;
  Real dVol,mtot=0.0,my_mtot,mratio;

  /* Loop over the root domain */
  if (pM->Domain[0][0].Grid != NULL){
    pG = pM->Domain[0][0].Grid;
    pD = (DomainS*)&(pM->Domain[0][0]);

    dVol = 1.0;
    if (pG->dx1 > 0.0) dVol *= pG->dx1;
    if (pG->dx2 > 0.0) dVol *= pG->dx2;
    if (pG->dx3 > 0.0) dVol *= pG->dx3;

    for (k=pG->ks; k<=pG->ke; k++) {
    for (j=pG->js; j<=pG->je; j++) {
    for (i=pG->is; i<=pG->ie; i++) {
      mtot += pG->U[k][j][i].d;
    }}}
    mtot *= dVol;

#ifdef MPI_PARALLEL
    ierr = MPI_Reduce(&(mtot), &(my_mtot), 1,
           MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
#else
    my_mtot = mtot;
#endif

    if(myID_Comm_world==0){  /* I'm the parent */

      dVol = (pD->MaxX[0]-pD->MinX[0])
            *(pD->MaxX[1]-pD->MinX[1])*(pD->MaxX[2]-pD->MinX[2]);
      my_mtot /= dVol;

      /* calculate the mass conversion ratio */
      mratio = mass / my_mtot;
			printf("%s %0.9G \n", "Mass Cons", mratio);
		
    }
  }

#ifdef MPI_PARALLEL
  ierr = MPI_Bcast(&(mratio), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  return mratio;
}

#ifdef RESISTIVITY
static void density_profile(MeshS *pM)
{ 
  GridS *pG;
  DomainS *pD;
  int i,j,k,kg,kdisp;
  double darea,Lx,Ly;

#ifdef MPI_PARALLEL
  double *my_rho;
  double *g_rho;
  int ierr,myID_Comm_Domain;
#endif
// At level=0, there is only one domain 
  pG = pM->Domain[0][0].Grid;

  if (pG != NULL) {
    pD = (DomainS*)&(pM->Domain[0][0]);

    for (k=0; k<pM->Nx[2]; k++) {
      rho_avg[k] = 0.0;
    } 
    kdisp=pG->Disp[2];
// Compute 1d averaged density 
    for (k=pG->ks; k<=pG->ke; k++) {
      kg=k+kdisp-nghost;
      for (j=pG->js; j<=pG->je; j++) {
        for (i=pG->is; i<=pG->ie; i++) {
          rho_avg[kg] += pG->U[k][j][i].d;
        }}
    }
    Lx = pM->RootMaxX[0] - pM->RootMinX[0];
    Ly = pM->RootMaxX[1] - pM->RootMinX[1];
    darea = (double)(pG->dx1*pG->dx2/(Lx*Ly));
    for (k=0; k<pM->Nx[2]; k++) {
      rho_avg[k] *= darea;
	  //printf("%i  %i  %i  %i  %i  %0.9G  marker in density_profile \n", myID_Comm_world, pG->Disp[0], pG->Disp[1], pG->Disp[2], k, rho_avg[k]);
    }
// The parent sums the density array. 
#ifdef MPI_PARALLEL
    ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
    if(ierr != MPI_SUCCESS)
      ath_error("[density_profile]: MPI_Comm_rank error = %d\n",ierr);
 
    my_rho = (double*) calloc_1d_array(pM->Nx[2],sizeof(double));
    g_rho  = (double*) calloc_1d_array(pM->Nx[2],sizeof(double));

    for (k=0; k<pM->Nx[2]; k++) {
      my_rho[k] = rho_avg[k];

    }

    ierr = MPI_Reduce(my_rho, g_rho, pM->Nx[2],
                    MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);

    if(ierr)
      ath_error("[density_profile]: MPI_Reduce returned error = %d\n",ierr);

    if (myID_Comm_Domain == 0){ /* I'm the parent */
      for (k=0; k<pM->Nx[2]; k++) {
        rho_avg[k] = g_rho[k];
		//printf("%0.9G \n", rho_avg[k]);
      }
    }
    //printf("%i  marker 2 in density_profile \n", myID_Comm_world); 
    MPI_Bcast(rho_avg, pM->Nx[2], MPI_DOUBLE, 0, pD->Comm_Domain);
    //printf("%i  marker 3 in density_profile \n", myID_Comm_world);
#endif
  } /* if (pG != NULL) */

#ifdef MPI_PARALLEL
#ifdef STATIC_MESH_REFINEMENT
  MPI_Bcast(rho_avg, pM->Nx[2], MPI_DOUBLE, 0, MPI_Comm_World);
#endif

  free_1d_array(my_rho);
  free_1d_array(g_rho);

#endif
  return;
}

static void ionization_rate(MeshS *pM)
{
  int k,ks,ke;
  GridS *pG;
  DomainS *pD;

  ks = nghost;
  ke = ks + pM->Nx[2]-1;

/*------------------- calculate column densities --------------
 */
  for (k=0; k<=ke+nghost; k++) {
    Sigma_top[k] = 0.0;
    Sigma_bot[k] = 0.0;
  }
  Sigma_top[ks] = 0.5*rho_avg[0];
  Sigma_bot[ke] = 0.5*rho_avg[pM->Nx[2]-1];

  for (k=ks+1;k<=ke;k++) {
    Sigma_top[k] = Sigma_top[k-1] + 0.5*(rho_avg[k-ks-1]+rho_avg[k-ks]);
  }
  for (k=ke+1;k<=ke+nghost;k++) {
    Sigma_top[k] = Sigma_top[ke] + Sigma_bot[ke];
  }
  for (k=0;k<=ke-1;k++) {
    Sigma_bot[k] = Sigma_top[ke+1]-Sigma_top[k];
  }

  /* convert to physical unit */
  for (k=0;k<=ke+nghost;k++) {
    Sigma_top[k] *= mySigma*pM->dx[2];
    Sigma_bot[k] *= mySigma*pM->dx[2];
  }


/*-------------- determine location for FUV ionization front --*/

  int check_bot = 1;
  int check_top = 1;

/* At level=0, there is only one domain */
  pG = pM->Domain[0][0].Grid;

  if (pG != NULL) {
    pD = (DomainS*)&(pM->Domain[0][0]);
  }
  for (k=ks;k<=ke;k++){
    if (Sigma_bot[k] <= Sig_FUV && check_bot == 1) {
      zib = zbtm + ((Real)(k-nghost)+0.5)*pM->dx[2];
      check_bot = 0;
    }
    if (Sigma_top[k] >= Sig_FUV && check_top == 1) {
      zit = zbtm + ((Real)(k-nghost)+0.5)*pM->dx[2];
      check_top = 0;
    }
  }
/* Make zit and zib both positive numbers */
  zib = fabs(zib);
  zit = fabs(zit);


/*-------------- calculate the ionization rate --------------
 */
// DELETED A BUNCH OF STUFF
  return;
}
#endif /* RESISTIVITY */



void Userwork_before_loop(MeshS *pM)
{
  //printf("%i  marker 1 in Userwork_before_loop \n", myID_Comm_world);
  GridS *pG;
  DomainS *pD;
/* At level=0, there is only one domain */
  pG = pM->Domain[0][0].Grid;
  if (pG != NULL) {
    pD = (DomainS*)&(pM->Domain[0][0]);
  }
  nxgrid = pD->NGrid[0];
  nygrid = pD->NGrid[1];
  nzgrid = pD->NGrid[2];
  gridnumx1 = pG->Disp[0]/pG->Nx[0];
  gridnumy1 = pG->Disp[1]/pG->Nx[1];
  gridnumz1 = pG->Disp[2]/pG->Nx[2];

  gdens1D = (double*)calloc_1d_array((pM->Nx[2]+2*nghost*nzgrid)*(pM->Nx[1]+2*nghost*nygrid)*(pM->Nx[0]+2*nghost*nxgrid),sizeof(double));
  ldens   = (double***)calloc_3d_array(pG->Nx[2]+2*nghost,pG->Nx[1]+2*nghost,pG->Nx[0]+2*nghost,sizeof(double));
  gdens   = (double***)calloc_3d_array(pM->Nx[2]+2*nghost*nzgrid,pM->Nx[1]+2*nghost*nygrid,pM->Nx[0]+2*nghost*nxgrid,sizeof(double));

#ifdef RESISTIVITY
  int i,j,mynz=pM->Nx[2]+2*nghost;
  Real Mstar,Tdisk,R_AU;
  Real Sigma,myT,myCs,myH,myeta,myQ_H,myQ_A;
  char *fname=NULL;
  char line[MAXLEN];
  FILE *fp;
  //printf("%i  marker 2 in Userwork_before_loop \n", myID_Comm_world);

  if (par_getd_def("problem","CASE",1) == 3) {

/* initialize 1d (horizontally averaged) for ionization rate calculation */
    rho_avg  = (Real*)calloc_1d_array(pM->Nx[2],sizeof(Real));
    Sigma_top= (Real*)calloc_1d_array(mynz,sizeof(Real));
    Sigma_bot= (Real*)calloc_1d_array(mynz,sizeof(Real));
    ion_rate = (Real*)calloc_1d_array(mynz,sizeof(Real));
	
    //printf("%s \n", "MARKER 7 inside Userwork_before_loop");
/* calculate the horizontally averaged density profile */
    //printf("%i  marker in Userwork_before_loop before density_profile \n", myID_Comm_world);
    density_profile(pM);
    //printf("%i  marker in Userwork_before_loop before density_profile \n", myID_Comm_world);
/* claculate the ionization rate profile */
    ionization_rate(pM);
    //printf("%s \n", "MARKER 9 inside Userwork_before_loop");
    //combine_densities(pM);
    //printf("%s \n", "MARKER 10 inside Userwork_before_loop");
  }

/* Final thing to do is to calculate eta so that it can be ready for the data_out routine in the main loop */

  int k,nl,nd;
  int il, iu, is, ie;
  int jl, ju, js, je;
  int kl, ku, ks, ke;

/* Calculate the magnetic diffusivity array */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {

        pG=pM->Domain[nl][nd].Grid;
        is = pG->is;
        ie = pG->ie;
        js = pG->js;
        je = pG->je;
        ks = pG->ks;
        ke = pG->ke;

       il = is - nghost;
       iu = ie + nghost;
       if (pG->Nx[1] > 1){
         jl = js - nghost;
         ju = je + nghost;
       } else {
         jl = js;
         ju = je;
       }
       if (pG->Nx[2] > 1){
         kl = ks - nghost;
         ku = ke + nghost;
       } else {
         kl = ks;
         ku = ke;
       }

       for (k=kl; k<=ku; k++) {
       for (j=jl; j<=ju; j++) {
       for (i=il; i<=iu; i++) {

          get_eta_user(pG, i,j,k, &(pG->eta_Ohm[k][j][i]),
                               &(pG->eta_Hall[k][j][i]), &(pG->eta_AD[k][j][i]));

       }}}

      }
    }
  }

#endif

  return;
}

void Userwork_in_loop(MeshS *pM)
{
  GridS *pG;
  DomainS *pD;
  int nl,nd,i,j,k;
  Real mratio, press, temp;

/* Enforce mass conservation */
  mratio = mass_cons(pM);

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){

        pG = pM->Domain[nl][nd].Grid;
//        for (k=pG->ks; k<=pG->ke; k++) {
//          for (j=pG->js; j<=pG->je; j++) {
//            for (i=pG->is; i<=pG->ie; i++) {
//              pG->U[k][j][i].d = MAX(pG->U[k][j][i].d * mratio, D_FLOOR);
        for (k=pG->ks; k<=pG->ke; k++) {
          for (j=pG->js; j<=pG->je; j++) {
            for (i=pG->is; i<=pG->ie; i++) {
              pG->U[k][j][i].d = MAX(pG->U[k][j][i].d * mratio, D_FLOOR);
              //printf("%s %i %i %i %i %i %i \n",     "pG->ks pG->ke pG->js pG->je pG->is pG->ie", pG->ks, pG->ke, pG->js, pG->je, pG->is, pG->ie);
            }}}
  }}}

#ifdef RESISTIVITY
	if (count1 == 0) {
		lastRecalcTime = pG->time;	
	}
	count1 = count1 + 1;
  /* user defined diffusivities from look-up tables */
  if (par_geti_def("problem","CASE",1) == 3) { /* requires CASE = 3 */
		count1 = count1 + 1;
		//if (count1 % 10 == 0) {
		if (lastRecalcTime + dtRecalc < pG->time) {

			//clock_t start1 = clock(), diff1;
			density_profile(pM);
			//diff1 = clock() - start1;
			//int msec1 = diff1 * 1000 / CLOCKS_PER_SEC;
			//printf("density_profile Time taken %d seconds %d milliseconds \n", msec1/1000.0, msec1%1000);

			//clock_t start2 = clock(), diff2;
		  ionization_rate(pM);
			//diff2 = clock() - start2;
			//int msec2 = diff2 * 1000 / CLOCKS_PER_SEC;
			//printf("ionization_rate Time taken %d seconds %d milliseconds \n", msec2/1000, msec2%1000);

			//clock_t start3 = clock(), diff3;
			//combine_densities(pM);
			//diff3 = clock() - start3;
			//int msec3 = diff3 * 1000 / CLOCKS_PER_SEC;
			//printf("combine_densities Time taken %d seconds %d milliseconds \n", msec3/1000, msec3%1000);

    	//density_profile(pM);
    	//ionization_rate(pM);
    	//combine_densities(pM);
		}
  }
#endif


/* Apply a pressure and/or temperature floor */
#ifdef ADIABATIC
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){
        pG = pM->Domain[nl][nd].Grid;
        for (k=pG->ks; k<=pG->ke; k++) {
          for (j=pG->js; j<=pG->je; j++) {
            for (i=pG->is; i<=pG->ie; i++) {
              // calculate pressure = (Etot - KE - ME)*Gamma_1 
	    	  press = pG->U[k][j][i].E 
					   - 0.5*(SQR(pG->U[k][j][i].M1)
                            + SQR(pG->U[k][j][i].M2)
                            + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d
			           - 0.5*(SQR(pG->U[k][j][i].B1c)
                            + SQR(pG->U[k][j][i].B2c)
                            + SQR(pG->U[k][j][i].B3c));
			  press *= Gamma_1;
			  // enforce pressure floor
			  press = MAX(press,T_FLOOR*pG->U[k][j][i].d);
			  press = MAX(press,P_FLOOR);
			  // re calculate total energy
			  pG->U[k][j][i].E = 0.5*(SQR(pG->U[k][j][i].M1)
                                    + SQR(pG->U[k][j][i].M2)
                                    + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d
							   + 0.5*(SQR(pG->U[k][j][i].B1c)
                                    + SQR(pG->U[k][j][i].B2c)
                                    + SQR(pG->U[k][j][i].B3c))
						       + press/Gamma_1;
            }}}
  }}}

#endif /* ADIABATIC */


/* check for nans and quit if they exist */
  int isNans=0;
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){
        pG = pM->Domain[nl][nd].Grid;
		for (k=pG->ks-nghost; k<=pG->ke+nghost; k++) {
		  for (j=pG->js-nghost; j<=pG->je+nghost; j++) {
			for (i=pG->is-nghost; i<=pG->ie+nghost; i++) {			
			  if (pG->U[k][j][i].d != pG->U[k][j][i].d){printf("%s %i %i %i %i \n",     "NAN FOUND IN d   at ", i, j, k, nd); isNans=1;}
			  if (pG->U[k][j][i].d != pG->U[k][j][i].d){printf("%s %i %i %i %i \n",     "NAN FOUND IN d   at ", i, j, k, nd); isNans=1;}
			  if (pG->U[k][j][i].M1 != pG->U[k][j][i].M1){printf("%s %i %i %i %i \n",   "NAN FOUND IN M1  at ", i, j, k, nd); isNans=1;}
			  if (pG->U[k][j][i].M2 != pG->U[k][j][i].M2){printf("%s %i %i %i %i \n",   "NAN FOUND IN M2  at ", i, j, k, nd); isNans=1;}
              if (pG->U[k][j][i].M3 != pG->U[k][j][i].M3){printf("%s %i %i %i %i \n",   "NAN FOUND IN M3  at ", i, j, k, nd); isNans=1;}
			  #ifdef adiabatic
			  if (pG->U[k][j][i].E != pG->U[k][j][i].E){printf("%s %i %i %i %i \n",     "NAN FOUND IN E   at ", i, j, k, nd); isNans=1;}
			  #endif //adiabatic
			  if (pG->U[k][j][i].B1c != pG->U[k][j][i].B1c){printf("%s %i %i %i %i \n", "NAN FOUND IN B1c at ", i, j, k, nd); isNans=1;}
			  if (pG->U[k][j][i].B2c != pG->U[k][j][i].B2c){printf("%s %i %i %i %i \n", "NAN FOUND IN B2c at ", i, j, k, nd); isNans=1;}
              if (pG->U[k][j][i].B3c != pG->U[k][j][i].B3c){printf("%s %i %i %i %i \n", "NAN FOUND IN B3c at ", i, j, k, nd); isNans=1;}
            }}}
  }}}
  if (isNans==1) {exit(0);}

  return;
}




void Userwork_in_loop2(MeshS *pM)
{
  GridS *pG;
  DomainS *pD;
  int nl,nd,i,j,k;
  Real mratio, press, temp, P_FLOOR2;

/* Apply a pressure and/or temperature floor */
#ifdef ADIABATIC
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){
        pG = pM->Domain[nl][nd].Grid;
        //for (k=pG->ks; k<=pG->ke; k++) {
		for (k=pG->ks-nghost; k<=pG->ke+nghost; k++) {
          for (j=pG->js-nghost; j<=pG->je+nghost; j++) {
            for (i=pG->is-nghost; i<=pG->ie+nghost; i++) {
              // calculate pressure = (Etot - KE - ME)*Gamma_1 
	    	  press = pG->U[k][j][i].E 
					   - 0.5*(SQR(pG->U[k][j][i].M1)
                            + SQR(pG->U[k][j][i].M2)
                            + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d
			           - 0.5*(SQR(pG->U[k][j][i].B1c)
                            + SQR(pG->U[k][j][i].B2c)
                            + SQR(pG->U[k][j][i].B3c));
			  press *= Gamma_1;
			  // enforce pressure floor
			  press = MAX(press,T_FLOOR*pG->U[k][j][i].d);
			  press = MAX(press,P_FLOOR);
			  // re calculate total energy
			  pG->U[k][j][i].E = 0.5*(SQR(pG->U[k][j][i].M1)
                                    + SQR(pG->U[k][j][i].M2)
                                    + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d
							   + 0.5*(SQR(pG->U[k][j][i].B1c)
                                    + SQR(pG->U[k][j][i].B2c)
                                    + SQR(pG->U[k][j][i].B3c))
						       + press/Gamma_1;
            }}}
  }}}
#endif /* ADIABATIC */

  return;
}


void combine_densities(MeshS *pM)
{

#ifdef MPI_PARALLEL
  int i,j,k,m,n,p;
  int is,ie,js,je,ks,ke;
  long int num_ldens;
  long int count,koffset,joffset,ioffset;
  GridS *pG;
  DomainS *pD;
  int ierr,myID_Comm_Domain;

// At level=0, there is only one domain 
  pG = pM->Domain[0][0].Grid; 
  if (pG != NULL) {
    pD = (DomainS*)&(pM->Domain[0][0]);
  }

  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;

  num_ldens = (pG->Nx[0]+2*nghost)*(pG->Nx[1]+2*nghost)*(pG->Nx[2]+2*nghost);

  ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
  if(ierr != MPI_SUCCESS)
    ath_error("[combine_densities]: MPI_Comm_rank error = %d\n",ierr);

// Set up density array on each grid 
  for (k=ks-nghost; k<=ke+nghost; k++) {
  for (j=js-nghost; j<=je+nghost; j++) {
  for (i=is-nghost; i<=ie+nghost; i++) {
      ldens[k][j][i] = pG->U[k][j][i].d;
  }}}

  ierr = MPI_Allgather(**ldens,num_ldens,MPI_DOUBLE,gdens1D,num_ldens,MPI_DOUBLE,pD->Comm_Domain);
  if(ierr)
    ath_error("[combine_densities]: MPI_Allgather call returned error = %d\n",ierr);

  count = 0;
  koffset = 0;
  for (p=0; p<=nzgrid-1; p++) {
    joffset = 0;
    for (n=0; n<=nygrid-1; n++) {
      ioffset = 0;
      for (m=0; m<=nxgrid-1; m++) {

        for (k=ks-nghost; k<=ke+nghost; k++) {
        for (j=js-nghost; j<=je+nghost; j++) {
        for (i=is-nghost; i<=ie+nghost; i++) {
          count++;
          gdens[k+koffset][j+joffset][i+ioffset] = gdens1D[count];
        }}}

        ioffset = ioffset+2*nghost+pG->Nx[0];
      }
    joffset = joffset+2*nghost+pG->Nx[1];
    }
  koffset = koffset+2*nghost+pG->Nx[2];
  }
#endif
  return;
}















void Userwork_after_loop(MeshS *pM)
{
#ifdef RESISTIVITY
  if (rho_avg  != NULL)  free_1d_array(rho_avg);
  if (Sigma_top != NULL) free_1d_array(Sigma_top);
  if (Sigma_bot != NULL) free_1d_array(Sigma_bot);
  if (ion_rate != NULL)  free_1d_array(ion_rate);
  if (rho_tab  != NULL)  free_1d_array(rho_tab);
  if (zeta_tab != NULL)  free_1d_array(zeta_tab);
  if (etaO_tab != NULL)  free_2d_array(etaO_tab);
  if (Q_H_tab != NULL)  free_2d_array(Q_H_tab);
  if (Q_A_tab != NULL)  free_2d_array(Q_A_tab);
#endif
  if (ldens != NULL)    free_3d_array(ldens);
  if (gdens != NULL)    free_3d_array(gdens);
}

/*------------------------------------------------------------------------------
 */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/*! \fn double ran2(long int *idum)
 *  \brief Extracted from the Numerical Recipes in C (version 2) code.  Modified
 *   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003 
 *
 * Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1.
 */

double ran2(long int *idum)
{
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX

/*----------------------------------------------------------------------------*/
/*! \fn static Real UnstratifiedDisk(const Real x1, const Real x2,const Real x3)
 *  \brief tidal potential in 3D shearing box */
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
#ifndef FARGO
  phi -= qshear*Omega_0*Omega_0*x1*x1;
#endif
  return phi;
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real VertGrav(const Real x1, const Real x2, const Real x3)
 *  \brief potential for vertical component of gravity */
static Real VertGrav(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0,z;

/* If outflow boundaries are used in z, we just use the normal
   z potential.  Otherwise, we ensure periodicity and also
   smooth the potential near the vertical boundaries */

  if (zbc_out == 1) {
    z = x3;
    phi += 0.5*Omega_0*Omega_0*z*z;
  } else {
    if(x3 > ztop)
      z=x3-ztop+zbtm;
    else if (x3 < zbtm)
      z=x3-zbtm+ztop;
    else
      z=x3;
    phi += 0.5*Omega_0*Omega_0*
     (SQR(fabs(ztop)-sqrt(SQR(fabs(ztop)-fabs(z)) + 0.01)));
  }
  return phi;
}











































/*! \fn static void strat_ix3(GridS *pG)
 *  \brief  Here is the lower z outflow boundary. 
            The basic idea is that the pressure and density
            are exponentially extrapolated in the ghost zones 
            assuming a constant temperature there (i.e., an
            isothermal atmosphere). The z velocity (NOT the
            momentum) are set to zero in the ghost zones in the
            case of the last lower physical zone having an inward
            flow.  All other variables are extrapolated into the 
            ghost zones with zero slope.
*/
static void strat_ix3(GridS *pG)
{
  int ks = pG->ks;
  int ie = pG->ie;
  int je = pG->je;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */
  Real x1,x2,x3;
  Real press,pressks,Tks;
  static Real x3b;

  Real maxB=0.0;
  int maxBi = 0;
  int maxBj = 0;
  int maxBk = 0;

  x3b = zbtm+0.5*pG->dx3;

  if (pG->Nx[0] > 1){
    iu = pG->ie + nghost;
    il = pG->is - nghost;
  } else {
    iu = pG->ie;
    il = pG->is;
  }
  if (pG->Nx[1] > 1){
    ju = pG->je + nghost;
    jl = pG->js - nghost;
  } else {
    ju = pG->je;
    jl = pG->js;
  }

#ifdef MHD
/* Copy field components from last physical zone */
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,ks-k,&x1,&x2,&x3);
        pG->B1i[ks-k][j][i] = pG->B1i[ks][j][i]*exp(-(x3*x3-x3b*x3b));
        pG->U[ks-k][j][i].B1c = pG->U[ks][j][i].B1c*exp(-(x3*x3-x3b*x3b));
        pG->B2i[ks-k][j][i] = pG->B2i[ks][j][i]*exp(-(x3*x3-x3b*x3b));
        pG->U[ks-k][j][i].B2c = pG->U[ks][j][i].B2c*exp(-(x3*x3-x3b*x3b));
        pG->B3i[ks-k][j][i] = pG->B3i[ks][j][i];
      }
    }
  }
/* Cell centered Bz is calculated from average of interface fields */
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pG->U[ks-k][j][i].B3c = 0.5*(pG->B3i[ks-k+1][j][i]+pG->B3i[ks-k][j][i]);
      }
    }
  }
#endif /* MHD */

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,ks-k,&x1,&x2,&x3);
/* First calculate the effective gas temperature in the last physical zone */
#ifdef ADIABATIC
		pressks = pG->U[ks][j][i].E - 0.5*(SQR(pG->U[ks][j][i].M1)
                                         + SQR(pG->U[ks][j][i].M2)
                                         + SQR(pG->U[ks][j][i].M3))/pG->U[ks][j][i].d;
#ifdef MHD
        pressks -= 0.5*(SQR(pG->U[ks][j][i].B1c)
                      + SQR(pG->U[ks][j][i].B2c)
                      + SQR(pG->U[ks][j][i].B3c));

#endif /* MHD */
        pressks *= Gamma_1;
	    //pressks = MAX(pressks,T_FLOOR*pG->U[ks][j][i].d);
		//pressks = MAX(pressks,P_FLOOR);
		Tks = pressks/pG->U[ks][j][i].d;
#else
        Tks = 0.5*Omega_0*Omega_0;
#endif /* ADIABATIC */

/* Now extrapolate the density to balance gravity assuming a constant temperature in the ghost zones */
        pG->U[ks-k][j][i].d = pG->U[ks][j][i].d;
        //pG->U[ks-k][j][i].d = pG->U[ks][j][i].d*exp(-(x3*x3-x3b*x3b)/(2.0*Tks/(Omega_0*Omega_0)));
		//pG->U[ks-k][j][i].d = pG->U[ks][j][i].d*exp(-(x3-x3b)/(2.0*Tks/(Omega_0*Omega_0)));
		//pG->U[ks-k][j][i].d = (1-(ks-k)/(2*nghost))*pG->U[ks][j][i].d;
/* Copy the velocities, but not the momenta --- important because of the density extrapolation above */
        pG->U[ks-k][j][i].M1 = pG->U[ks][j][i].M1/pG->U[ks][j][i].d*pG->U[ks-k][j][i].d;
        pG->U[ks-k][j][i].M2 = pG->U[ks][j][i].M2/pG->U[ks][j][i].d*pG->U[ks-k][j][i].d;
/* If there's inflow into the grid, set the normal velocity to zero */
        if (pG->U[ks][j][i].M3 >= 0.0) {
          pG->U[ks-k][j][i].M3 = 0.0;
        } else {
		  pG->U[ks-k][j][i].M3 = pG->U[ks][j][i].M3/pG->U[ks][j][i].d*pG->U[ks-k][j][i].d;
        }
#ifdef ADIABATIC
		press = MAX(P_FLOOR,pG->U[ks-k][j][i].d*Tks);
		pG->U[ks-k][j][i].E = press/Gamma_1
        + 0.5*(SQR(pG->U[ks-k][j][i].M1) 
             + SQR(pG->U[ks-k][j][i].M2)
             + SQR(pG->U[ks-k][j][i].M3))/pG->U[ks-k][j][i].d;
#ifdef MHD
        pG->U[ks-k][j][i].E += 0.5*(SQR(pG->U[ks-k][j][i].B1c) 
								  + SQR(pG->U[ks-k][j][i].B2c) 
								  + SQR(pG->U[ks-k][j][i].B3c));
#endif /* MHD */
#endif /* ADIABATIC */












// CHECKS AND SUCH, DOESN'T ACTUALLY DO ANYTING CODE-WISE
#ifdef adiabatic
		Real tempPress=Gamma_1*(pG->U[ks-k][j][i].E 
    							- 0.5*(SQR(pG->U[ks-k][j][i].M1) 
									 + SQR(pG->U[ks-k][j][i].M2) 
                                     + SQR(pG->U[ks-k][j][i].M3))/pG->U[ks-k][j][i].d 
								- 0.5*(SQR(pG->U[ks-k][j][i].B1c) 
						             + SQR(pG->U[ks-k][j][i].B2c) 
                                     + SQR(pG->U[ks-k][j][i].B3c)));
		if (tempPress<0) {printf("%s \n", "NEGATIVE PRESSURE AT LOC2 in Gz routine");}

		// find the largest B in the gzs
		if (0.5*(SQR(pG->U[ks-k][j][i].B1c)+SQR(pG->U[ks-k][j][i].B2c)+SQR(pG->U[ks-k][j][i].B3c))/pG->U[ks-k][j][i].E > maxB) {
		maxB = 0.5*(SQR(pG->U[ks-k][j][i].B1c)+SQR(pG->U[ks-k][j][i].B2c)+SQR(pG->U[ks-k][j][i].B3c))/pG->U[ks-k][j][i].E;
        maxBi = i; maxBj = j; maxBk = ks-k;}
#endif //adiabatic		
      }
    }
  }
  return;
}

















/*! \fn static void strat_ox3(GridS *pG)
 *  \brief  Here is the upper z outflow boundary. 
            The basic idea is that the pressure and density
            are exponentially extrapolated in the ghost zones 
            assuming a constant temperature there (i.e., an
            isothermal atmosphere). The z velocity (NOT the
            momentum) are set to zero in the ghost zones in the
            case of the last upper physical zone having an inward
            flow.  All other variables are extrapolated into the 
            ghost zones with zero slope.
*/

static void strat_ox3(GridS *pG)
{
  int ke = pG->ke;
  int ie = pG->ie;
  int je = pG->je;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */
  Real x1,x2,x3;
  Real press,presske,Tke;
  static Real x3t;

  Real maxB=0.0;
  int maxBi = 0;
  int maxBj = 0;
  int maxBk = 0;

  x3t = ztop-0.5*pG->dx3;

  if (pG->Nx[0] > 1){
    iu = pG->ie + nghost;
    il = pG->is - nghost;
  } else {
    iu = pG->ie;
    il = pG->is;
  }
  if (pG->Nx[1] > 1){
    ju = pG->je + nghost;
    jl = pG->js - nghost;
  } else {
    ju = pG->je;
    jl = pG->js;
  }

#ifdef MHD
/* Copy field components from last physical zone */
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,ke+k,&x1,&x2,&x3);
        pG->B1i[ke+k][j][i] = pG->B1i[ke][j][i]*exp(-(x3*x3-x3t*x3t));
        pG->U[ke+k][j][i].B1c = pG->U[ke][j][i].B1c*exp(-(x3*x3-x3t*x3t));
        pG->B2i[ke+k][j][i] = pG->B2i[ke][j][i]*exp(-(x3*x3-x3t*x3t));
        pG->U[ke+k][j][i].B2c = pG->U[ke][j][i].B2c*exp(-(x3*x3-x3t*x3t));
       
/* ke+1 interface field is already handled via CT */
        if (k > 1) {pG->B3i[ke+k][j][i] = pG->B3i[ke+1][j][i];}
      }
    }
  }
/* Update cell centered Bz by averaging interface fields */
  for (k=1;k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        if (k != nghost) {
          pG->U[ke+k][j][i].B3c = 0.5*(pG->B3i[ke+k+1][j][i]+pG->B3i[ke+k][j][i]);
        } else {
          pG->U[ke+k][j][i].B3c = pG->B3i[ke+k][j][i];
        }
      }
    }
  }
#endif /* MHD */

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,ke+k,&x1,&x2,&x3);
#ifdef ADIABATIC
        presske = pG->U[ke][j][i].E - 0.5*(SQR(pG->U[ke][j][i].M1)
                                         + SQR(pG->U[ke][j][i].M2)
                                         + SQR(pG->U[ke][j][i].M3))/pG->U[ke][j][i].d;
#ifdef MHD
        presske -= 0.5*(SQR(pG->U[ke][j][i].B1c)
                      + SQR(pG->U[ke][j][i].B2c)
                      + SQR(pG->U[ke][j][i].B3c));
#endif /* MHD */
        presske *= Gamma_1;
	    //presske = MAX(presske,T_FLOOR*pG->U[ke][j][i].d);
		//presske = MAX(presske,P_FLOOR);
        Tke = presske/pG->U[ke][j][i].d;
#else
        Tke = 0.5*Omega_0*Omega_0;
#endif /* ADIABATIC */
/* Now extrapolate the density to balance gravity assuming a constant temperature in the ghost zones */
		pG->U[ke+k][j][i].d = pG->U[ke][j][i].d;
        //pG->U[ke+k][j][i].d = pG->U[ke][j][i].d*exp(-(x3*x3-x3t*x3t)/(2.0*Tke/(Omega_0*Omega_0)));
		//pG->U[ke+k][j][i].d = pG->U[ke][j][i].d*exp(-(x3-x3t)/(2.0*Tke/(Omega_0*Omega_0)));
		//pG->U[ke+k][j][i].d = (1-(k-ke)/(2*nghost))*pG->U[ke][j][i].d;
/* Copy the velocities, but not the momenta --- important because of the density extrapolation above */
        pG->U[ke+k][j][i].M1 = pG->U[ke][j][i].M1/pG->U[ke][j][i].d*pG->U[ke+k][j][i].d;
        pG->U[ke+k][j][i].M2 = pG->U[ke][j][i].M2/pG->U[ke][j][i].d*pG->U[ke+k][j][i].d;
/* If there's inflow into the grid, set the normal velocity to zero */
        if (pG->U[ke][j][i].M3 <= 0.0) {
          pG->U[ke+k][j][i].M3 = 0.0;
        } else {
          pG->U[ke+k][j][i].M3 = pG->U[ke][j][i].M3/pG->U[ke][j][i].d*pG->U[ke+k][j][i].d;
        }
#ifdef ADIABATIC
        press = MAX(P_FLOOR, pG->U[ke+k][j][i].d*Tke);
        pG->U[ke+k][j][i].E = press/Gamma_1
        + 0.5*(SQR(pG->U[ke+k][j][i].M1) 
             + SQR(pG->U[ke+k][j][i].M2)
             + SQR(pG->U[ke+k][j][i].M3))/pG->U[ke+k][j][i].d;
#ifdef MHD
        pG->U[ke+k][j][i].E += 0.5*(SQR(pG->U[ke+k][j][i].B1c)
		                         + SQR(pG->U[ke+k][j][i].B2c) 
						         + SQR(pG->U[ke+k][j][i].B3c));
#endif /* MHD */
#endif /* ADIABATIC */












// CHECKS AND SUCH, DOESN'T ACTUALLY DO ANYTING CODE-WISE
#ifdef adiabatic
		Real tempPress=Gamma_1*(pG->U[ke+k][j][i].E 
								- 0.5*(SQR(pG->U[ke+k][j][i].M1) 
									 + SQR(pG->U[ke+k][j][i].M2) 
                                     + SQR(pG->U[ke+k][j][i].M3))/pG->U[ke+k][j][i].d 
								- 0.5*(SQR(pG->U[ke+k][j][i].B1c) 
						             + SQR(pG->U[ke+k][j][i].B2c) 
                                     + SQR(pG->U[ke+k][j][i].B3c)));
		if (tempPress<0) {printf("%s \n", "NEGATIVE PRESSURE AT LOC4 in Gz routine");}

		// find the largest B in the gzs
		if (0.5*(SQR(pG->U[ke+k][j][i].B1c)+SQR(pG->U[ke+k][j][i].B2c)+SQR(pG->U[ke+k][j][i].B3c))/pG->U[ke+k][j][i].E > maxB) {
		maxB = 0.5*(SQR(pG->U[ke+k][j][i].B1c)+SQR(pG->U[ke+k][j][i].B2c)+SQR(pG->U[ke+k][j][i].B3c))/pG->U[ke+k][j][i].E;
        maxBi = i; maxBj = j; maxBk = ke+k;}
#endif //adiabatic
      }
    }
  }

  //printf("%s %.2i %.2i %.2i %.9G %.9G %.9G \n", "LARGEST B IN GZs: ", maxBi, maxBj, maxBk, 
  //0.5*(SQR(pG->U[maxBk][maxBj][maxBi].B1c)+SQR(pG->U[maxBk][maxBj][maxBi].B2c)+SQR(pG->U[maxBk][maxBj][maxBi].B3c))/pG->U[maxBk][maxBj][maxBi].E,
  //0.5*(SQR(pG->U[maxBk][maxBj][maxBi].M1)+SQR(pG->U[maxBk][maxBj][maxBi].M2)+SQR(pG->U[maxBk][maxBj][maxBi].M3))/pG->U[maxBk][maxBj][maxBi].E,
  //0.5*((SQR(pG->U[maxBk][maxBj][maxBi].M1)+SQR(pG->U[maxBk][maxBj][maxBi].M2)+SQR(pG->U[maxBk][maxBj][maxBi].M3))+(SQR(pG->U[maxBk][maxBj][maxBi].B1c)+SQR(pG->U[maxBk][maxBj][maxBi].B2c)+SQR(pG->U[maxBk][maxBj][maxBi].B3c)))/pG->U[maxBk][maxBj][maxBi].E);	


  return;

}
































/*----------------------------------------------------------------------------*/
/*! \fn static Real expr_dV2(const GridS *pG, const int i, const int j, 
 *			     const int k)
 *  \brief Computes delta(Vy) 
 */
static Real expr_dV2(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  return (pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  return (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real expr_beta(const GridS *pG, const int i, const int j, 
 *			      const int k)
 *  \brief Computes beta=P/(B^2/8pi)  
 */
static Real expr_beta(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,B2;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef MHD
  B2=pG->U[k][j][i].B1c*pG->U[k][j][i].B1c;
  B2+=pG->U[k][j][i].B2c*pG->U[k][j][i].B2c;
  B2+=pG->U[k][j][i].B3c*pG->U[k][j][i].B3c;

#ifdef ISOTHERMAL
  return (2.0*Iso_csound2*pG->U[k][j][i].d/B2);
#else
  return 0.0;
#endif

#else
  return 0.0;
#endif /* MHD */
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real expr_ME(const GridS *pG, const int i, const int j, 
 *			    const int k)
 *  \brief  Computes B^2/8pi
 */
static Real expr_ME(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,B2;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef MHD
  B2=pG->U[k][j][i].B1c*pG->U[k][j][i].B1c;
  B2+=pG->U[k][j][i].B2c*pG->U[k][j][i].B2c;
  B2+=pG->U[k][j][i].B3c*pG->U[k][j][i].B3c;
  return (B2/2.0);
#else
  return NULL;
#endif
}
/*----------------------------------------------------------------------------*/
/*! \fn static Real expr_KE(const GridS *pG, const int i, const int j, 
 *			    const int k)
 *  \brief Computes dens*(Vx^2+Vy^2+Vz^2)/2
 */
static Real expr_KE(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,Vy,Vx,Vz;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  Vy = (pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  Vy = (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
  Vx = pG->U[k][j][i].M1/pG->U[k][j][i].d;
  Vz = pG->U[k][j][i].M3/pG->U[k][j][i].d;

  return pG->U[k][j][i].d*(Vx*Vx + Vy*Vy + Vz*Vz)/2.0;

}

/*------------------------------------------------------------------------------
 * Hydro history variables:
 * hst_rho_Vx_dVy: Reynolds stress, added as history variable.
 * hst_rho_dVy2: KE in y-velocity fluctuations
 * hst_E_total: total energy (including tidal potential).
 */
/*! \fn static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j, 
 *				  const int k)
 *  \brief Reynolds stress, added as history variable. */
static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  return pG->U[k][j][i].M1*(pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  return pG->U[k][j][i].M1*
    (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
}

/*! \fn static Real hst_rho_dVy2(const GridS *pG, const int i, const int j, 
 *				const int k)
 *  \brief KE in y-velocity fluctuations */
static Real hst_rho_dVy2(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,dVy;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  dVy = (pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  dVy = (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
  return pG->U[k][j][i].d*dVy*dVy;
}

#ifdef ADIABATIC
/*! \fn static Real hst_E_total(const GridS *pG, const int i, const int j, 
 *				const int k)
 *  \brief total energy (including tidal potential). */
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,phi;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  phi = UnstratifiedDisk(x1, x2, x3);

  return pG->U[k][j][i].E + pG->U[k][j][i].d*phi;
}
#endif /* ADIABATIC */

/*------------------------------------------------------------------------------
 * MHD history variables
 * hst_Bx, etc.: Net flux, and Maxwell stress, added as history variables
 */

#ifdef MHD
/*! \fn static Real hst_Bx(const GridS *pG, const int i,const int j,const int k)
 *  \brief x-component of magnetic field */
static Real hst_Bx(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B1c;
}

/*! \fn static Real hst_By(const GridS *pG,const int i,const int j,const int k) 
 *  \brief y-component of magnetic field */
static Real hst_By(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B2c;
}

/*! \fn static Real hst_Bz(const GridS *pG, const int i, const int j, 
 *			  const int k)
 *  \brief z-component of magnetic field */
static Real hst_Bz(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B3c;
}

/*! \fn static Real hst_BxBy(const GridS *pG, const int i, const int j, 
 *			     const int k)
 *  \brief Maxwell stress */
static Real hst_BxBy(const GridS *pG, const int i, const int j, const int k)
{
  return -pG->U[k][j][i].B1c*pG->U[k][j][i].B2c;
}

#endif /* MHD */
#ifdef RESISTIVITY
/*----------------------------------------------------------------------------*/
/*! \fn static Real expr_Lambda(const GridS *pG, const int i, const int j, 
 *                              const int k)
 *  \brief  Computes the Ohmic Elsasser number
 */
static Real expr_Lambda(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,B2;
  B2= SQR(pG->U[k][j][i].B1c)+SQR(pG->U[k][j][i].B2c)+SQR(pG->U[k][j][i].B3c);

  return (B2/(pG->eta_Ohm[k][j][i]*pG->U[k][j][i].d*Omega_0));
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real expr_chi(const GridS *pG, const int i, const int j, 
 *                              const int k)
 *  \brief  Computes the Hall Elsasser number
 */
static Real expr_chi(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,B2;
  B2= SQR(pG->U[k][j][i].B1c)+SQR(pG->U[k][j][i].B2c)+SQR(pG->U[k][j][i].B3c);

  return (B2/(pG->eta_Hall[k][j][i]*pG->U[k][j][i].d*Omega_0));
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real expr_Am(const GridS *pG, const int i, const int j, 
 *                          const int k)
 *  \brief  Computes the ambipolar diffusion Elsasser number
 */
static Real expr_Am(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,B2;
  B2= SQR(pG->U[k][j][i].B1c)+SQR(pG->U[k][j][i].B2c)+SQR(pG->U[k][j][i].B3c);

  return (B2/(pG->eta_AD[k][j][i]*pG->U[k][j][i].d*Omega_0));
}
#endif

/* Here is the output routine to calculate 1D horizontally
   averaged quantities.  Currently, only works without SMR */

/*! \fn static void output_1d(MeshS *pM, OutputS *pOut)
 *  \brief output routine to calculate 1D horizontally
    averaged quantities.  Currently, only outputs at lowest
    refinement level */

static void output_1d(MeshS *pM, OutputS *pOut)
{
  GridS *pGrid;
  DomainS *pD;
  int i,j,k;
  int tot1d,i1d,nzmx,my_nz,kg,kdisp;
  int dnum = pOut->num,nl,nd;
  static int FIRST = 0;
  double darea,**out1d;
  double x1,x2,x3,Lx,Ly,press;
  static double *out_x3;

  FILE *p_1dfile;
  char *fname;
  double area_rat; /* (Grid Volume)/(dx1*dx2*dx3) */

#ifdef MPI_PARALLEL
  double *my_out1d;
  double *g_out1d;
  int zproc;
  int ierr,myID_Comm_Domain;
#endif

#ifdef MHD
#ifdef RESISTIVITY  
  tot1d=18;
#else
  tot1d=15;
#endif /* RESISTIVITY */
#else
  tot1d=7;
#endif /* MHD */
#ifdef ADIABATIC
  tot1d=tot1d+3;
#endif /* ADIABATIC */

  Lx = pM->RootMaxX[0] - pM->RootMinX[0];
  Ly = pM->RootMaxX[1] - pM->RootMinX[1];
  nzmx = pM->Nx[2];

/* At level=0, there is only one domain */

  pGrid = pM->Domain[0][0].Grid;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  pD = (DomainS*)&(pM->Domain[0][0]);

#ifdef MPI_PARALLEL
  int nproc = pD->NGrid[0]*pD->NGrid[1]*pD->NGrid[2];
#endif

#ifdef MPI_PARALLEL
  ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
  if(ierr != MPI_SUCCESS)
    ath_error("[change_rundir]: MPI_Comm_rank error = %d\n",ierr);
#endif
  if (FIRST == 0){
#ifdef MPI_PARALLEL
    if (myID_Comm_Domain == 0) {
#endif
      out_x3 = (double *) calloc_1d_array(nzmx,sizeof(double));
#ifdef MPI_PARALLEL
    }
#endif
  }

  out1d = (double **) calloc_2d_array(nzmx,tot1d,sizeof(double));
#ifdef MPI_PARALLEL
  my_out1d = (double *) calloc_1d_array(nzmx,sizeof(double));
  g_out1d = (double *) calloc_1d_array(nzmx,sizeof(double));
#endif
  for (k=0; k<nzmx; k++) {
    for (i1d=0; i1d<tot1d; i1d++) {
      out1d[k][i1d] = 0.0;
    }
  }
  kdisp=pGrid->Disp[2];

/* First calculate the x3 coordinate and save it to be dumped
   by root in every 1d file */
  if (FIRST == 0) {
#ifdef MPI_PARALLEL
  if (myID_Comm_Domain == 0) {
#endif
    for (k=0; k<nzmx; k++) {
      x3 = pM->RootMinX[2] + (k + 0.5)*pGrid->dx3;
      out_x3[k] = x3;
    }
#ifdef MPI_PARALLEL
  }
#endif
  }

/* Compute 1d averaged variables */
  for (k=ks; k<=ke; k++) {
    kg=k+kdisp-nghost;
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        i1d=0;
        out1d[kg][i1d] += pGrid->U[k][j][i].d;
        i1d++;
#ifdef ISOTHERMAL
        out1d[kg][i1d] += pGrid->U[k][j][i].d*Iso_csound2;
#else
        press           = MAX(Gamma_1*(pGrid->U[k][j][i].E - expr_KE(pGrid,i,j,k)
#ifdef MHD
                                 - expr_ME(pGrid,i,j,k)
#endif
                                ),TINY_NUMBER);
        out1d[kg][i1d] += press;
#endif
#ifdef ADIABATIC
        i1d++;
        out1d[kg][i1d] += press/(Gamma_1*pGrid->U[k][j][i].d);
        i1d++;
        out1d[kg][i1d] += pGrid->U[k][j][i].E;
        i1d++;
        out1d[kg][i1d] += hst_E_total(pGrid,i,j,k);
#endif
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M1)/pGrid->U[k][j][i].d;
        i1d++;
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M2)/pGrid->U[k][j][i].d;
#else
        out1d[kg][i1d] += 0.5*pGrid->U[k][j][i].d*SQR(pGrid->U[k][j][i].M2/pGrid->U[k][j][i].d + qshear*Omega_0*x1);
#endif
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M3)/pGrid->U[k][j][i].d;
        i1d++;
        out1d[kg][i1d] += expr_KE(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_rho_Vx_dVy(pGrid,i,j,k);
#ifdef MHD
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].B1c);
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].B2c);
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].B3c);
        i1d++;
        out1d[kg][i1d] += expr_ME(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_Bx(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_By(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_Bz(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_BxBy(pGrid,i,j,k);
#endif
#ifdef RESISTIVITY
        i1d++;
        out1d[kg][i1d] += expr_Lambda(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += expr_chi(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += expr_Am(pGrid,i,j,k);
#endif
      }
    }
  }

  /* Calculate the (Grid Volume) / (Grid Cell Volume) Ratio */
  area_rat = Lx*Ly/(pGrid->dx1*pGrid->dx2);

/* The parent sums the scal[] array.
 * Note that this assumes (dx1,dx2,dx3) = const. */

#ifdef MPI_PARALLEL 
  for(i1d=0; i1d<tot1d; i1d++){
    for (k=0; k<nzmx; k++) {
      my_out1d[k] = out1d[k][i1d];
    }
    ierr = MPI_Reduce(my_out1d, g_out1d, nzmx,
                      MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
    if(ierr)
      ath_error("[output_1d]: MPI_Reduce call returned error = %d\n",ierr);
    for (k=0; k<nzmx; k++) {
      out1d[k][i1d] = g_out1d[k];
    }
  }
#endif

/* For parallel calculations, only the parent computes the average
 * and writes the output. */
#ifdef MPI_PARALLEL
  if(myID_Comm_Domain == 0){ /* I'm the parent */
#endif

  darea = 1.0/(double)area_rat;
  for (k=0; k<nzmx; k++) {
    for (i1d=0; i1d<tot1d; i1d++) {
      out1d[k][i1d] *= darea;
    }
  }

/* Generate filename */
#ifdef MPI_PARALLEL
  fname = ath_fname("../",pM->outfilename,NULL,NULL,num_digit,dnum,NULL,"1d");
#else
  fname = ath_fname(NULL,pM->outfilename,NULL,NULL,num_digit,dnum,NULL,"1d");
#endif
  if (fname == NULL) {
    ath_error("[output_1d]: Error constructing output filename\n");
    return;
  }

/* open filename */
  p_1dfile = fopen(fname,"w");
  if (p_1dfile == NULL) {
    ath_error("[output_1d]: Unable to open 1d average file %s\n",fname);
    return;
  }

/* Write out data */

  for (k=0; k<nzmx; k++) {
#ifdef ISOTHERMAL
#ifdef MHD
#ifdef RESISTIVITY
      if (k == 0) {
        fprintf(p_1dfile,"# x3     dens    pressure    KEx         KEy         KEz         KE          Reynolds    MEx         MEy         MEz         ME          Bx           By           Bz          Maxwell     Lambda      Chi         Am\n");
      }
      fprintf(p_1dfile,"%G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],
            out1d[k][3],out1d[k][4],out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9],out1d[k][10],out1d[k][11],
            out1d[k][12],out1d[k][13],out1d[k][14],out1d[k][15],out1d[k][16],out1d[k][17]);
#else
      if (k == 0) {
        fprintf(p_1dfile,"# x3     dens    pressure    KEx         KEy         KEz         KE          Reynolds    MEx         MEy         MEz         ME          Bx           By           Bz          Maxwell\n");
      }
      fprintf(p_1dfile,"%G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],
            out1d[k][3],out1d[k][4],out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9],out1d[k][10],out1d[k][11],
            out1d[k][12],out1d[k][13],out1d[k][14]);

#endif /* RESISTIVITY */
#else
    if (k == 0) {
      fprintf(p_1dfile,"# x3     dens  pressure    KEx         KEy         KEz         KE          Reynolds\n");
    }
    fprintf(p_1dfile,"%G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],out1d[k][3],out1d[k][4],
            out1d[k][5],out1d[k][6]);
#endif /* MHD */
#else
#ifdef MHD
    if (k == 0) {
      fprintf(p_1dfile,"# x3     dens    pressure    temperature  E     Etot     KEx         KEy         KEz         KE          Reynolds    MEx         MEy         MEz         ME          Bx           By           Bz          Maxwell\n");
    }
    fprintf(p_1dfile,"%G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],
            out1d[k][3],out1d[k][4],out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9],out1d[k][10],out1d[k][11],
            out1d[k][12],out1d[k][13],out1d[k][14],out1d[k][15],out1d[k][16],out1d[k][17]);
#else
    if (k == 0) {
      fprintf(p_1dfile,"# x3     dens    pressure    temperature  E     Etot     KEx         KEy         KEz         KE          Reynolds\n");
    }
    fprintf(p_1dfile,"%G %G %G %G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],out1d[k][3],out1d[k][4],
            out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9]);
#endif /* MHD */
#endif /* ISOTHERMAL */
  }

  fclose(p_1dfile);
  free(fname);
#ifdef MPI_PARALLEL
  }
#endif

  free_2d_array(out1d); /* Free the memory we malloc'd */
#ifdef MPI_PARALLEL
  free_1d_array(my_out1d); /* Free the memory we malloc'd */
  free_1d_array(g_out1d); /* Free the memory we malloc'd */
#endif
  if (FIRST == 0) {
    FIRST = 1;
  }

return;
}
