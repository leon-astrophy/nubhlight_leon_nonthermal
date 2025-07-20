/******************************************************************************
 *                                                                            *
 * BOUNDS.C                                                                   *
 *                                                                            *
 * PHYSICAL BOUNDARY CONDITIONS                                               *
 *                                                                            *
 ******************************************************************************/

/* include headers */
#include "decs.h"

/* define function used below */
void inflow_check(double *Pr, int ii, int jj, int type); 

/*########################################################################################################*/

/************************************************************************/
/* apply boundary conditions to primitive variables */
void bound_prim(grid_prim_type prim) {

  // start timer //
  timer_start(TIMER_BOUND);

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  // globally inner-x boundary //
  if (global_start[1] == 0) {
#pragma omp parallel for collapse(2)
    JSLOOP(0, N2 - 1) {
      KSLOOP(0, N3 - 1) {
        ISLOOP(-NG, -1) {
          #if N1 < NG
          int   iactive           = NG;
          PLOOP prim[i][j][k][ip] = prim[iactive][j][k][ip];
          pflag[i][j][k]          = pflag[iactive][j][k];
          #else // N1 < NG
          {
            #if X1L_GAS_BOUND == BC_OUTFLOW
            int   iactive           = NG;
            PLOOP prim[i][j][k][ip] = prim[iactive][j][k][ip];
            pflag[i][j][k]          = pflag[iactive][j][k];
            double rescale_fac = ggeom[iactive][j][CENT].g / ggeom[i][j][CENT].g;
            prim[i][j][k][B1] *= rescale_fac;
            prim[i][j][k][B2] *= rescale_fac;
            prim[i][j][k][B3] *= rescale_fac;
            #elif X1L_GAS_BOUND == BC_PROB // X1L_GAS_BOUND == BC_OUTFLOW
            bound_gas_prob_x1l(i, j, k, prim);
            #elif X1L_GAS_BOUND != BC_PERIODIC // X1L_GAS_BOUND == BC_OUTFLOW
            printf("X1L_GAS_BOUND choice %i not supported\n", X1L_GAS_BOUND);
            exit(-1);
            #endif // X1L_GAS_BOUND == BC_OUTFLOW
          }
          #endif // N1 < NG
        }
      }
    }
  } // global_start[1] == 0

  /* fill ghost zone across MPI processes */
  sync_mpi_boundaries_X1L(prim);

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  // globally outer-x boundary //  
  if (global_stop[1] == N1TOT) {
#pragma omp parallel for collapse(2)
    JSLOOP(0, N2 - 1) {
      KSLOOP(0, N3 - 1) {
        ISLOOP(N1, N1 - 1 + NG) {
          #if N1 < NG
          int   iactive           = N1 - 1 + NG;
          PLOOP prim[i][j][k][ip] = prim[iactive][j][k][ip];
          pflag[i][j][k]          = pflag[iactive][j][k];
          #else // N1 < NG
          {
            #if X1R_GAS_BOUND == BC_OUTFLOW
            int   iactive           = N1 - 1 + NG;
            PLOOP prim[i][j][k][ip] = prim[iactive][j][k][ip];
            pflag[i][j][k]          = pflag[iactive][j][k];
            double rescale_fac = ggeom[iactive][j][CENT].g / ggeom[i][j][CENT].g;
            prim[i][j][k][B1] *= rescale_fac;
            prim[i][j][k][B2] *= rescale_fac;
            prim[i][j][k][B3] *= rescale_fac;
            #elif X1R_GAS_BOUND == BC_PROB // X1R_GAS_BOUND == BC_OUTFLOW
            bound_gas_prob_x1r(i, j, k, prim);
            #elif X1R_GAS_BOUND != BC_PERIODIC // X1R_GAS_BOUND == BC_OUTFLOW
            printf("X1R_GAS_BOUND choice %i not supported\n", X1R_GAS_BOUND);
            exit(-1);
            #endif // X1R_GAS_BOUND == BC_OUTFLOW
          }
          #endif // N1 < NG
        }
      }
    }
  } // global_stop[1] == N1TOT

  /* fill ghost zone across MPI processes */
  sync_mpi_boundaries_X1R(prim);

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  // globally inner-y boundary //  
  if (global_start[2] == 0) {
#pragma omp parallel for collapse(2)
    ILOOPALL {
      KSLOOP(0, N3 - 1) {
        JSLOOP(-NG, -1) {
          #if N2 < NG
          int   jactive           = NG;
          PLOOP prim[i][j][k][ip] = prim[i][jactive][k][ip];
          pflag[i][j][k]          = pflag[i][jactive][k];
          #else // N2 < NG
          {
            #if X2L_GAS_BOUND == BC_OUTFLOW
            int   jactive           = NG;
            PLOOP prim[i][j][k][ip] = prim[i][jactive][k][ip];
            pflag[i][j][k]          = pflag[i][jactive][k];
            #elif X2L_GAS_BOUND == BC_POLAR
            int   jactive           = -j + 2 * NG - 1;
            PLOOP prim[i][j][k][ip] = prim[i][jactive][k][ip];
            pflag[i][j][k]          = pflag[i][jactive][k];
            prim[i][j][k][U2] *= -1.;
            prim[i][j][k][B2] *= -1.;
            #elif X2L_GAS_BOUND == BC_PROB
            printf("X2L_GAS_BOUND choice BC_PROB not supported\n");
            exit(-1);
            #elif X2L_GAS_BOUND != BC_PERIODIC
            printf("X2L_GAS_BOUND choice %i not supported\n", X2L_GAS_BOUND);
            exit(-1);
            #endif // X2L_GAS_BOUND 
          }
          #endif // N2 < NG
        }
      }
    }
  } // global_start[2] == 0

  /* fill ghost zone across MPI processes */
  sync_mpi_boundaries_X2L(prim);

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  // globally outer-y boundary //  
  if (global_stop[2] == N2TOT) {
#pragma omp parallel for collapse(2)
    ILOOPALL {
      KSLOOP(0, N3 - 1) {
        JSLOOP(N2, N2 - 1 + NG) {
          #if N2 < NG
          int   jactive           = N2 - 1 + NG;
          PLOOP prim[i][j][k][ip] = prim[i][jactive][k][ip];
          pflag[i][j][k]          = pflag[i][jactive][k];
          #else // N2 < NG
          {
            #if X2R_GAS_BOUND == BC_OUTFLOW
            int   jactive           = N2 - 1 + NG;
            PLOOP prim[i][j][k][ip] = prim[i][jactive][k][ip];
            pflag[i][j][k]          = pflag[i][jactive][k];
            #elif X2R_GAS_BOUND == BC_POLAR
            int   jactive           = -j + 2 * (N2 + NG) - 1;
            PLOOP prim[i][j][k][ip] = prim[i][jactive][k][ip];
            pflag[i][j][k]          = pflag[i][jactive][k];
            prim[i][j][k][U2] *= -1.;
            prim[i][j][k][B2] *= -1.;
            #elif X2R_GAS_BOUND == BC_PROB
            printf("X2R_GAS_BOUND choice BC_PROB not supported\n");
            exit(-1);
            #elif X2R_GAS_BOUND != BC_PERIODIC
            printf("X2R_GAS_BOUND choice %i not supported\n", X2R_GAS_BOUND);
            exit(-1);
            #endif // X2R_GAS_BOUND
          }
          #endif // N2 < NG
        }
      }
    }
  } // global_stop[2] == N2TOT

  /* fill ghost zone across MPI processes */
  sync_mpi_boundaries_X2R(prim);

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  // globally inner-z boundary // 
  if (global_start[3] == 0) {
#pragma omp parallel for collapse(2)
    ILOOPALL {
      JLOOPALL {
        KSLOOP(-NG, -1) {
          #if N3 < NG
          int   kactive           = NG;
          PLOOP prim[i][j][k][ip] = prim[i][j][kactive][ip];
          pflag[i][j][k]          = pflag[i][j][kactive];
          #else // N3 < NG
          {
            #if X3L_GAS_BOUND == BC_OUTFLOW
            int   kactive           = NG;
            PLOOP prim[i][j][k][ip] = prim[i][j][kactive][ip];
            pflag[i][j][k]          = pflag[i][j][kactive];
            #elif X3L_GAS_BOUND == BC_PROB
            printf("X3L_GAS_BOUND choice BC_PROB not supported\n");
            exit(-1);
            #elif X3L_GAS_BOUND != BC_PERIODIC
            printf("X3L_GAS_BOUND choice %i not supported\n", X3L_GAS_BOUND);
            exit(-1);
            #endif // X3L_GAS_BOUND
          }
        #endif // N3 < NG
        }
      }
    }
  } // global_start[3] == 0

  /* fill ghost zone across MPI processes */
  sync_mpi_boundaries_X3L(prim);

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  // globally outer-z boundary // 
  if (global_stop[3] == N3TOT) {
#pragma omp parallel for collapse(2)
    ISLOOP(-NG, N1 - 1 + NG) {
      JSLOOP(-NG, N2 - 1 + NG) {
        KSLOOP(N3, N3 - 1 + NG) {
          #if N3 < NG
          int   kactive           = N3 - 1 + NG;
          PLOOP prim[i][j][k][ip] = prim[i][j][kactive][ip];
          pflag[i][j][k]          = pflag[i][j][kactive];
          #else // N3 < NG
          {
            #if X3R_GAS_BOUND == BC_OUTFLOW
            int   kactive           = NG;
            PLOOP prim[i][j][k][ip] = prim[i][j][kactive][ip];
            pflag[i][j][k]          = pflag[i][j][kactive];
            #elif X3R_GAS_BOUND == BC_PROB
            printf("X3R_GAS_BOUND choice BC_PROB not supported\n");
            exit(-1);
            #elif X3R_GAS_BOUND != BC_PERIODIC
            printf("X3R_GAS_BOUND choice %i not supported\n", X3R_GAS_BOUND);
            exit(-1);
            #endif // X3R_GAS_BOUND
          }
        #endif // N3 < NG
        }
      }
    }
  } // global_stop[3] == N3TOT

  /* fill ghost zone across MPI processes */
  sync_mpi_boundaries_X3R(prim);

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

  #if METRIC == MKS
  // Make sure there is no inflow at the inner boundary // 
  if (global_start[1] == 0 && X1L_INFLOW == 0) {
#pragma omp parallel for collapse(2)
    ISLOOP(-NG, -1) {
      JLOOPALL {
        KLOOPALL { inflow_check(prim[i][j][k], i, j, 0); }
      }
    }
  }

  // Make sure there is no inflow at the outer boundary //
  if (global_stop[1] == N1TOT && X1R_INFLOW == 0) {
#pragma omp parallel for collapse(2)
    ISLOOP(N1, N1 - 1 + NG) {
      JLOOPALL {
        KLOOPALL { inflow_check(prim[i][j][k], i, j, 1); }
      }
    }
  }
  #endif

  // stop the timer //
  timer_stop(TIMER_BOUND);

}

#if METRIC == MKS
/************************************************************************/
/* make sure no inflow at inner/outer boundaries */
void inflow_check(double *Pr, int ii, int jj, int type) {

  // variables //
  struct of_geom *geom;
  double          ucon[NDIM];
  double          alpha, beta1, gamma, vsq;

  // metric and 4-velocity //
  geom = get_geometry(ii, jj, 0, CENT);
  ucon_calc(Pr, geom, ucon);

  // check for u1 > 0 or u1 < 0 //
  if (((ucon[1] > 0.) && (type == 0)) || ((ucon[1] < 0.) && (type == 1))) {
    // Find gamma and remove it from primitives //
    if (mhd_gamma_calc(Pr, geom, &gamma)) {
      fprintf(stderr, "\ninflow_check(): gamma failure\n");
      fail(FAIL_GAMMA);
    }
    Pr[U1] /= gamma;
    Pr[U2] /= gamma;
    Pr[U3] /= gamma;
    alpha = geom->alpha;
    beta1 = geom->gcon[0][1] * alpha * alpha;

    // Reset radial velocity so radial 4-velocity is zero //
    Pr[U1] = beta1 / alpha;

    // Now find new gamma and put it back in //
    vsq = 0.;
    for (int mu = 1; mu < NDIM; mu++) {
      for (int nu = 1; nu < NDIM; nu++) {
        vsq += geom->gcov[mu][nu] * Pr[U1 + mu - 1] * Pr[U1 + nu - 1];
      }
    }
    if (fabs(vsq) < 1.e-13)
      vsq = 1.e-13;
    if (vsq >= 1.) {
      vsq = 1. - 1. / (GAMMAMAX * GAMMAMAX);
    }
    gamma = 1. / sqrt(1. - vsq);
    Pr[U1] *= gamma;
    Pr[U2] *= gamma;
    Pr[U3] *= gamma;
  }
}

/************************************************************************/
/* fix the Riemann flux at zone boundaries to prevent inflow */
void fix_flux(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3) {

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* globally inner-x boundary */
  if (global_start[1] == 0 && X1L_INFLOW == 0) {
#pragma omp parallel for simd collapse(2)
    JLOOPALL {
      KLOOPALL {
        //////////////////////
        // JSLOOP(0, N2-1) {
        //  KSLOOP(0, N3-1) {
        //////////////////////
        F1[0 + NG][j][k][RHO] = MY_MIN(F1[0 + NG][j][k][RHO], 0.);
      }
    }
  }

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* globally outer-x boundary */
  if (global_stop[1] == N1TOT && X1R_INFLOW == 0) {
#pragma omp parallel for simd collapse(2)
    JLOOPALL {
      KLOOPALL {
        //////////////////////
        // JSLOOP(0, N2-1) {
        //  KSLOOP(0, N3-1) {
        //////////////////////
        F1[N1 + NG][j][k][RHO] = MY_MAX(F1[N1 + NG][j][k][RHO], 0.);
      }
    }
  }

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* globally inner-y boundary */
  if (global_start[2] == 0 && X2L_INFLOW == 0) {
#pragma omp parallel for simd collapse(2)
    ILOOPALL {
      KLOOPALL {
        //////////////////////
        // ISLOOP(0, N1-1) {
        //  KSLOOP(0, N3-1) {
        //////////////////////
        F1[i][-1 + NG][k][B2]      = -F1[i][0 + NG][k][B2];
        F3[i][-1 + NG][k][B2]      = -F3[i][0 + NG][k][B2];
        PLOOP F2[i][0 + NG][k][ip] = 0.;
      }
    }
  }

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* globally outer-y boundary */
  if (global_stop[2] == N2TOT && X2R_INFLOW == 0) {
#pragma omp parallel for simd collapse(2)
    ILOOPALL {
      KLOOPALL {
        //////////////////////
        // ISLOOP(0, N1-1) {
        //  KSLOOP(0, N3-1) {
        //////////////////////
        F1[i][N2 + NG][k][B2]       = -F1[i][N2 - 1 + NG][k][B2];
        F3[i][N2 + NG][k][B2]       = -F3[i][N2 - 1 + NG][k][B2];
        PLOOP F2[i][N2 + NG][k][ip] = 0.;
      }
    }
  }
}
#endif // METRIC

#if RADIATION
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

/************************************************************************/
/* apply boundary conditions for superphotons */
void bound_superphotons(grid_prim_type P, double t, double dt) {

  // start timer //
  timer_start(TIMER_BOUND);

  // initialize //
  int step_lost_local  = 0;
  int step_tot_local   = 0;
  int tracer_tot_local = 0;

  ///////////////////////////////////////////////////////////////////
  /*double ephot_beg = 0.;
  int ns_beg = 0;
  #pragma omp parallel reduction(+:ephot_beg) reduction(+:ns_beg)
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    while (ph != NULL) {
      ephot_beg -= ph->w*kphys_to_num*ph->Kcov[2][0];
      if (ph->w > 0.)
        ns_beg += 1;
      ph = ph->next;
    }
  }*/
  ///////////////////////////////////////////////////////////////////

  // initialize //
  int n_to_send = 0;

#pragma omp parallel reduction(+:step_lost_local) reduction(+:step_tot_local) \
                     reduction(+:n_to_send) reduction(+:tracer_tot_local)
  {

    // variables //
    double            X[NDIM], Kcov[NDIM], Kcon[NDIM];
    struct of_photon *ph, *head, *prev;

    // superphoton list //
    ph   = photon_lists[omp_get_thread_num()];
    prev = NULL;
    head = ph;

    // loop over non null photon list //
    while (ph != NULL) {

      // Use interpolation to get X^{\mu}, K_{\mu} at time t + dt //
      int status = get_X_K_interp(ph, t + dt, P, X, Kcov, Kcon);

      // Test whether superphoton inactive at t = t + dt //
      int active_bc = bound_rad_isactive(X, ph);

      // transport superphotons across boundaries //
      bound_rad_transport(X, ph, 0);

      // check for error //
      int active_error = rad_error_check(&ph);

      // Update global log variables //
      //////////////////////////
      // rad_boundary_log(ph);
      //////////////////////////
      if (!active_bc || !active_error || status == SPH_INTERP_FAIL) {
        #if RADIATION == RADTYPE_NEUTRINOS
        /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
        record_lepton_flux(ph);
        /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
        #endif
        step_lost_local++;
        step_tot_local--;
        if (ph->type == TYPE_TRACER)
          tracer_tot_local--;
        list_remove(&ph, &head, &prev);
        continue;
      }

      // record mpi transport //
      if (rad_mpi_transport(&ph, &head, &prev, X, active_bc)) {
        n_to_send++;
        continue;
      }

      // continue to the next photon //
      prev = ph;
      ph   = ph->next;

    }

    // photon list //
    photon_lists[omp_get_thread_num()] = head;

  } // omp parallel

  // update //
  step_lost += step_lost_local;
  step_tot += step_tot_local;
  tracer_tot += tracer_tot_local;

  // Reduce MPI-ready superphotons across threads //
  struct of_photon *ph_mpi = NULL, *ph = NULL;
  for (int n = 0; n < nthreads; n++) {
    if (photon_mpi_lists[n] != NULL) {
      ph = photon_mpi_lists[n];
      while (ph != NULL) {
        swap_ph(&ph, &ph_mpi);
      }
      photon_mpi_lists[n] = NULL;
    }
  }

  // synchronization //
  sync_mpi_photons(&ph_mpi, P, t, dt);
  sync_radG();
  sync_Jrad();
  sync_radtype_vec(Nem_phys);
  sync_radtype_vec(Nabs_phys);

  // stop timer //
  timer_stop(TIMER_BOUND);

}

/************************************************************************/
/* check for superphoton error */
int rad_error_check(struct of_photon **ph) {
  int active = 1;
  if ((*ph)->w < SMALL || (*ph)->Kcov[2][0] > 0.)
    active = 0;
  return active;
}

/************************************************************************/
// Return 1 if superphoton should be communicated over MPI //
int rad_mpi_transport(struct of_photon **ph, struct of_photon **head,
                      struct of_photon **prev, double X[NDIM], int active) {

  // Only transport superphotons active at X[0] = t + dt //
  if (!active)
    return 0;

  // check if the photon is within the domain //
  if (((X[1] < startx_proc[1] || X[1] > stopx_proc[1]) && N1CPU > 1) ||
      ((X[2] < startx_proc[2] || X[2] > stopx_proc[2]) && N2CPU > 1) ||
      ((X[3] < startx_proc[3] || X[3] > stopx_proc[3]) && N3CPU > 1)) {
    struct of_photon *next = (*ph)->next;

    (*ph)->next = photon_mpi_lists[omp_get_thread_num()];
    photon_mpi_lists[omp_get_thread_num()] = *ph;

    if (*prev != NULL) {
      (*prev)->next = next;
      *ph           = (*prev)->next;
    } else {
      (*head) = next;
      *ph     = *head;
    }

    return 1;
  } else {
    return 0;
  }

}

/************************************************************************/
/* define constant used below */
#define DX2MIN (0.05)

/************************************************************************/
// Replace superphoton with one headed away from polar boundary //
void polar_fix(double X[NDIM], struct of_photon *ph) {
  #if METRIC == MKS
  if (X[2] < DX2MIN || X[2] > 1. - DX2MIN) {
    printf("POLAR FIXUP!\n");
    printf("BEFORE:\n");
    for (int mu = 0; mu < NDIM; mu++) {
      printf("[%i] X[][] = %e %e %e\n", mu, ph->X[0][mu], ph->X[1][mu], ph->X[2][mu]);
      printf("[%i] Kcon[][] = %e %e %e\n", mu, ph->Kcon[0][mu], ph->Kcon[1][mu], ph->Kcon[2][mu]);
      printf("[%i] Kcov[][] = %e %e %e\n", mu, ph->Kcov[0][mu], ph->Kcov[1][mu], ph->Kcov[2][mu]);
      printf("KDOTK = %e\n", ph->Kcon[2][0] * ph->Kcov[2][0] +
                             ph->Kcon[2][1] * ph->Kcov[2][1] +
                             ph->Kcon[2][2] * ph->Kcov[2][2] +
                             ph->Kcon[2][3] * ph->Kcov[2][3]);
    }

    ph->t0 = ph->X[1][0];
    for (int mu = 1; mu < NDIM; mu++) {
      ph->Kcon[1][mu] *= -1.;
      ph->Kcov[1][mu] *= -1.;
    }
    for (int mu = 0; mu < NDIM; mu++) {
      ph->X[2][mu]    = ph->X[1][mu];
      ph->Kcon[2][mu] = ph->Kcon[1][mu];
      ph->Kcov[2][mu] = ph->Kcov[1][mu];
      ph->X[0][mu]    = ph->X[1][mu];
      ph->Kcon[0][mu] = ph->Kcon[1][mu];
      ph->Kcov[0][mu] = ph->Kcov[1][mu];
    }

    printf("AFTER:\n");
    for (int mu = 0; mu < NDIM; mu++) {
      printf("[%i] X[][] = %e %e %e\n", mu, ph->X[0][mu], ph->X[1][mu], ph->X[2][mu]);
      printf("[%i] Kcon[][] = %e %e %e\n", mu, ph->Kcon[0][mu], ph->Kcon[1][mu], ph->Kcon[2][mu]);
      printf("[%i] Kcov[][] = %e %e %e\n", mu, ph->Kcov[0][mu], ph->Kcov[1][mu], ph->Kcov[2][mu]);
      printf("KDOTK = %e\n", ph->Kcon[2][0] * ph->Kcov[2][0] +
                             ph->Kcon[2][1] * ph->Kcov[2][1] +
                             ph->Kcon[2][2] * ph->Kcov[2][2] +
                             ph->Kcon[2][3] * ph->Kcov[2][3]);
    }
  }
#endif // METRIC
}

/************************************************************************/
/* define constant */
#undef DX2MIN

/************************************************************************/
// If X[] (at t+dt) calls for boundary transport, modify both current and
// previous X^{\mu}, K_{\mu} //
void bound_rad_transport(double X[NDIM], struct of_photon *ph, int is_transported) {

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* inner-x boundary */
  if (X[1] < startx_rad[1]) {
    #if X1L_RAD_BOUND == BC_PERIODIC
    if (N1CPU == 1 || is_transported) {
      for (int n = 0; n < NSUP; n++) {
        ph->X[n][1] = stopx_rad[1] - (startx_rad[1] - ph->X[n][1]);
      }
    }
    #endif
  }

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* outer-x boundary */
  if (X[1] > stopx_rad[1]) {
    #if X1R_RAD_BOUND == BC_PERIODIC
    if (N1CPU == 1 || is_transported) {
      for (int n = 0; n < NSUP; n++) {
        ph->X[n][1] = startx_rad[1] + (ph->X[n][1] - stopx_rad[1]);
      }
    }
    #endif
  }

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* inner-y boundary */
  if (X[2] < startx_rad[2]) {
    #if X2L_RAD_BOUND == BC_PERIODIC
    if (N2CPU == 1 || is_transported) {
      for (int n = 0; n < NSUP; n++) {
        ph->X[n][2] = stopx_rad[2] - (startx_rad[2] - ph->X[n][2]);
      }
    }
      #endif
  }

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* outer-y boundary */
  if (X[2] > stopx_rad[2]) {
    #if X2R_RAD_BOUND == BC_PERIODIC
    if (N2CPU == 1 || is_transported) {
      for (int n = 0; n < NSUP; n++) {
        ph->X[n][2] = startx_rad[2] + (ph->X[n][2] - stopx_rad[2]);
      }
    }
    #endif
  }

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* inner-z boundary */
  if (X[3] < startx_rad[3]) {
    #if X3L_RAD_BOUND == BC_PERIODIC
    if (N3CPU == 1 || is_transported) {
      for (int n = 0; n < NSUP; n++) {
        ph->X[n][3] = stopx_rad[3] - (startx_rad[3] - ph->X[n][3]);
      }
    }
    #endif
  }

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* outer-z boundary */
  if (X[3] > stopx_rad[3]) {
    #if X3R_RAD_BOUND == BC_PERIODIC
    if (N3CPU == 1 || is_transported) {
      for (int n = 0; n < NSUP; n++) {
        ph->X[n][3] = startx_rad[3] + (ph->X[n][3] - stopx_rad[3]);
      }
    }
    #endif
  }

}

/************************************************************************/
/* check if a photon is active */
int bound_rad_isactive(double X[NDIM], struct of_photon *ph) {

  // initialize //
  int active = 1;

  // Polar fix, to prevent challenging MPI transports and geodesic explosions //
  #if METRIC == MKS
  double th = th_of_X(X) * 180. / M_PI;
  /////////////////////////////////////////////////
  // if (N3TOT > 1 && (th < 3. || th > 177.)) {
  /////////////////////////////////////////////////
  if (th < 3. || th > 177.) {
    active = 0;
  }
  #endif

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* inner-x boundary */
  if (X[1] < startx_rad[1]) {
    #if X1L_RAD_BOUND == BC_PERIODIC
    #elif X1L_RAD_BOUND == BC_ESCAPE
    active = 0;
    #elif X1L_RAD_BOUND == BC_EQUILIB
    printf("X1L_RAD_BOUND BC_EQUILIB not supported\n");
    exit(-1);
    #elif X1L_RAD_BOUND == BC_CAMERA
    printf("X1L_RAD_BOUND BC_CAMERA not supported\n");
    exit(-1);
    #else
    printf("X1L_RAD_BOUND %i not supported\n", X1L_RAD_BOUND);
    exit(-1);
    #endif
  }

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* outer-x boundary */
  if (X[1] > stopx_rad[1] && ph->type != TYPE_TRACER) {
    #if X1R_RAD_BOUND == BC_PERIODIC
    #elif X1R_RAD_BOUND == BC_ESCAPE
    active = 0;
    #elif X1R_RAD_BOUND == BC_REQUILIB
    printf("X1R_RAD_BOUND BC_EQUILIB not supported\n");
    exit(-1);
    #elif X1R_RAD_BOUND == BC_CAMERA
    #if METRIC == MINKOWSKI
    printf("BC_CAMERA not supported with MINKOWSKI METRIC\n");
    exit(-1);
    #endif
    record_superphoton(X, ph);
    active = 0;
    // record_super_photon -- put energy into relevant bins, dont worry about //
    // destroying superphoton //
    ///////////////////////////////////////////////////////
    // printf("X1R_RAD_BOUND BC_CAMERA not supported\n");
    // exit(-1);
    ///////////////////////////////////////////////////////
    #else
    printf("X1R_RAD_BOUND %i not supported\n", X1R_RAD_BOUND);
    exit(-1);
    #endif
  }

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  // Special case for tracers since their boundary is at
  // Rout, not Rout_rad
  if (X[1] > stopx[1] && ph->type == TYPE_TRACER) {
    #if X1R_RAD_BOUND == BC_PERIODIC
    #elif X1R_RAD_BOUND == BC_ESCAPE
    active = 0;
    #elif X1R_RAD_BOUND == BC_REQUILIB
    printf("X1R_RAD_BOUND BC_EQUILIB not supported\n");
    exit(-1);
    #elif X1R_RAD_BOUND == BC_CAMERA
    #if METRIC == MINKOWSKI
    printf("BC_CAMERA not supported with MINKOWSKI METRIC\n");
    exit(-1);
    #endif
    active = 0;
    #else
    printf("X1R_RAD_BOUND %i not supported\n", X1R_RAD_BOUND);
    exit(-1);
    #endif
  }

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* inner-y boundary */
  if (X[2] < startx_rad[2]) {
    #if X2L_RAD_BOUND == BC_PERIODIC
    #elif X2L_RAD_BOUND == BC_ESCAPE
    active = 0;
    #elif X2L_RAD_BOUND == BC_EQUILIB
    printf("X2L_RAD_BOUND BC_EQUILIB not supported\n");
    exit(-1);
    #elif X2L_RAD_BOUND == BC_CAMERA
    printf("X2L_RAD_BOUND BC_CAMERA not supported\n");
    exit(-1);
    #else
    printf("X2L_RAD_BOUND %i not supported\n", X2L_RAD_BOUND);
    exit(-1);
    #endif
  }

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* outer-y boundary */
  if (X[2] > stopx_rad[2]) {
    #if X2R_RAD_BOUND == BC_PERIODIC
    #elif X2R_RAD_BOUND == BC_ESCAPE
    active = 0;
    #elif X2R_RAD_BOUND == BC_EQUILIB
    printf("X2R_RAD_BOUND BC_EQUILIB not supported\n");
    exit(-1);
    #elif X2R_RAD_BOUND == BC_CAMERA
    printf("X2R_RAD_BOUND BC_CAMERA not supported\n");
    exit(-1);
    #else
    printf("X2R_RAD_BOUND %i not supported\n", X2R_RAD_BOUND);
    exit(-1);
    #endif
  }

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* inner-z boundary */
  if (X[3] < startx_rad[3]) {
    #if X3L_RAD_BOUND == BC_PERIODIC
    #elif X3L_RAD_BOUND == BC_ESCAPE
    active = 0;
    #elif X3L_RAD_BOUND == BC_EQUILIB
    printf("X3L_RAD_BOUND BC_EQUILIB not supported\n");
    exit(-1);
    #elif X3L_RAD_BOUND == BC_CAMERA
    printf("X3L_RAD_BOUND BC_CAMERA not supported\n");
    exit(-1);
    #else
    printf("X3L_RAD_BOUND %i not supported\n", X3L_RAD_BOUND);
    exit(-1);
    #endif
  }

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* outer-z boundary */
  if (X[3] > stopx_rad[3]) {
    #if X3R_RAD_BOUND == BC_PERIODIC
    #elif X3R_RAD_BOUND == BC_ESCAPE
    active = 0;
    #elif X3R_RAD_BOUND == BC_EQUILIB
    printf("X3R_RAD_BOUND BC_EQUILIB not supported\n");
    exit(-1);
    #elif X3R_RAD_BOUND == BC_CAMERA
    printf("X3R_RAD_BOUND BC_CAMERA not supported\n");
    exit(-1);
    #else
    printf("X3R_RAD_BOUND %i not supported\n", X3R_RAD_BOUND);
    exit(-1);
    #endif
  }

  // return //
  return active;
  
}

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
#endif // RADIATION
