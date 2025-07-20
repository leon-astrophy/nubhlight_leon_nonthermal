/******************************************************************************
 *                                                                            *
 * STEP.C                                                                     *
 *                                                                            *
 * ADVANCES FLUID QUANTITIES BY ONE TIMESTEP                                  *
 *                                                                            *
 ******************************************************************************/

// headers //
#include "decs.h"

// define function used here // 
double advance(grid_prim_type Pi, grid_prim_type Pb, double Dt, grid_prim_type Pf, int stage);
double fluxcalc(grid_prim_type Pr);
void   flux_ct(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3);
void   lr_to_flux(double p_l[NVAR], double p_r[NVAR], struct of_geom *geom, int dir, double Flux[NVAR], double *maxSpeed);
void apply_perturbation(grid_prim_type Pf);
#if RADIATION
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
void apply_rad_force(grid_prim_type Pr, double Dt);
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
#endif

// DEBUG, record minimum dt //
#if RECORD_DT_MIN
double dt_min = INFINITY;
#endif

/*#################################################################################################*/

/*************************************************************************/
/* march foward in RK substep */
void step() {

  // time step //
  double ndt;
  #if RADIATION
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  double dt_cool;
  #if RADIATION == RADTYPE_NEUTRINOS && LOCAL_ANGULAR_DISTRIBUTIONS && RAD_NUM_TYPES >= 4 && NEUTRINO_OSCILLATIONS
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  // not used for timestep control. Used to turn oscillations on or off //
  double dt_osc = get_dt_oscillations();
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  #endif // oscillations
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  #endif // radiation

  // backup time step //
  dtsave = dt;

  // Need both P_n and P_n+1 to calculate current
  // Leon: added OMP here //
#pragma omp parallel for simd collapse(4)
  ZSLOOP(-NG, N1 - 1 + NG, -NG, N2 - 1 + NG, -NG, N3 - 1 + NG) {
    PLOOP { Psave[i][j][k][ip] = P[i][j][k][ip]; }
  }

  // Estimate electron temperature //
  #if RADIATION && ESTIMATE_THETAE
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  estimate_Thetae(P, extra, t, dt);
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  #endif

  /* perform perturbation */
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  #if NONTHERMAL && PERTURB_NTE
  apply_perturbation(P);
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  #endif

  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  // Predictor step //
  ndt = advance(P, P, 0.5 * dt, Ph, 0);

  /* Now that we have primitive Ph after explicit step */
  /* perform operator split (simple split) */
  #if NONTHERMAL && !(SKIP_ADIAB)
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* First do adiabatic step for nonthermal */
  nonthermal_adiab(P, P, Ph, 0.5*dt);
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  #endif
  #if NONTHERMAL && CONST_INJECTION
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* then, for constant injection model*/
  const_inject(P, P, Ph, 0.5*dt);
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  #endif
  /* then viscous heat */
  #if ELECTRONS
  heat_electrons(P, P, Ph, 0.5 * dt);
  #endif
  #if NONTHERMAL && !(SKIP_COOLING)
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* finally, apply cooling to nonthermal */
  cool_nonthermal(P, P, Ph, 0.5*dt);
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  #endif

  /* fixup */
  fixup(Ph, extra);
  fixup_utoprim(Ph, extra);
  #if ELECTRONS
  fixup_electrons(Ph);
  #if COULOMB
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  coulomb(P, P, Ph, 0.5 * dt);
  fixup_electrons(Ph);
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  #endif // COULOMB
  #endif // ELECTRONS

  /* boundary condition */
  bound_prim(Ph);

  // Radiation step //
  #if RADIATION
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  // get fluid qunatities //
  precompute_microphysics();
  // emit superphoton //
  make_superphotons(Ph, extra, t, dt);
  /////////////////////////////////////////
  // check_nu_type("after make"); // DEBUG
  /////////////////////////////////////////
  // evolve superphoton along geodesic //
  push_superphotons(P, Ph, dt);
  /////////////////////////////////////////
  // check_nu_type("after push"); // DEBUG
  /////////////////////////////////////////
  // perform scattering/absorption //
  interact(Ph, extra, t, dt);
  /////////////////////////////////////////////
  // check_nu_type("after interact"); // DEBUG
  /////////////////////////////////////////////
  // boundary condition //
  bound_superphotons(Ph, t, dt);
  ///////////////////////////////////////////
  // check_nu_type("after bound"); // DEBUG
  ///////////////////////////////////////////
  #if RADIATION == RADTYPE_NEUTRINOS && LOCAL_ANGULAR_DISTRIBUTIONS && RAD_NUM_TYPES >= 4 && NEUTRINO_OSCILLATIONS
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  int oscillations_active = (dt_osc <= dt);
  if (mpi_io_proc()) {
    printf("\t[Oscillations] Active? %d dt_osc = %.14e\n", oscillations_active, dt_osc);
  }
  if (oscillations_active) { // TOOD(JMM): Some safety factor?
    accumulate_local_angles();
    oscillate(local_moments, Gnu);
  }
  ///////////////////////////////////////////////////
  // check_nu_type("after oscillate"); // DEBUG
  ///////////////////////////////////////////////////
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  #endif // OSCILLATIONS
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  #endif

  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/ 
  // Corrector step //
  ndt = advance(P, Ph, dt, P, 1);

  /* Now that we have primitive P after explicit step */
  /* perform operator split (simple split) */
  #if NONTHERMAL && !(SKIP_ADIAB)
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* First do adiabatic step for nonthermal */
  /* Leon: shouldn't Psave be Ph? */
  nonthermal_adiab(P, Ph, P, dt);
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  #endif
  #if NONTHERMAL && CONST_INJECTION
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  //* then, for constant injection model*/
  const_inject(P, Ph, P, dt);
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  #endif
  /* then viscous heat */
  #if ELECTRONS
  heat_electrons(P, Ph, P, dt);
  #endif
  #if NONTHERMAL && !(SKIP_COOLING)
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* finally, apply cooling to nonthermal */
  cool_nonthermal(P, Ph, P, dt);
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  #endif

  /* fixup */
  fixup(P, extra);
  fixup_utoprim(P, extra);
  #if ELECTRONS
  fixup_electrons(P);
  #if COULOMB
  /*@@@@@@@@@@@@@@@@@@@@@@@@*/
  coulomb(P, Ph, P, dt);
  fixup_electrons(P);
  /*@@@@@@@@@@@@@@@@@@@@@@@@*/
  #endif // COULOMB
  #endif // ELECTRONS

  /* boundary condition */
  bound_prim(P);

  #if RADIATION
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  // Apply radiation four-force to fluid //
  apply_rad_force(P, dt);
  fixup(P, extra);
  fixup_utoprim(P, extra);
  #if ELECTRONS
  apply_rad_force_e(Ph, P, radG, dt);
  fixup_electrons(P);
  #endif // ELECTRONS 
  bound_prim(P);

  // Reset radG by memset, much master //
  memset((void *)&radG[0][0][0][0], 0, (N1 + 2 * NG) * (N2 + 2 * NG) * (N3 + 2 * NG) * (NDIM + NRADCOMP) * sizeof(double));

  #if TRACERS
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  prune_tracers();
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  #endif // TRACERS
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  #endif // RADIATION

  // Reset fixup mask // 
  memset((void *)&fixup_required[0][0][0], 0, (N1 + 2 * NG) * (N2 + 2 * NG) * (N3 + 2 * NG) * sizeof(int));

  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

  // Increment time //
  t += dt;

  // Get dt for the next step //
  #if RADIATION
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  dt_cool = get_min_dt_cool(P, extra);
  ndt     = MY_MIN(cour * dt_light_min, cour_cool * dt_cool);
  #if RECORD_DT_MIN
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  if (ndt < dt_min)
    dt_min = ndt;
  // DEBUG //
  printf("dt_cool       = %g\n"
         "dt_light_min  = %g\n"
         "dt_global_min = %g\n",
          dt_cool, dt_light_min, dt_min);
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  #endif // RECORD_DT_MIN
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  #endif // RADIATION

  // Limit the  next timestep //
  if (ndt > SAFE * dt) {
    ndt = SAFE * dt;
  }

  //////////////////////////
  // ndt = 0.01; // DEBUG
  //////////////////////////
  dt = ndt;

  /////////////
  // dt = ndt;
  /////////////
  // MPI synchronize //
  dt = mpi_min(dt);

  // if dt is too small, try to advance anyway. Also complain. //
  #if RADIATION
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  if (dt < SMALL) { // || dt_cool/dt_light_min < 1e-2) {
    fprintf(stderr,
        "ERROR: dt too smalL! Trying to advance anyway.\n"
        "\tdt     = %g\n"
        "\ttrying = %g\n",
        dt, SMALL + dtsave / SAFE);
    dt = SMALL + dtsave / SAFE;
  }
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  #endif

  // Don't step beyond end of run //
  if (t + dt > tf) {
    dt = tf - t;
  }

}

/**********************************************************************************/
/* march foward in time in a RK sub-step */
double advance(grid_prim_type Pi, grid_prim_type Pb, double Dt, grid_prim_type Pf, int stage) {

  // variables //
  double ndt, U[NVAR], dU[NVAR];
  struct of_state qi;

  // backup the primitives //
  #pragma omp parallel for simd collapse(4)
  ZLOOP PLOOP Pf[i][j][k][ip] = Pi[i][j][k][ip];

  // set timer //
  timer_start(TIMER_FLUXCALC);

  // get flux and timestep //
  ndt = fluxcalc(Pb);

  // fix the flux at the inner boundary //
  #if METRIC == MKS
  fix_flux(F1, F2, F3);
  #endif

  // flux-ct scheme to preserve divB = 0 //
  flux_ct(F1, F2, F3);

  // stop the timer // 
  timer_stop(TIMER_FLUXCALC);

  // Evaluate diagnostics based on fluxes //
  timer_start(TIMER_DIAG);
  diag_flux(F1, F2, F3);
  timer_stop(TIMER_DIAG);

  // trivial exit with no update //
  #if NO_GRMHD_UPDATE
  #pragma omp parallel for simd collapse(3)
  ZLOOP pflag[i][j][k] = 0;
  return ndt;
  #endif

  // Update Pi to Pf //
  timer_start(TIMER_UPDATE);
  #pragma omp parallel for private(dU, qi, U) collapse(3) // schedule(guided)
  ZLOOP {

    // source term and Riemann fluxes //
    source(Pb[i][j][k], &(ggeom[i][j][CENT]), i, j, dU, Dt, extra[i][j][k]);
    get_state(Pi[i][j][k], &(ggeom[i][j][CENT]), &qi);
    primtoflux(Pi[i][j][k], &qi, 0, 0, &(ggeom[i][j][CENT]), U);

    // loop over equations //
    PLOOP {
      U[ip] += Dt * (-(F1[i + 1][j][k][ip] - F1[i][j][k][ip]) / dx[1] -
                      (F2[i][j + 1][k][ip] - F2[i][j][k][ip]) / dx[2] -
                      (F3[i][j][k + 1][ip] - F3[i][j][k][ip]) / dx[3] + dU[ip]);
    }

    // convert conservative to primitive variables //  
    pflag[i][j][k] = Utoprim(U, &(ggeom[i][j][CENT]), Pf[i][j][k]);

    // record failure mode //
    if (pflag[i][j][k])
      fail_save[i][j][k] = 1;

  } // ZLOOP

  // set timer // 
  timer_stop(TIMER_UPDATE);

  // return with time step //
  return (ndt);

}

#if RADIATION
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

/*****************************************************************************/
/* Apply radiation 4-force */
void apply_rad_force(grid_prim_type Pr, double Dt) {

  // variables //
  double          U[NVAR];
  struct of_state q;

  // start the timer //
  timer_start(TIMER_UPDATE);

  // loop over //
#pragma omp parallel for private(q, U) collapse(3) // schedule(guided)
  ZLOOP {

    /////////////////////////////////////////////////////////////////
    // DEBUG
    /*
    double zonevol = dV*L_unit*L_unit*L_unit*ggeom[i][j][CENT].g;
    printf("Dt = %g\n",Dt);
    printf("dU/dt = %g\n",Dt*radG[i][j][k][0]*U_unit);
    printf("dE/dt = %g\n",zonevol*Dt*radG[i][j][k][0]*U_unit);
    */
    /////////////////////////////////////////////////////////////////

    // Store primitive variables before cooling for supercooling diagnostics //
    PLOOP psupersave[i][j][k][ip] = Pr[i][j][k][ip];

    // get conservative variables //
    get_state(Pr[i][j][k], &(ggeom[i][j][CENT]), &q);
    primtoflux(Pr[i][j][k], &q, 0, 0, &(ggeom[i][j][CENT]), U);

    // update conservative variables //
    for (int ip = 1; ip < 5; ip++) {
      U[ip] += Dt * radG[i][j][k][ip - 1];
    }

    // update 4-force //
    DLOOP1 { radG_int[i][j][k][mu] += Dt * radG[i][j][k][mu]; }

    // TODO generalize this if need be //
    #if RADIATION == RADTYPE_NEUTRINOS
    /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
    {
      U[YE] += Dt * radG[i][j][k][RADG_YE];
      U[YE_EM] += Dt * radG[i][j][k][RADG_YE_EM];
      radG_int[i][j][k][RADG_YE] += Dt * radG[i][j][k][RADG_YE];
      radG_int[i][j][k][RADG_YE_EM] += Dt * radG[i][j][k][RADG_YE_EM];
    }
    /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
    #endif

    // convert from conservative to primitive //
    pflag[i][j][k] = Utoprim(U, &(ggeom[i][j][CENT]), Pr[i][j][k]);

    ///////////////////////////////////////
    /*if (pflag[i][j][k] == 5) {
      Pr[i][j][k][UU] = 0.;
      pflag[i][j][k] = 0;
      fixup1zone(i, j, k, Pr[i][j][k]);
    }*/
    ///////////////////////////////////////

    // check for inversion error //
    if (pflag[i][j][k]) {
      fail_save[i][j][k] = 1;
    }

  } // ZLOOP

  // stop timer //
  timer_stop(TIMER_UPDATE);

}

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
#endif // RADIATION //

/*****************************************************************************/
/* calculate the numerical flux */
double fluxcalc(grid_prim_type Pr) {

  // define variables //
  double P_l[NMAX + 2 * NG][NVAR], P_r[NMAX + 2 * NG][NVAR];
  double cij, cmax1, cmax2, cmax3;
  double Ptmp[NMAX + 2 * NG][NVAR];

  // initialize //
  cmax1 = cmax2 = cmax3 = 0.;

#pragma omp parallel private(Ptmp, P_l, P_r, cij) reduction(max:cmax1) reduction(max:cmax2) reduction(max:cmax3)
  {

    /*########################################################################################*/
    // loop over the first dimension //
#pragma omp for collapse(2) nowait
    JSLOOP(-1, N2) {
      KSLOOP(-1, N3) {

        // backup primitive //
        ISLOOP(-NG, N1 - 1 + NG) PLOOP Ptmp[i][ip] = Pr[i][j][k][ip];

        // do reconstruction //
        reconstruct(Ptmp, N1, P_l, P_r);

        // given primitive variables at cell interface, solve the Riemann problem //
        ISLOOP(0, N1) {
          lr_to_flux(P_r[i - 1], P_l[i], &(ggeom[i][j][FACE1]), 1, F1[i][j][k], &cij);
          cmax1 = (cij > cmax1 ? cij : cmax1);
        } // ISLOOP

      }   // KSLOOP
    }     // JSLOOP

    /*########################################################################################*/
    // loop over the second dimension //
#pragma omp for collapse(2) nowait
    ISLOOP(-1, N1) {
      KSLOOP(-1, N3) {

        // backup primitive //
        JSLOOP(-NG, N2 - 1 + NG) PLOOP Ptmp[j][ip] = Pr[i][j][k][ip];
  
        // reconstruct //
        reconstruct(Ptmp, N2, P_l, P_r);

        // given primitive variables at cell interface, solve the Riemann problem //
        JSLOOP(0, N2) {
          lr_to_flux(P_r[j - 1], P_l[j], &(ggeom[i][j][FACE2]), 2, F2[i][j][k], &cij);
          cmax2 = (cij > cmax2 ? cij : cmax2);
        } // JSLOOP

      }   // KSLOOP
    }     // ISLOOP

    /*########################################################################################*/
    // loop over the third dimension //
#pragma omp for collapse(2)
    ISLOOP(-1, N1) {
      JSLOOP(-1, N2) {

        // backup primitive //
        KSLOOP(-NG, N3 - 1 + NG) PLOOP Ptmp[k][ip] = Pr[i][j][k][ip];

        // reconstruct //
        reconstruct(Ptmp, N3, P_l, P_r);

        // given primitive variables at cell interface, solve the Riemann problem //
        KSLOOP(0, N3) {
          lr_to_flux(P_r[k - 1], P_l[k], &(ggeom[i][j][FACE3]), 3, F3[i][j][k], &cij);
          cmax3 = (cij > cmax3 ? cij : cmax3);
        } // KSLOOP

      }   // JSLOOP
    }     // ISLOOP

  }       // omp parallel

  // Otherwise timestep changes with MPI! //
  cmax1 = mpi_max(cmax1);
  cmax2 = mpi_max(cmax2);
  cmax3 = mpi_max(cmax3);

  // time step //
  double ndt1 = cour * dx[1] / cmax1;
  double ndt2 = cour * dx[2] / cmax2;
  double ndt3 = cour * dx[3] / cmax3;

  // return //
  return (1. / (1. / ndt1 + 1. / ndt2 + 1. / ndt3));
}

/*****************************************************************************/
/* flux-ct scheme, see toth et al */
void flux_ct(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3) {

  // edge-centered emf //
  static double emf1[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];
  static double emf2[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];
  static double emf3[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];

  // Calculate EMFs via average to corners (Toth approach) //
#pragma omp parallel
  {

    // get emf from fluxes //
#pragma omp for collapse(3)
    ZSLOOP(0, N1, 0, N2, 0, N3) {
      emf3[i][j][k] = 0.25 *  (F1[i][j][k][B2] + F1[i][j - 1][k][B2] -
                               F2[i][j][k][B1] - F2[i - 1][j][k][B1]);
      emf2[i][j][k] = -0.25 * (F1[i][j][k][B3] + F1[i][j][k - 1][B3] -
                               F3[i][j][k][B1] - F3[i - 1][j][k][B1]);
      emf1[i][j][k] = 0.25 *  (F2[i][j][k][B3] + F2[i][j][k - 1][B3] -
                               F3[i][j][k][B2] - F3[i][j - 1][k][B2]);
    }

    // Rewrite EMFs as fluxes, after Toth //
#pragma omp for collapse(3) nowait
    ZSLOOP(0, N1, 0, N2 - 1, 0, N3 - 1) {
      F1[i][j][k][B1] = 0.;
      F1[i][j][k][B2] = 0.5 * (emf3[i][j][k] + emf3[i][j + 1][k]);
      F1[i][j][k][B3] = -0.5 * (emf2[i][j][k] + emf2[i][j][k + 1]);
    }
#pragma omp for collapse(3) nowait
    ZSLOOP(0, N1 - 1, 0, N2, 0, N3 - 1) {
      F2[i][j][k][B1] = -0.5 * (emf3[i][j][k] + emf3[i + 1][j][k]);
      F2[i][j][k][B2] = 0.;
      F2[i][j][k][B3] = 0.5 * (emf1[i][j][k] + emf1[i][j][k + 1]);
    }
#pragma omp for collapse(3)
    ZSLOOP(0, N1 - 1, 0, N2 - 1, 0, N3) {
      F3[i][j][k][B1] = 0.5 * (emf2[i][j][k] + emf2[i + 1][j][k]);
      F3[i][j][k][B2] = -0.5 * (emf1[i][j][k] + emf1[i][j + 1][k]);
      F3[i][j][k][B3] = 0.;
    }

  } // omp parallel
}

/*****************************************************************************/
/* given primitive variables at cell interfaces, solve riemann problem */
void lr_to_flux(double P_l[NVAR], double P_r[NVAR], struct of_geom *geom, int dir, double Flux[NVAR], double *maxSpeed) {

  // variables //
  struct of_state state_l, state_r;
  double F_l[NVAR], F_r[NVAR], U_l[NVAR], U_r[NVAR];
  double cmax_l, cmax_r, cmin_l, cmin_r, cmax, cmin, ctop;

  // sanity check //
  if (geom->g < SMALL) {
    PLOOP Flux[ip] = 0.;
    *maxSpeed      = 0.;
    return;
  }

  // get state vector //
  get_state(P_l, geom, &state_l);
  get_state(P_r, geom, &state_r);

  // conservative variables //
  primtoflux(P_l, &state_l, dir, 0, geom, F_l);
  primtoflux(P_r, &state_r, dir, 0, geom, F_r);

  // primitive variables //
  primtoflux(P_l, &state_l, 0, 0, geom, U_l);
  primtoflux(P_r, &state_r, 0, 0, geom, U_r);

  // mhd charateristic speed //
  mhd_vchar(P_l, &state_l, geom, dir, &cmax_l, &cmin_l);
  mhd_vchar(P_r, &state_r, geom, dir, &cmax_r, &cmin_r);

  // signal speed //
  cmax = fabs(MY_MAX(MY_MAX(0., cmax_l), cmax_r));
  cmin = fabs(MY_MAX(MY_MAX(0., -cmin_l), -cmin_r));
  ctop = MY_MAX(cmax, cmin);

  // Riemann flux by LF scheme //
  PLOOP Flux[ip] = 0.5 * (F_l[ip] + F_r[ip] - ctop * (U_r[ip] - U_l[ip]));

  // assign maximum signal speed //
  *maxSpeed = ctop;
  
}

/***********************************************************************************
 * 
 * Apply perturbation to the primitive variables, copy from koral 
 * 
 ***********************************************************************************/
void apply_perturbation(grid_prim_type Pf) {
    
  /* loop over domain */
#pragma omp parallel for collapse(3)
  ZLOOP {

    // initialize changes to 0 ;
    double dvx=0.;
    double dvy=0.;
      
    // Generate a random number between 0 and 1
    // Scale and shift to get a number between -1 and 1
    double rand_float = (double)rand() / RAND_MAX;
    dvx = 2 * rand_float - 1;
    Pf[i][j][k][U1] += 1e-5*dvx;

    // Generate a random number between 0 and 1
    // Scale and shift to get a number between -1 and 1
    rand_float = (double)rand() / RAND_MAX;
    dvy = 2 * rand_float - 1;
    Pf[i][j][k][U2] += 1e-5*dvy;
  
  }

  /* fixup, apply floors */
  fixup(Pf, extra);

  /* apply boundary condition */
  bound_prim(Pf);
    
}

