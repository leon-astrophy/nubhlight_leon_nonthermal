/******************************************************************************
 *                                                                            *
 * MAIN.C                                                                     *
 *                                                                            *
 * ESTABLISH MPI COMMUNICATION, LOOP OVER TIME, COMPLETE                      *
 *                                                                            *
 ******************************************************************************/
 
/* include headers */
#include "decs.h"
#include "defs.h"
#include <time.h>

/* define function used here */
void init_first();
void init_core();
void init_final();

/*##################################################################################################################*/

/**************************************************************************************/
/* Main chunck of the program */
int main(int argc, char *argv[]) {

  /* nonthermal need to be setup with electrons */
  #if NONTHERMAL && !ELECTRONS
  fprintf(stderr, "Error! NONTHERMAL need to setup with ELECTRONS!\n");
  exit(1);
  #endif

  /* set OPENMP threads to 1 */
  #if !OPENMP
  omp_set_num_threads(1);
  #endif

  // Check for minimal required MPI thread support //
  int threadSafety;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &threadSafety);
  if (threadSafety < MPI_THREAD_FUNNELED) {
    fprintf(stderr, "Thread support < MPI_THREAD_FUNNELED. Unsafe.\n");
    exit(1);
  }

  /* initialize MPI stuff */
  init_mpi();

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

  /* print out for debug */
  if (mpi_myrank() == 0) {
    ////////////////////
    // clang-format off
    ////////////////////
    fprintf(stdout, "\n          ************************************************************\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *                          BHLIGHT                         *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *          Ryan, Dolence & Gammie ApJ 807:31, 2015         *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *  B R Ryan                                                *\n");
    fprintf(stdout, "          *  J C Dolence                                             *\n");
    fprintf(stdout, "          *  C F Gammie                                              *\n");
    fprintf(stdout, "          *  S M Ressler                                             *\n");
    fprintf(stdout, "          *  J M Miller                                              *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *                         RESOLUTION                       *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *  N1TOT : %04d (N1CPU: %04d)                              *\n", N1TOT, N1CPU);
    fprintf(stdout, "          *  N2TOT : %04d (N2CPU: %04d)                              *\n", N2TOT, N2CPU);
    fprintf(stdout, "          *  N3TOT : %04d (N3CPU: %04d)                              *\n", N3TOT, N3CPU);
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *                          OPTIONS                         *\n");
    fprintf(stdout, "          *                                                          *\n");
    #if OPENMP
    fprintf(stdout, "          *  OPENMP                                                  *\n");
    #endif
    #if EOS == EOS_TYPE_TABLE
    fprintf(stdout, "          *  TABULATED EOS  (O'connor and Ott. CQG 27:114103, 2010)  *\n");
    #endif
    #if EOS == EOS_TYPE_GAMMA
    fprintf(stdout, "          *  IDEAL GAS EOS                                           *\n");
    #endif
    #if RADIATION == RADTYPE_LIGHT
    fprintf(stdout, "          *  RADIATION: PHOTON TRANSPORT                             *\n");
    #elif RADIATION == RADTYPE_NEUTRINOS
    fprintf(stdout, "          *  RADIATION: NEUTRINO TRANSPORT                           *\n");
    #endif
    #if TRACERS
    fprintf(stdout, "          *  TRACERS                                                 *\n");
    #endif
    #if ELECTRONS
    fprintf(stdout, "          *  ELECTRONS      (Ressler et al. MNRAS 454:1848, 2015)    *\n");
    #endif
    #if EXIT_ON_INIT
    fprintf(stdout, "          *  EXIT ON INIT                                            *\n");
    #endif
    #if NVAR_PASSIVE > 0
    fprintf(stdout, "          *  PASSIVE VARIABLES : %02d                                  *\n", NVAR_PASSIVE);
    #endif
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *                          SYNTAX                          *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *  -p /path/to/param_file                                  *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          ************************************************************\n\n");
    /////////////////////
    // clang-format on
    /////////////////////
  }

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

  // Read command line arguments //
  strcpy(outputdir, ""); // Default value
  char pfname[STRLEN];
  for (int n = 0; n < argc; n++) {
    // Check for argv[n] of the form '-*' //
    if (*argv[n] == '-' && *(argv[n] + 1) != '\0' && *(argv[n] + 2) == '\0' &&
        n < argc - 1) {
      if (*(argv[n] + 1) == 'p') { // Set parameter file path
        strcpy(pfname, argv[++n]);
      }
    }
  }

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  
  /* Set number of threads if OPENMP is turn on */
  #pragma omp parallel
  {
  #pragma omp master
    { nthreads = omp_get_num_threads(); } // omp master
  }                                       // omp parallel

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

  // Perform initializations, either directly, via checkpoint, or initialize //
  // with GRMHD data //
  init_params(pfname);
  init_first();

  /* initialize nonthermal electrons gammas */
  #if NONTHERMAL
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  set_units();
  set_nonthermal_gammas();
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  #endif

  // initialize input-output //
  init_io();

  /* check if we are doing restart */
  is_restart = restart_init();

  /* if not, initialize the data */
  if (!is_restart) {

    // call the initialization subroutine 
    init_core();

    // Set primitive variables //
    if (mpi_io_proc()) {
      printf("Init from GRMHD? %i\n", strcmp(init_from_grmhd, "No") != 0);
    }

    // Initialize from problem.c if we are not restarting from a GRMHD snapshot //
    if (strcmp(init_from_grmhd, "No") == 0) { 
      init_prob();
      #if ELECTRONS
      init_electrons();
      #endif
    // Otherwise, initialize fluid from GRMHD restart file //
    } else { 
      init_fluid_restart();
    }

    // initialize every thing //
    init_final();

    // print to debug //
    if (mpi_myrank() == 0)
      fprintf(stdout, "Initial conditions generated\n\n");
  }

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

  // Initial diagnostics //
  diag(DIAG_INIT);

  // Print to debug //
  if (mpi_io_proc())
    fprintf(stdout, "t = %e tf = %e\n", t, tf);

  // Output dump files //
  if (!is_restart)
    diag(DIAG_DUMP);

  // Set up timers //
  time_init();

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

  // Print out for debugging //
  if (EXIT_ON_INIT) { // TODO: make this a runtime parameter
    if (mpi_io_proc()) {
      fprintf(stdout, "Finished!\n");
    }
    MPI_Finalize();
    return 0;
  }

  // Print out for debugging //
  if (mpi_io_proc())
    fprintf(stdout, "\nEntering main loop\n");

  // Set flag //
  int dumpThisStep = 0;

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

  // Loop over time //
  while (t < tf) {

    // set flag //
    dumpThisStep = 0;

    // set timer //
    timer_start(TIMER_ALL);

    // Step variables forward in time //
    step();

    // increase step size // 
    nstep++;

    // Per-step output //
    #if RADIATION
    /*##################################################################*/
    step_made_all  = mpi_reduce_int(step_made);
    step_abs_all   = mpi_reduce_int(step_abs);
    step_scatt_all = mpi_reduce_int(step_scatt);
    step_lost_all  = mpi_reduce_int(step_lost);
    step_rec_all   = mpi_reduce_int(step_rec);
    step_sent_all  = mpi_reduce_int(step_sent);
    step_rcvd_all  = mpi_reduce_int(step_rcvd);
    step_fail_all  = mpi_reduce_int(step_fail);
    step_tot_all   = mpi_reduce_int(step_tot);
    tracer_tot_all = mpi_reduce_int(tracer_tot);
    step_tot_max   = mpi_max(step_tot);
    step_tot_min   = mpi_min(step_tot);
    // Check load balancing //
    // 0 -> well balanced   //
    // 1 -> poorly balanced //
    load_imbalance = 1.0 - ((float)step_tot_min / (float)step_tot_max);
    /*##################################################################*/
    #endif // RADIATION

    // MPI synchronization //
    mpi_sync_output();

    /////////////////////////////////////////////////////////////////////////////
    /*  fprintf(stdout, "[%i] made = %d abs = %d scatt = %d lost = %d rec = %d
    sent = %d rcvd = %d fail = %d tot = %d\n", mpi_myrank(), step_made,
    step_abs, step_scatt, step_lost, step_rec, step_sent, step_rcvd, step_fail,
    step_tot); mpi_sync_output();*/
    /////////////////////////////////////////////////////////////////////////////

    // print out for debug //
    if (mpi_io_proc()) {
      fprintf(stdout, "t = %10.5g dt = %10.5g n = %8d\n", t, dt, nstep);
      #if RADIATION
      /*##################################################################*/
      fprintf(stdout,
          "[%i] made = %d abs = %d scatt = %d lost = %d rec = %d sent = %d "
          "rcvd = %d fail = %d tot = %d target = %g\n",
          mpi_myrank(), step_made, step_abs, step_scatt, step_lost, step_rec,
          step_sent, step_rcvd, step_fail, step_tot, nph_per_proc);
      fprintf(stdout,
          "ALL made = %d abs = %d scatt = %d lost = %d rec = %d sent = %d rcvd "
          "= %d fail = %d tot = %d target = %g\n",
          step_made_all, step_abs_all, step_scatt_all, step_lost_all,
          step_rec_all, step_sent_all, step_rcvd_all, step_fail_all,
          step_tot_all, nph_per_proc * mpi_nprocs());
      #if TRACERS
      /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
      fprintf(stdout, "ALL tracers:     %i\n", tracer_tot_all);
      fprintf(stdout, "ALL non-tracers: %i\n", step_tot_all - tracer_tot_all);
      /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
      #endif // TRACERS
      if (step_sent_all != step_rcvd_all) {
        printf("PHOTON MPI ANOMALY!\n");
        printf("sent = %i rcvd = %i\n", step_sent_all, step_rcvd_all);
      }
      if (step_tot > 1.e8) {
        printf("Too many superphotons on this node!\n");
        exit(-1);
      }
      /*##################################################################*/
      #endif // RADIATION
    }

    #if RADIATION
    /*###################################################*/
    update_superphoton_resolution(P, extra);
    #if RADIATION == RADTYPE_NEUTRINOS
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    lepton_lost_step = mpi_reduce(lepton_lost_local);
    count_leptons(P, dt, nstep);
    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
    #endif // RADIATION == RADTYPE_NEUTRINOS
    /*###################################################*/
    #endif // RADIATION

    // File I/O with set frequencies //
    if (t < tf) {
      if (t >= tdump) {
        dumpThisStep = 1;
        diag(DIAG_DUMP);
        #if LOGDUMPING
        DTd = DTd*10;
        tdump = DTd;
        #else // LOGDUMPING
        tdump += DTd;
        #endif // LOGDUMPING
      }
      if (t >= tlog) {
        diag(DIAG_LOG);
        tlog += DTl;
      }
      if (nstep % DNr == 0) {
        restart_write(RESTART_TEMP);
      }
      if (t >= trestart) {
        restart_write(RESTART_PERM);
        trestart += DTr;
      }
    }

    // Reset all log variables to zero //
    reset_log_variables();

    // stop the timer //
    timer_stop(TIMER_ALL);

    // reset timers at the first step //
    if (nstep == 1)
      timers_reset();

    // Do some report //
    if (nstep % DTp == 0) {
      report_performance();
      #if EOS == EOS_TYPE_TABLE
      /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
      print_root_fcounts();
      /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
      #endif // EOS == EOS_TYPE_TABLE
      #if RADIATION
      /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
      report_load_imbalance();
      #if RADIATION == RADTYPE_NEUTRINOS
      /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
      print_rad_types();
      /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
      #endif // NEUTRINOS
      /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
      #endif // RADIATION
    }

  }

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

  // Do the final diagonsis //
  if (dumpThisStep == 0)
    diag(DIAG_FINAL);

  // Finalize MPI //
  MPI_Finalize();

  // return //
  return 0;
}

/**************************************************************************************/
/* Do some initialization */
void init_first() {

  // Initialize some global variables //
  reset_log_variables();

  // Set step to zero //
  nstep = 0;

  // Set random numbers // 
  init_random(300 * mpi_myrank());

  // Initialize EOS //
  init_EOS();

  // EOS Table //
  #if EOS == EOS_TYPE_TABLE
  /*#########################*/
  initialize_root_fcounts();
  /*#########################*/
  #endif

  // Passive scalar initialize //
  #if NVAR_PASSIVE > 0
  /*@@@@@@@@@@@@@@@@@@@@@@@@@*/
  init_passives();
  name_passives();
  /*@@@@@@@@@@@@@@@@@@@@@@@@@*/
  #endif

}

/**************************************************************************************/
/* Do some initialization */
void init_core() {

  // zero the arrays //
  zero_arrays();

  // reset the dump variables //
  reset_dump_variables();

  // Set the grid and geometry //
  set_grid();

  // initialize some quantities //
  failed    = 0;
  t         = 0.;
  dump_cnt  = 0;
  rdump_cnt = 0;

  // Set units //
  #if NEED_UNITS
  /*@@@@@@@@@@@@@@@@@@@@*/
  set_units();
  /*@@@@@@@@@@@@@@@@@@@@*/
  #endif // NEED_UNITS

  // Eos Table // 
  #if EOS == EOS_TYPE_TABLE
  /*@@@@@@@@@@@@@@@@@@@@@@@@*/
  EOS_SC_print_table_mins();
  /*@@@@@@@@@@@@@@@@@@@@@@@@*/
  #endif

  // Electron fraction as passive scalar //
  #if EOS == EOS_TYPE_TABLE
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  PASSTYPE(YE)    = PASSTYPE_NUMBER;
  PASSTYPE(YE_EM) = PASSTYPE_NUMBER;
  #if METRIC == MKS
  /*#################################*/
  PASSTYPE(ATM) = PASSTYPE_NUMBER;
  /*#################################*/
  #endif // METRIC
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  #endif // EOS_TYPE_TABLE

  // Initialize radiation //
  #if RADIATION
  /*@@@@@@@@@@@@@@@@@@@@@@@@*/
  init_rad(P);
  dt = cour * dt_light_min;
  /*@@@@@@@@@@@@@@@@@@@@@@@@*/
  #endif

}

/**************************************************************************************/
/* Do some initialization */
void init_final() {

  // Set superphoton weight //
  #if RADIATION
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  set_weight(P, extra);
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  #endif

  // ceiling and floors //
  fixup(P, extra);

  // boundary condition //
  bound_prim(P);

  // set output frequencies //
  tdump    = DTd;
  trestart = DTr;
  tlog     = DTl;
  dtsave   = dt;

  // dump the grid //
  dump_grid();
  
}
