/******************************************************************************
 *                                                                            *
 * PASSIVE.C                                                                  *
 *                                                                            *
 * PASSIVE SCALARS                                                            *
 *                                                                            *
 ******************************************************************************/

/* include headers */
#include "decs.h"

/*#################################################################################*/

// compile only if we decided to advect passive scalar //
#if NVAR_PASSIVE > 0
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

/**********************************************************************/
/* fixup passive variables */
void fixup_passive(int i, int j, int k, double pv[NVAR], double pv_prefloor[NVAR]) {
  PASSLOOP { // advected numbers should not be changed by fixup //
    if (PASSTYPE(ipass) == PASSTYPE_NUMBER) {
      pv[ipass] = pv_prefloor[ipass];
    }
    if (do_passive_fixup[PASSELEM(ipass)] != NULL) {
      do_passive_fixup[PASSELEM(ipass)](i, j, k, pv, pv_prefloor);
    }
  }
}

/**********************************************************************/
/* initialize pointers to the fixup subroutine */
void __attribute__((weak)) init_passives() {
  PASSLOOP { do_passive_fixup[PASSELEM(ipass)] = NULL; }

  #if EOS == EOS_TYPE_TABLE
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  do_passive_fixup[PASSELEM(YE)]    = do_ye_fixup;
  do_passive_fixup[PASSELEM(YE_EM)] = do_ye_em_fixup;
  #if METRIC == MKS
  do_passive_fixup[PASSELEM(ATM)] = do_atm_fixup;
  #endif // METRIC == MKS
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  #endif // EOS == EOS_TYPE_TABLE
}

/**********************************************************************/
/* initialize pointers to the passive scalar name */
void __attribute__((weak)) name_passives() {
  char name[STRLEN];
  PASSLOOP {
    sprintf(name, "var%d", PASSELEM(ipass));
    strcpy(PASSNAME(ipass), name);
  }

  #if EOS == EOS_TYPE_TABLE
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  // The first passive scalar is always Ye, the electron/proton fraction
  strcpy(PASSNAME(YE), "Ye");
  strcpy(PASSNAME(YE_EM), "Ye_em");
  #if METRIC == MKS
  strcpy(PASSNAME(ATM), "ATM");
  #endif // METRIC
  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  #endif // EOS_TYPE_TABLE
}

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
#endif // NVAR_PASSIVE > 0
