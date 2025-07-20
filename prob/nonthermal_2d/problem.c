/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR ADVECTION OF A PASSIVE SCALAR in 2D                 *
 *                                                                            *
 ******************************************************************************/

// headers //
#include "decs.h"

/*##############################################################################################*/

/*********************************************************************/
/* set problem runtime parameters */
void set_problem_params()
{
  // nothing to set //
}

/*********************************************************************/
/* initialize the problem */
void init_prob()
{

  // sanity check //
  #if METRIC != MINKOWKI
  fprintf(stderr, "ERROR: This problem must be run with the Minkowski metric!\n");
  exit(1);
  #endif // METRIC != MINKOWSKI

  // print out //
  printf("Beginning problem setup.\n");

  // nonthermal and thermal electrons handled in init_electons //
  ZLOOP {
    PLOOP P[i][j][k][ip] = 0.;
    P[i][j][k][RHO] = 1.0;
    P[i][j][k][UU] = (8.6e-4)*(8.6e-4)*P[i][j][k][RHO]/gam/(gam - 1);
  }

  // print out //
  printf("Problem setup complete.\n");
  
}
