/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR NONTHERMAL ELECTRONS PROBLEM                        *
 *                                                                            *
 ******************************************************************************/

// headers //
#include "decs.h"

/*#######################################################################################*/

/********************************************************************/
/* read the problem parameter */
void set_problem_params()
{
  // no parameters //
}

/********************************************************************/
/* assign initial condition */
void init_prob()
{

  // sanity check //
  #if METRIC != MINKOWKI
  fprintf(stderr, "ERROR: This problem must be run with the Minkowski metric!\n");
  exit(1);
  #endif // METRIC != MINKOWSKI

  // print out //
  printf("Beginning problem setup.\n");

  // Loop over //
  // Leon: initial condition assign once only, no OMP //
  ZLOOP {
    PLOOP P[i][j][k][ip] = 0.;
  }

  // print out //
  printf("Problem setup complete.\n");

}
