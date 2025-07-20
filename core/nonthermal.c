/******************************************************************************
 *                                                                            *
 * NONTHERMAL.C                                                               *
 *                                                                            *
 * Quantities and functions for nonthermal electron evolution                 *
 *                                                                            *
 * Leon: should make reference to KORAL and Andrew's paper I guess            *
 * Leon: all ME here should be converted to code unit !!!!!!!!!!!!            *
 * Leon: migrate from JAD's ebhlight                                          *
 *                                                                            *
 ******************************************************************************/

/* header */
#include "decs.h"

#if NONTHERMAL
/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

/************************************************************************************
 * 
 * @brief Populate log10nteGammas and nteGammas with appropriate values 
 * based on NTGAMMAMAX and NTGAMMAMIN
 * 
 ************************************************************************************/
void set_nonthermal_gammas()
{

  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  // TODO: It's not so obvious what I should do about bins. 
  // i.e, should the values be centered in the bins and I set seperate edge values or something?
  // TODO: May want to use ln instead of log10
  // TODO: when I add injection, there may be some check here to make sure injection min/max is valid
  // Leon: the conversion from ln to log10 space is trivial
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /* gamma bin size in the log10 scale */
  log10BinSpace = (log10(NTGAMMAMAX) - log10(NTGAMMAMIN))/((double)(NTEBINS-1));

  /* Leon's patch, natural log bin space */
  logeBinSpace = log10BinSpace*log(10);

  /* lower bounds */
  log10nteGammas[0] = log10(NTGAMMAMIN);
  nteGammas[0] = NTGAMMAMIN;

  /* for integrating the distribution */
  double normterms[NTEBINS];
  double currgamma;

  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  /* Assign gamma. Note: bin number runs from 0 to NTEBINS - 1 */
  for (int ig = 1; ig<NTEBINS; ig++) {
    log10nteGammas[ig] = log10nteGammas[ig-1] + log10BinSpace;
    nteGammas[ig] = pow(10.0, log10nteGammas[ig]);
  }

  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  /* now assign normalization term */
  NTEGAMMALOOP {
    currgamma = nteGammas[ig];
    if((currgamma >= gammainjmin) && (currgamma <= gammainjmax)){ 
      normterms[ig] = (currgamma-1)*pow(currgamma, -PLAW); 
    } else {
      normterms[ig] = 0.0;
    }
  }

  /* this is just m_e*int(gam-1)*gam^-p */
  normterm = ME_CODE*pow(CL_CODE, 2.0)*gamma_integral(normterms);

  /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  /* Loop over */
  NTEGAMMALOOP {
    if((nteGammas[ig] >= gammainjmin) && (nteGammas[ig] <= gammainjmax)){ 
      normterms[ig] = pow(nteGammas[ig], -PLAW);
    } else {
      normterms[ig] = 0.0;
    }
  }

  /* integrating the powerlaw distribution */
  powlaw_int = gamma_integral(normterms);

}

/*********************************************************************************
 * 
 * @brief Computes an integral over gamma for a matrix of size NTEBINS
 * 
 * @param ureal Matrix to integrate
 * @return double Integral result
 * 
 *********************************************************************************/
double gamma_integral(double *ureal){

  //////////////////////////////////////////////////////////////////////
  // TODO: Would logspace or a simpson's integration be better here?
  // TODO: Is the edge handled properly?
  //////////////////////////////////////////////////////////////////////

  // variables 
  double utot = 0;

  //////////////////////////////////////////////////////////////////////
  //double dg;
  //NTEGAMMALOOP{
  //  if(ig == NTEBINS-1)
  //    dg = pow(10,log10nteGammas[ig]+log10BinSpace) - nteGammas[ig];
  //  else 
  //    dg = nteGammas[ig+1] - nteGammas[ig];
  //    utot += ureal[ig]*dg;
  //}
  //////////////////////////////////////////////////////////////////////

  /* Leon's patch: natural log space integration */
  NTEGAMMALOOP {
    utot += ureal[ig]*nteGammas[ig]*logeBinSpace;
  }

  // return //
  return utot;
    
}

/******************************************************************************/
/* I guess this is adiabatic gas compression and expansion */
/* Pi is the primitive before explicit step */
/* Pf is the primitive after explicit step */
/* Pf is the variable we wish to update */
void nonthermal_adiab(grid_prim_type Pi, grid_prim_type Ps, grid_prim_type Pf, double Dt){
    
  //////////////////////////////////////////////////////////////////////////////////
  // I worry since the expasion depends on the primitive left and right values 
  // that the values may be overwritten as the expansion is calculated...
  // #ifdef ART_ADIAB
  // #define SEMIART_ADIAB (ART_ADIAB)
  // #endif
  // Leon: no need to worry, just simple central differences
  // But might fails close to strong shocks!
  //////////////////////////////////////////////////////////////////////////////////

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* Leon: not sure what SEMIART_ADIAB means here */
  /* wait, is it artificial adiabatic rate?       */
  /* Oh I get it! Chael section 4.2               */
  #ifdef SEMIART_ADIAB
  double X[NDIM];
  // Leon: added OMP here //
  #pragma omp parallel for collapse(3) // schedule(dynamic)
  ZLOOPALL {
    coord(i, j, k, CENT, X);
    Pf[i][j][k][U1] = X[1]*SEMIART_ADIAB/3;
    Pf[i][j][k][U2] = X[2]*SEMIART_ADIAB/3;
    Pf[i][j][k][U3] = X[3]*SEMIART_ADIAB/3;
  }
  #endif 

  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  /* apply nonthermal heating/cooling */
  #pragma omp parallel for collapse(3) // schedule(dynamic)
  ZLOOP {
    nonthermal_adiab_zone(i, j, k, Pi, Ps, Pf, Dt);
  }
}

/********************************************************************************************
 * 
 * @brief Applies the cooling/heating caused by adiabatic expansion/compression in one zone
 * 
 * @param i Dimension 1 index
 * @param j Dimension 2 index
 * @param k Dimension 3 index
 * @param Pi, Ps, Pf Full primitives matrix (need neighboring zones as well)
 * 
 * Pi is the primitive before explicit step 
 * Pf is the primitive after explicit step 
 * Pf is the variable we wish to update 
 * 
 ********************************************************************************************/
void nonthermal_adiab_zone(int i, int j, int k, grid_prim_type Pi, grid_prim_type Ps, grid_prim_type Pf, double Dt){

  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  // TODO: Definitely needs testing...
  // TODO: Include the electrons returning to the thermal distribution (Chael section 3.2 ii)
  // TODO: Viscous dissipation rate and injection terms into thermal and nonthermal pops (Chael section 3.2 iii and eq. 26/29)
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /*#############################################################################################*/
  /* check if we defined scaling option */
  #ifndef ADIABTIC_SCALING
  #define ADIABTIC_SCALING (0)
  #endif

  /* adiabatic rate, depends on both Ps and Pf */
  double adiab = calc_expansion(i, j, k, Pi, Ps, Pf, Dt);
    
  /* artifical adiabatic rate */
  #ifdef ART_ADIAB
  adiab = ART_ADIAB;
  #endif

  /*#############################################################################################*/
  /* variables */
  double nprime[NTEBINS], deltan[NTEBINS], ngamma[NTEBINS];

  /* loop over gamma, assign the number density distribution */
  /* ngamma is previous step to get the source term properly */
  NTEGAMMALOOP ngamma[ig] = Ps[i][j][k][ig+NTESTART];

  /////////////////////////////////////////////////////////////////////////
  // Find the initial values of n and u to compare to the final values
  // double n_tot_start = gamma_integral(ngamma);
  /////////////////////////////////////////////////////////////////////////

  /*#############################################################################################*/
  /* Leon's patch, upwind in log space */
  nonthermal_adiab_lf(adiab, ngamma, nprime);
  
  /*#############################################################################################*/
  // state and geometric variables //
  struct of_state q;
  struct of_geom *geom;

  // assign geometric variables //
  geom = &ggeom[i][j][CENT];

  // 4-vectors //
  get_state(Ps[i][j][k], geom, &q);

  // Find dtau, convert from coordinate to proper time //
  // Leon: are we sure the proper time is scaled by ucon //
  // computed by primitive after explicit step ? //
  double dtau = Dt/(q.ucon[0]);

  /*#############################################################################################*/

  #if ADIABTIC_SCALING
  // Find change in n and the resulting real energy change (Chael eq. 47 and 11)
  double uexpected[NTEBINS], utotexpected;
  double ureal[NTEBINS], utot;
  #endif

  /* gamma dot */
  double gdot; 

  /*#############################################################################################*/
  /* loop over */
  NTEGAMMALOOP {
      
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Version with semi-analytic derivative
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // deltan[ig] = dtau*(adiab/3.)*( (1+pow(nteGammas[ig],-2))*ngamma[ig] + (nteGammas[ig]-(1/nteGammas[ig]))*nprime[ig] );
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /* find change in the number density distriubtion */
    deltan[ig] = dtau*(adiab/3.)*nprime[ig];
    
    /* Leon: Should there be a dtau here? the dimension looks correct but Chael 2017 doesn't have */
    /* Leon: I remove the mec^2 here beacuse it will cancelled out during divison */
    #if ADIABTIC_SCALING
    gdot = (-1.0/3.0)*(adiab)*(nteGammas[ig] - 1./nteGammas[ig]);
    uexpected[ig] = dtau*gdot*ngamma[ig];
    ureal[ig] = (nteGammas[ig]-1.)*deltan[ig];
    #endif
      
  }

  /*#############################################################################################*/

  #if ADIABTIC_SCALING
  // Rescale change by expected vs actual energies to preserve energy
  utot = ME_CODE*pow(CL_CODE, 2.0)*gamma_integral(ureal);
  utotexpected = ME_CODE*pow(CL_CODE, 2.0)*gamma_integral(uexpected);
  #endif

  /*#############################################################################################*/
    
  /* right and side */
  double rhs;

  /* loop over */
  NTEGAMMALOOP {

    /* Leon's comment: I also doubt if this is needed in ebhlight */
    /* The scaling is to ensure the energy is computed consistent with Eq.49 in Chael 2017 */
    /* So that the estimation of viscous heating is correct. But in ebhlight, the */
    /* viscous heating is computed directly by the total entropy equation */
    #if ADIABTIC_SCALING
    if ((fabs(utotexpected) > SMALL) && (fabs(utot) > SMALL)) deltan[ig] *= utotexpected/utot;
    #endif

    // Apply the change to the bins
    // note that deltan[ig] is calculate based on 
    // variables at the previous substep (Ps) 
    // limite the change to 20% of the primitive 
    rhs = deltan[ig];
    #if LIMIT_NTE
    /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
    if(Pf[ig + NTESTART] > 0) {
      rhs = MY_SIGN(rhs)*MY_MIN(fabs(rhs), fabs(0.2*Pf[i][j][k][ig+NTESTART]));
    }
    /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
    #endif
    Pf[i][j][k][ig+NTESTART] += rhs;

    // Floor bins (ie. no negative n)
    if(Pf[i][j][k][ig+NTESTART] < 0) Pf[i][j][k][ig+NTESTART] = 0.;

  }

  /*#############################################################################################*/
  /* Leon: I am quite unsure about the implementation in the following section */
  /* Catch energy coming out of the nonthermal electrons */
  double ncool, qcool;

  /////////////////////////////////////////////////////////////////////////
  // Expanding
  // Leon: I think this is incorrect. The change in density in a cell 
  // can be due to advection to lower/higher gamma, not necessary escape
  // out of the distribution 
  /////////////////////////////////////////////////////////////////////////
  //if(adiab > 0){
  //  ncool = -deltan[0];
  //  qcool = ncool*ME_CODE*(nteGammas[0]-1);
  //}
  // Compressing, small amount can escape out the top
  //else{
  //  ncool = -deltan[NTEBINS-1];
  //  qcool = ncool*ME_CODE*(nteGammas[NTEBINS-1]-1);
  //}
  /////////////////////////////////////////////////////////////////////////

  /* Leon's patch, account for particle and energy flux connecting to the thermal path */
  /* Note: In Koral, ndot is computed as per the TOTAL change of thermal number density */
  /* Whereas here we calculate by the fluxes through the gamma boundary */
  if(adiab > 0){
    gdot = (-1.0/3.0)*(adiab)*(nteGammas[0] - 1./nteGammas[0]);
    ncool = gdot*ngamma[0];
    qcool = ME_CODE*pow(CL_CODE,2)*(nteGammas[0]-1)*ncool;
  } else {
    gdot = (-1.0/3.0)*(adiab)*(nteGammas[NTEBINS-1] - 1./nteGammas[NTEBINS-1]);
    ncool = gdot*ngamma[NTEBINS-1];
    qcool = ME_CODE*pow(CL_CODE,2)*(nteGammas[NTEBINS-1]-1)*ncool;
  }
  
  /* take absolute value because particle joining thermal distribution always carries energy */
  ncool = fabs(ncool);
  qcool = fabs(qcool);

  /* heating thermal electrons */
  /* Leon: make sure the heat applied is in code unit */
  /* Also these are source terms computed with the previous substep */
  /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
  //double chempot = calc_potential(Ps[i][j][k]);
  /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
    
  /* we are updating Pf here */
  /* Leon: need discussion here */
  /* Should we also modify the total entropy? */ 
  //apply_thermal_heating(Pf[i][j][k], qcool, ncool, dtau);
  
  /*#############################################################################################*/
    
}

/************************************************************************************************************************
 * 
 * @brief Calculate the divergence of the four-velocity (expansion/contraction term). 
 * Uses {u^\{alpha}}_{;\alpha} = (sqrt(-g)*u^\alpha)_,\alpha / sqrt(-g)
 * 
 * @param i grid coordinate 1
 * @param j grid coordinate 2
 * @param k grid coordinate 3
 * @return double {u^\{alpha}}_{;\alpha}
 *
 * Ps is the primitive before explicit step
 * Pf is the primitive after explicit step
 * Leon: make reference to Koral to set how to calculate the spatial divergence
 * 
 ************************************************************************************************************************/
double calc_expansion(int i, int j, int k, grid_prim_type Pi, grid_prim_type Ps, grid_prim_type Pf, double Dt){

  // variables //
  struct of_state ql, qc, qr;
  struct of_geom *geoml, *geomc, *geomr;
  double Xl[NDIM], Xc[NDIM], Xr[NDIM];
  double du[NDIM];
  double result = 0;

  // Center: //
  geomc = &ggeom[i][j][CENT];
  get_state(Pf[i][j][k], geomc, &qc);
  coord(i,j,k,CENT,Xc);

  // Dimension 0: //
  double pucon[NDIM];
  ucon_calc(Ps[i][j][k], geomc, pucon);
  du[0] = ( (geomc->g)*(qc.ucon[0])-(geomc->g)*(pucon[0]) )/Dt;

  // Dimension 1: //
  geoml = &ggeom[i-1][j][CENT];
  get_state(Ps[i-1][j][k], geoml, &ql);
  coord(i-1,j,k,CENT,Xl);
  geomr = &ggeom[i+1][j][CENT];
  get_state(Ps[i+1][j][k], geomr, &qr);
  coord(i+1,j,k,CENT,Xr);

  // Could add the option later but I'll just do the center derivative for now //
  ///////////////////////////////////////////////////////////////////////////////
  du[1] = ( (geomr->g)*(qr.ucon[1])-(geoml->g)*(ql.ucon[1]) )/(Xr[1]-Xl[1]);
  ///////////////////////////////////////////////////////////////////////////////
  //du[1] = ( (geomr->g)*(qr.ucon[1])-(geomc->g)*(qc.ucon[1]) )/(Xr[1]-Xc[1]);

  // Dimension 2: //
  geoml = &ggeom[i][j-1][CENT];
  get_state(Ps[i][j-1][k], geoml, &ql);
  coord(i,j-1,k,CENT,Xl);
  geomr = &ggeom[i][j+1][CENT];
  get_state(Ps[i][j+1][k], geomr, &qr);
  coord(i,j+1,k,CENT,Xr);

  // Could add the option later but I'll just do the center derivative for now //
  ///////////////////////////////////////////////////////////////////////////////
  du[2] = ( (geomr->g)*(qr.ucon[2])-(geoml->g)*(ql.ucon[2]) )/(Xr[2]-Xl[2]);
  ///////////////////////////////////////////////////////////////////////////////
  //du[2] = ( (geomr->g)*(qr.ucon[2])-(geomc->g)*(qc.ucon[2]) )/(Xr[2]-Xc[2]);

  // Dimension 3: //
  geoml = &ggeom[i][j][CENT];
  get_state(Ps[i][j][k-1], geoml, &ql);
  coord(i,j,k-1,CENT,Xl);
  geomr = &ggeom[i][j][CENT];
  get_state(Ps[i][j][k+1], geomr, &qr);
  coord(i,j,k+1,CENT,Xr);

  // Could add the option later but I'll just do the center derivative for now //
  ///////////////////////////////////////////////////////////////////////////////
  du[3] = ( (geomr->g)*(qr.ucon[3])-(geoml->g)*(ql.ucon[3]) )/(Xr[3]-Xl[3]);
  ///////////////////////////////////////////////////////////////////////////////
  //du[3] = ( (geomr->g)*(qr.ucon[3])-(geomc->g)*(qc.ucon[3]) )/(Xr[3]-Xc[3]);

  // Sum and divide: //
  DLOOP1 result += du[mu];

  // return //
  return result/(geomc->g);
    
}

/************************************************************************************************
 * 
 * Leon's patch: evolving the hyperbolic conservation law in the natural log space
 * for the nonthermal electron population (Chael 2017 Eq 26)
 * Here, I use the simple Lax–Friedrichs method
 * This is to get d/dgamma (n(gamma)(gamma - 1/gamma)), but in the log space
 * 
 * Input - dtau (timestep), adiab (4-velocity divergence), ngamma (number density distriubtion)
 * Output - deltan (the change in the number density distribution)
 * 
 * Note that I solve the hyperbolic system in the natural log space
 * Leon's update: OK seems like upwind scheme is better ... 
 * 
 ************************************************************************************************/
void nonthermal_adiab_lf(double adiab, double *ngamma, double *nprime) {

  /* sign convention */
  double sign;

  /* gamma and left and right handside */
  double gamma_c, gamma_up;

  /* number density at left and right hand side */
  double ngamma_c, ngamma_up;

  /* flux and conservative variables */
  double fluxes_c, fluxes_up;

  /* loop over bins */
  NTEGAMMALOOP {

    /* assign current cell */
    ngamma_c = ngamma[ig];
    gamma_c = nteGammas[ig];
    fluxes_c = ngamma_c*(gamma_c - 1.0/gamma_c);

    /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
    // div V > 0, - (div V) < 0, -ve speed //
    if(adiab > 0){
      sign = 1.0;
      if(ig == (NTEBINS-1)) {
        fluxes_up = 0.0;
      } else{
        ngamma_up = ngamma[ig+1];
        gamma_up = nteGammas[ig+1];
        fluxes_up = ngamma_up*(gamma_up - 1.0/gamma_up);
      } 
    }
    /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
    // div V < 0, - (div V) > 0, +ve speed //
    else{
      sign = -1.0;
      if(ig == 0){
        fluxes_up = 0.0;
      }
      else{
        ngamma_up = ngamma[ig-1];
        gamma_up = nteGammas[ig-1];
        fluxes_up = ngamma_up*(gamma_up - 1.0/gamma_up);
      } 
    }
    /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
    /* calculate change in gamma*ngamma   */
    nprime[ig] = sign*(fluxes_up - fluxes_c)/logeBinSpace;
    nprime[ig] = nprime[ig]/gamma_c;
    
  }
}

/******************************************************************************
 * 
 * Constant injection model 
 * Pf here is the primitive variable we wish to upate 
 * 
 ******************************************************************************/
void const_inject(grid_prim_type Pi, grid_prim_type Ps, grid_prim_type Pf, double Dt) {
    
  /* variables */
  double norms = const_inj*normterm/powlaw_int;
  struct of_state q;
  struct of_geom *geom;
    
  /* Leon: shouldn't Dt here be Dtau? */
  // Leon: added OMP here //
  #pragma omp parallel for collapse(3) // schedule(dynamic)
  ZSLOOP(-NG, N1-1+NG, -NG, N2-1+NG, -NG, N3-1+NG) {

    /* get proper time */
    geom = &ggeom[i][j][CENT];
    get_state(Ps[i][j][k], geom, &q);
    double dtau = Dt/(q.ucon[0]);
    inject_nonthermal(Pf[i][j][k], norms, dtau);
      
  }
    
}

/**********************************************************************************************************************/
/* viscously heat electrons, nonthermal version */
void heat_electrons_zone_nonthermal(int i, int j, int k, double Pi[NVAR], double Ps[NVAR], double Pf[NVAR], double Dt){
  
  // variables //
  double ktotharm, ktotadv, fel, felth, felnth;
    
  // fraction goes to nonthermal particles //
  felnth = get_felnth(i,j,k,Ps);

  // geometry and state //
  struct of_geom *geom = &ggeom[i][j][CENT];
  struct of_state qf;

  // ucov ucon bcov bcon //
  get_state(Pf, geom, &qf);

  // Calculated and advected entropy at final time //
  ktotharm = (gam-1.)*Pf[UU]/pow(Pf[RHO],gam);
  ktotadv = Pf[KTOT];

  // Electron heating fraction // 
  #ifdef FEL
  fel = FEL;
  #else
  fel = get_fel(i, j, k, Ps);
  #endif

  // thermal electron heating fraction //
  felth = fel*(1-felnth);

  // Update thermal electron entropy according to Ressler+ 2015 Eqn. 27, Use half-step variables //
  /* Leon: need discussion here */
  /* Leon: do we need to include adiabatic gain/loss by nonthermal electrons? */
  Pf[KEL] += (game-1.)/(gam-1.)*pow(Ps[RHO],gam-game)*felth*(ktotharm-ktotadv);

  // Update nonthermal electron entropy according to Chael 2017 Eqn. 30: //
  /* Leon: need discussion here */
  /* Leon: do we need to include adiabatic gain/loss by nonthermal electrons? */
  double Qtot = (pow(Ps[RHO],gam-1)/(gam-1)) * (Pf[RHO]*(qf.ucon[0])*(ktotharm-ktotadv)/Dt); // - duenth_adb[i][j][k]*qf.ucon[0]/Dt;
    
  // Diagnostics //
  struct of_state q;

  // ucon ucov bcon bcov //
  get_state(Ps, geom, &q);

  // proper time //
  double dtau = Dt/qf.ucon[0];

  /* Leon: make sure the heat injection is code unit consistent */
  /* Leon: shouldn't Dt here be Dtau? */
  inject_nonthermal(Pf, felnth*fel*Qtot, dtau);
  
  // heating rate //
  /* Leon: need to check what q it is */
  double uadv = ktotadv/(gam-1.)*pow(Pf[RHO],gam);
  double Qud = q.ucon[0]*q.ucov[0]*(Pf[UU] - uadv)*pow(Ps[RHO]/Pf[RHO],gam)/Dt;

  // viscous heating rate //
  // Leon: comeback to this later //
  Qvisc[i][j][k] = fel * Qud;

  // Reset total entropy //
  Pf[KTOT] = ktotharm;

}

/****************************************************************************************************
 *
 * @brief Calculates the electron heating fraction (fel) from eqn 48 in the Ressler 15 paper
 * 
 * @param i First grid index
 * @param j Second grid index
 * @param k Third grid index
 * @param P Primitive matrix at the specified gridpoint
 * @return double felnth
 * 
 ****************************************************************************************************/
double get_felnth(int i, int j, int k, double P[NVAR])
{

  /* define */
  double felnth;

  // constant fraction //
  #ifdef FELNTH
  felnth = FELNTH;
  return felnth;
  #endif

  // This function has access to all the primitives, so any custom fel_nth could be constructed here!
  felnth = 0.015; // This was the default constant value from Chael
    
  /* return */
  return felnth; 

}

/***************************************************************************************************
 * 
 * @brief General injection function. Can be easily modified to use custom injection distributions 
 * or even an injection calculated seperately for each bin
 * 
 * @param Pf - Primitives wish to udpate (by the operator +=)
 * @param Q - Heat to be injected (usually fel*felnth*Qtot)
 * @param Dt - Time step
 * Leon: the input dtau should be the proper time 
 * 
 ***************************************************************************************************/
void inject_nonthermal(double *Pf, double Q, double Dtau){
  // normterm is set in set_nonthermal_gammas and is = m_e*int(gam-1)*gam^-p //
  /* Make sure code unit consistent */
  //double plaw_norm = Q/normterm; 
  inject_nonthermal_plaw(Pf, Q, Dtau);
}

/***************************************************************************************************
 * 
 * @brief Injection of power law distributed nonthermal electrons (see Chael eq. 29)
 * 
 * @param Pf - primitive to be updated by the += operator
 * exit
 * Given timestep Dtau, powerlaw index p, and normalization constant C
 * Leon: the input dtau should be the proper time 
 *
 ***************************************************************************************************/
void inject_nonthermal_plaw(double *Pf, double normalization, double Dtau){

  /* variables */
  double rhs;
  double gammatemp;
  double Cfac = normalization/normterm;

  /* Make sure code unit consistent */
  NTEGAMMALOOP{

    /* Assign gamma */
    gammatemp = nteGammas[ig];

    /* within the range of gamma injection, also limited the rate of injection */
    if((gammatemp >= gammainjmin) && (gammatemp <= gammainjmax)) {
      rhs = Dtau*Cfac*pow(gammatemp,-PLAW);
      #if LIMIT_NTE
      /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
      if(Pf[ig + NTESTART] > 0) {
        rhs = MY_SIGN(rhs)*MY_MIN(fabs(rhs), fabs(0.2*Pf[ig + NTESTART]));
      }
      /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
      #endif
      Pf[ig + NTESTART] += rhs;
      if(Pf[ig + NTESTART] < 0) Pf[ig + NTESTART] = 0.0;
    }
      
  } 

}

/******************************************************************************/
/* do radiative cooling on nonthermal electron */
void cool_nonthermal(grid_prim_type Pi, grid_prim_type Ps, grid_prim_type Pf, double Dt){
  #pragma omp parallel for collapse(3) schedule(dynamic)
  ZLOOP {
    cool_nonthermal_zone(Pi[i][j][k], Ps[i][j][k], Pf[i][j][k], &(ggeom[i][j][CENT]), Dt);
  }
}

/******************************************************************************
 * 
 * @brief Applies radiative cooling to one zone
 * 
 * @param Pi Primitives to calculate the cooling from
 * @param Pf Primitives to apply the cooled result
 * @param geom Geometry in the zone being cooled
 * 
 ******************************************************************************/
void cool_nonthermal_zone(double *Pi, double *Ps, double *Pf, struct of_geom *geom, double Dt){
    
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // grid size in log gamma space //
  //double dg = log(nteGammas[1])-log(nteGammas[0]);
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  // gamma dot and number of gamma bins //
  double qdot_nth[NTEBINS], gdot_col[NTEBINS];
  double gdot[NTEBINS], ngammas[NTEBINS];
  double rhs, flux_up, flux_c; 
    
  // state vector //
  struct of_state q;

  // ucov ucon bcov bcon //
  get_state(Ps, geom, &q);

  // proper time //
  double dtau = Dt/q.ucon[0];
    
  // initialize //
  NTEGAMMALOOP{
    gdot[ig] = 0;
    ngammas[ig] = Ps[ig+NTESTART];
  } 

  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  // get gamma dot from radiation and coulomb collision //
  /* Leon: make sure gamma dot here is in code unit */
  /* gdot is negative, also, gdot inverse compton is omitted */
  calc_gdot_rad(Ps, geom, gdot, gdot_col);

  /* integrate, eq 38 in Chael 2017 */
  /* Leon: maybe use my formula of summation by part is better */
  /* coulomb collision between thermal and non-thermal electrons */
  NTEGAMMALOOP qdot_nth[ig] = -ME_CODE*pow(CL_CODE,2)*gdot_col[ig]*ngammas[ig];
  double qcnth = gamma_integral(qdot_nth); 
  qcnth = fabs(qcnth);
    
  /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  // Upwind derivative updates each n(gamma) bin
  NTEGAMMALOOP{
      
    /* centered flux */
    flux_c = gdot[ig]*ngammas[ig];

    /* upwind differencing */
    if (ig == NTEBINS-1) {
      flux_up = 0.0;
    } else {
      flux_up = gdot[ig+1]*ngammas[ig+1];
    }

    /* update, limit changes to 20% of the primitive */
    rhs = -dtau*(flux_up - flux_c)/logeBinSpace;
    rhs = rhs/nteGammas[ig];
    #if LIMIT_NTE
    /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
    if(Pf[ig + NTESTART] > 0) {
      rhs = MY_SIGN(rhs)*MY_MIN(fabs(rhs), fabs(0.2*Pf[ig+NTESTART]));
    }
    /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
    #endif
    Pf[ig+NTESTART] += rhs;
      
    /* floor the bins */
    if (Pf[ig+NTESTART] < 0.0) Pf[ig+NTESTART] = 0.0;
    
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //if(ig == NTEBINS-1){
    //  Pf[ig+NTESTART] -= dtau*(-gdot[ig]*ngammas[ig])/(dg*nteGammas[ig]);
    //}
    //else{
    //  Pf[ig+NTESTART] -= dtau*(gdot[ig+1]*ngammas[ig+1]-gdot[ig]*ngammas[ig])/(dg*nteGammas[ig]);
    //}
    //if (Pf[ig+NTESTART] < SMALL){
    //  Pf[ig+NTESTART] = 0;
    //}
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  }

  // Update thermal electrons with cooling nonthermal electrons //
  double ncool = -gdot[0]*ngammas[0];
  ncool = fabs(ncool);
  double qcool = ME_CODE*pow(CL_CODE,2)*(nteGammas[0]-1)*ncool;
  qcool = fabs(qcool);
    
  /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  /* chemical potential */
  //double chempot = calc_potential(Ps);
  /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
  // heat thermal electrons //
  /* Leon: make sure heat transfer through boundaries is in code unit */
  /* Leon: should also include radiative loss to the total entropy */
  //apply_thermal_heating(Pf, qcnth + qcool, ncool, dtau);
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
  //apply_thermal_heating(Pf, qcnth + qcool - chempot*ncool, dtau);
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
    
}

/***************************************************************************************************/
/* radiation injection */
/* Ps is the primitive we use to calculate the source term */
void calc_gdot_rad(double *Ps, struct of_geom *geom, double *gdot, double *gdot_col)
{

  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  // TODO: Does Bsq direction matter at all?
  // TODO: nion claculation
  // TODO: Inverse compton cooling is just a big ? not sure where to start...
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /* initialize just incase */
  NTEGAMMALOOP {
    gdot[ig] = 0.0;
    gdot_col[ig] = 0.0;
  }

  // Variable Declarations //
  struct of_state q;
  #if SYNCHROTRON_NTE
  double Bsq;
  #endif
  #if BREMSSTRAHLUNG_NTE || COULOMB_NTE
  double nion;
  #endif

  /* ucon ucov bcon bcov */
  get_state(Ps, geom, &q);
  
  /* Leon: all these rates should be in code unit */
  #if SYNCHROTRON_NTE
  /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
  #ifdef ART_SYNCH
  Bsq = pow(ART_SYNCH, 2.);
  #else
  Bsq = calc_bsq_cgs(Ps, geom);
  #endif
  //I'm temporarily hard coding in a constant B-field to test
  // Eq. 31 Chael + 17 
  NTEGAMMALOOP gdot[ig] += (-1.292e-11)*Bsq*pow(nteGammas[ig],2)/T_unit; 
  /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
  #endif

  #if BREMSSTRAHLUNG_NTE
  /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
  // This just assumes every ion is Hydrogen
  // Eq. 32 Chael + 17 
  nion = RHO_unit*Ps[RHO]/(MP + ME); // TODO: how to actually find nion... 
  NTEGAMMALOOP gdot[ig] += (-1.37e-16)*nion*nteGammas[ig]*(log(nteGammas[ig]) + 0.36)/T_unit; 
  /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
  #endif

  #if COULOMB_NTE
  /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
  /* here, ni = ne */
  nion = RHO_unit*Ps[RHO]/(MP+ME);
  NTEGAMMALOOP gdot_col[ig] = (-1.491e-14)*nion*(log(nteGammas[ig]) + log(nion) + 74.7)/T_unit; 
  NTEGAMMALOOP gdot[ig] += gdot_col[ig];
  /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
  #endif

  #if COMPTON_NTE
  /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
  // This one is a bit trickier... Need some help here on how to calulcate Trad and Ehat. 
  // See notes for thoughts
  /////////////////////////////////////////////////////////////////////////////
  // NTEGAMMALOOP gdot[ig] += (-3.25e-8)*(Ehat)*pow(nteGammas[ig],2)*FKN[ig];
  /////////////////////////////////////////////////////////////////////////////
  /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
  #endif

}

/*************************************************************************
 * 
 * @brief Finds the magnitude of the magnetic field in cgs units
 * 
 * @param Ps Primitives in the desired zone used to calculate magnetic field
 * @param geom Geomtry in the desired zone
 * @return double B^2 in cgs
 * 
 *************************************************************************/
double calc_bsq_cgs(double *Ps, struct of_geom *geom){

  // variables //
  double Bp[NDIM], Vcon[NDIM], Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
  double Vfac, VdotV, UdotBp;

  // primitive //
  Bp[1] = Ps[B1]*B_unit;
  Bp[2] = Ps[B2]*B_unit;
  Bp[3] = Ps[B3]*B_unit;
  Vcon[1] = Ps[U1];
  Vcon[2] = Ps[U2];
  Vcon[3] = Ps[U3];

  // Get Ucov //
  VdotV = 0.;
  for(int l = 1; l < NDIM; l++) {
    for(int m = 1; m < NDIM; m++) {
      VdotV += geom->gcov[l][m]*Vcon[l]*Vcon[m];
    }
  }

  // ucon ucov //
  Vfac = sqrt(-1./geom->gcon[0][0]*(1. + fabs(VdotV)));
  Ucon[0] = -Vfac*geom->gcon[0][0];
  for(int l = 1; l < NDIM; l++) 
    Ucon[l] = Vcon[l] - Vfac*geom->gcon[0][l];
  lower(Ucon, geom->gcov, Ucov);

  // Get Bcon, Bcov, and B //
  UdotBp = 0.;
  for(int l = 1; l < NDIM; l++)
    UdotBp += Ucov[l]*Bp[l];
  Bcon[0] = UdotBp;
  for(int l = 1; l < NDIM; l++)
    Bcon[l] = (Bp[l] + Ucon[l]*UdotBp)/Ucon[0];
  lower(Bcon, geom->gcov, Bcov);

  // return
  return dot(Bcon,Bcov);
}

/***************************************************************************************************/
/* heat the thermal electron */
/* note that Dtau here is the proper time */
/* Pf is the primitive wish to update by the += operator */
void apply_thermal_heating(double *Pf, double heat, double nrate, double Dtau){

  // get electron internal energy from entropy //
  double ue_f = Pf[KEL]*pow(Pf[RHO],game)/(game-1.);

  // Update electron internal energy via heat //
  /* Leon: make sure heat is in code unit */
  /* Leon: caution that number density is not changed */
  ue_f += heat*Dtau;

  // Update electron entropy
  Pf[KEL] = (game-1.)*ue_f*pow(Pf[RHO],-game);
    
}

/***************************************************************************************************/
/* Chemical potentials: Modified by leon to implement the newtonian form */
/* Leon: not useful right now */
double calc_potential(double *Ps){
  double rho = Ps[RHO];
  double kappae = Ps[KEL];
  double thetae = (MP/ME)*kappae*pow(rho,game-1.);
  double chem_po = ME_CODE*pow(CL_CODE,2)*thetae/(game - 1.0)*(game - log(kappae));
  return chem_po;
}

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
#endif