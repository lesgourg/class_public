/**
 * Small wrapper file for the hyrec code to be used in thermodynamics
 * Nils Schoeneberg Feb 2019
 * */

#include "common.h"
#include "thermodynamics.h"
#include "wrap_hyrec.h"

//(phy->hyrec_Alpha_inf_file)
/* Boundaries and number of elements of temperature tables */
#define TR_MIN 0.004            /* Tr parameters */
#define TR_MAX 0.4
#define NTR    100
#define TM_TR_MIN 0.1           /* Tm/Tr parameters */
#define TM_TR_MAX 1.0
#define NTM 40

int thermodynamics_hyrec_rec_1Hs_post_saha(struct thermohyrec* phy, int iz_out, double z_out, double xHeII,
                                           double H, double TR, double nH, double* xH1s);
int thermodynamics_hyrec_readrates(struct thermohyrec* phy);
int thermodynamics_hyrec_readtwogparams(struct thermohyrec* phy);

int thermodynamics_hyrec_init(struct precision* ppr, double Nnow, double T_cmb, double fHe, struct thermohyrec* phy){

  if(phy->thermohyrec_verbose > 0){
    printf(" -> Using the hyrec wrapper programmed by Nils Sch. (Feb2019)\n");
    printf("    implements HyRec version Oct2012 by Yacine Ali-Haimoud and Chris Hirata\n");
  }

  int index_virt,index_ly,iz;
  phy->N_LY = 3;
  phy->N_VIRT = NVIRT;

  strcpy(phy->alpha_file,ppr->hyrec_Alpha_inf_file);
  strcpy(phy->rr_file,ppr->hyrec_R_inf_file);
  strcpy(phy->twog_file,ppr->hyrec_two_photon_tables_file);

  phy->n_H_now = Nnow;
  phy->nH0 = phy->n_H_now; //TODO :: units?
  phy->T_cmb = T_cmb;
  phy->T0 = phy->T_cmb; //TODO :: units?
  phy->fHe = fHe;

  //At what values should the fminus and fplus histories be built?
  phy->zstart = 8060.; //Make sure hyrec filled some values before first call from class (at z=8050)
  phy->zend = 0.;
  phy->dlna = 8.49e-5;
  //phy->dlna = 1e-3;

  phy->nz = (long) floor(2.+log((1.+phy->zstart)/(1.+phy->zend))/phy->dlna);



  class_alloc(phy->rate_table,
              sizeof(HRATEEFF),
              phy->error_message);
  class_alloc(phy->twog_params,
              sizeof(TWO_PHOTON_PARAMS),
              phy->error_message);

  phy->rate_table->logTR_tab = create_1D_array(NTR);
  phy->rate_table->TM_TR_tab = create_1D_array(NTM);
  phy->rate_table->logAlpha_tab[0] = create_2D_array(NTM, NTR);
  phy->rate_table->logAlpha_tab[1] = create_2D_array(NTM, NTR);
  phy->rate_table->logR2p2s_tab = create_1D_array(NTR);

  class_call(thermodynamics_hyrec_readrates(phy),
             phy->error_message,
             phy->error_message);

  class_call(thermodynamics_hyrec_readtwogparams(phy),
             phy->error_message,
             phy->error_message);

  /* history of x_e and T_mat ( these are totally not needed, but convenience variables :D )*/
  phy->xe_output    = create_1D_array(phy->nz);
  //Tm_output         = create_1D_array(phy->nz);
  /* history of photon occupation numbers */
  phy->Dfnu_hist    = create_2D_array(NVIRT, phy->nz);
  phy->Dfminus_hist = create_2D_array(NVIRT, phy->nz);
  phy->Dfminus_Ly_hist[0] = create_1D_array(phy->nz);   /* Ly-alpha */
  phy->Dfminus_Ly_hist[1] = create_1D_array(phy->nz);   /* Ly-beta  */
  phy->Dfminus_Ly_hist[2] = create_1D_array(phy->nz);   /* Ly-gamma */

  /* Make sure the input spectrum is initialized at zero */
  /*for (iz=0; iz<phy->nz; iz++){
    phy->Dfminus_Ly_hist[0][iz] = 0;
    phy->Dfminus_Ly_hist[1][iz] = 0;
    phy->Dfminus_Ly_hist[2][iz] = 0;
  }*/


  //Ratios of fine structure constant and electron mass at recombination compared to now
  phy->fsR = 1.;
  phy->meR = 1.;
  phy->xHeIII = phy->fHe;   /* Delta_xe = xHeIII here */
  phy->xHeII = 0.; //Uninitialized

  phy->izH0 = (long) floor(1 + log(_k_B_/_eV_*phy->T_cmb/phy->fsR/phy->fsR/phy->meR*(1.+phy->zstart)/(TR_MAX*0.99))/phy->dlna);
  //When T_gamma reaches TR_MAX (just approximately!), 1% error margin
  phy->zH0  = (1.+phy->zstart)*exp(-phy->izH0 * phy->dlna) - 1.;

  phy->stage = 0;
  phy->saha_flag = _TRUE_;

  phy->filled_until_index_z = 0;
  /*for(index_virt=0;index_virt<phy->N_VIRT;++index_virt){
    phy->logfminus_hist[index_virt][0] = 0.0;
  }
  for(index_ly=0;index_ly<phy->N_LY;++index_ly){
    phy->logfminus_Ly_hist[index_ly][0] = 0.0;
  }*/


  phy->TR_prev = phy->T_cmb*(1+phy->zstart);
  phy->TM_prev = phy->TR_prev;
  phy->z_prev = (1+phy->zstart);
  phy->xHeIII_limit = 1e-8;


  return _SUCCESS_;
}

int thermodynamics_hyrec_free(struct thermohyrec* phy){

  int index_ly,index_virt;
  for(index_ly=0;index_ly<phy->N_LY;++index_ly){
    free(phy->Dfminus_Ly_hist[index_ly]);
  }
  for(index_virt=0;index_virt<phy->N_VIRT;++index_virt){
    free(phy->Dfminus_hist[index_virt]);
    free(phy->Dfnu_hist[index_virt]);
  }
  free(phy->Dfminus_hist);
  free(phy->Dfnu_hist);

  free(phy->rate_table);
  free(phy->twog_params);
  free(phy->xe_output);

  return _SUCCESS_;
}

int thermodynamics_hyrec_get_xe(struct thermohyrec * phy,
                                double z, double H, double T_b, double T_gamma,
                                double* x_e, double* dxe_dlna, double energy_injection) {

  /** Define local variables */
  int iz;
  /* For a step*/
  int iz_in,iz_out;
  double z_in,z_out;
  /* For the whole call */
  int iz_goal;
  double z_goal;
  /* For the final interpolation */
  double z_goalm1;

  /* Intermediate quantities that are used at some points locally */
  double xHeIISaha,dxHeIISaha_dlna,xH1s_p,xH1s_m,Dxe,DdxHeIIdlna_Dxe;
  double xH1s_in, xHeII_in;

  /* For any interpolation */
  double frac;

  /* Interpolated and exact quantities at EACH step */
  double nH,TR,TM,xe_in;

  /* Something related to switching off modes or not depending on pion */
  double Pion_TR_rescaled,Pion_RLya,Pion_four_betaB,Pion;

  /** Calculate the quantities until which the table should be extended */
  iz_goal = (int)ceil(-log((1+z)/(1.+phy->zstart))/phy->dlna);
  z_goal = (1.+phy->zstart)*exp(-phy->dlna*(iz_goal)) - 1.;
  if(z_goal<0.){z_goal=0.;}

  /** Only add new indices if that is really required */
  if(iz_goal>phy->filled_until_index_z){

    if(phy->thermohyrec_verbose > 2){printf("Filling [%i,%i]\n",phy->filled_until_index_z+1,iz_goal);}

    for(iz=phy->filled_until_index_z+1;iz<=iz_goal;++iz){
      iz_in = iz-1;
      iz_out = iz;
      z_in = (1.+phy->zstart)*exp(-phy->dlna*(iz_in)) - 1.;
      z_out = (1.+phy->zstart)*exp(-phy->dlna*(iz_out)) - 1.;

      nH = 1e-6*phy->nH0 * (1.+z_in)*(1.+z_in)*(1.+z_in);
      frac = ((1+z_in)-(1+phy->z_prev))/((1+z)-(1+phy->z_prev));
      TR = T_gamma * _k_B_/_eV_ * frac + (1.-frac)*phy->TR_prev;
      TM = T_b * _k_B_/_eV_ * frac + (1.-frac)*phy->TM_prev;





      if(phy->stage == 0){
        /**
         * Stage 0 : He III -> II
         *  This stage is calculated in saha equilibrium.
         *  ASSUMES T_gamma == T_matter.
         **/

        if(phy->xHeIII > phy->xHeIII_limit){
          phy->xe_output[iz_out] = rec_xesaha_HeII_III(phy->nH0, phy->T0, phy->fHe, z_out, &(phy->xHeIII), phy->fsR, phy->meR);
          if(phy->thermohyrec_verbose > 2){printf("[%i] : %.10e \n",iz_out,phy->xe_output[iz_out]);}
        }
        else{
          /* Switch to next stage, and initialize its variables. */
          phy->stage++;
          if(phy->thermohyrec_verbose > 2){printf("Next stage %i \n",phy->stage);}
          phy->dxHeIIdlna_prev2 = (phy->xe_output[iz_out-2] - phy->xe_output[iz_out-4])/(2.*phy->dlna);
          phy->dxHeIIdlna_prev  = (phy->xe_output[iz_out-1] - phy->xe_output[iz_out-3])/(2.*phy->dlna);
          phy->xHeII = rec_saha_xHeII(phy->nH0, phy->T0, phy->fHe, z_in, phy->fsR, phy->meR);
          phy->saha_flag = _TRUE_; /* Start with post-saha expansion */
        }

      }
      if(phy->stage == 1){
        /**
         * Stage 1 : He II -> I recombination.
         *  Hydrogen in Saha equilibrium with the free electrons.
         *  Tm fixed to steady state.
         *  Integrate until TR is low enough that can start integrating hydrogen recombination
         *  (this occurs at index izH0 computed in rec_get_cosmoparam).
         * Start with post-Saha expansion.
         **/
        if(iz_out<phy->izH0+1){
          //rec_get_xe_next1_He
          //////////////////////////////////
          phy->xH1s       = rec_saha_xH1s(phy->xHeII, phy->nH0, phy->T0, z_in, phy->fsR, phy->meR);
          phy->dxHeIIdlna = rec_helium_dxHeIIdlna(phy->xH1s, phy->xHeII, phy->nH0, phy->T0, phy->fHe, H, z_in, phy->fsR, phy->meR);

          /* Post-Saha approximation during the early phase of HeII->HeI recombination */
          if (phy->saha_flag == _TRUE_) {
              xHeIISaha = rec_saha_xHeII(phy->nH0, phy->T0, phy->fHe, z_out, phy->fsR, phy->meR);
              dxHeIISaha_dlna  = (1.+z_out)*(rec_saha_xHeII(phy->nH0, phy->T0, phy->fHe, z_out-0.5, phy->fsR, phy->meR)
                                            -rec_saha_xHeII(phy->nH0, phy->T0, phy->fHe, z_out+0.5, phy->fsR, phy->meR));

              Dxe    = 0.01*(phy->fHe - xHeIISaha);
              xH1s_p = rec_saha_xH1s(xHeIISaha+Dxe, phy->nH0, phy->T0, z_out, phy->fsR, phy->meR);
              xH1s_m = rec_saha_xH1s(xHeIISaha-Dxe, phy->nH0, phy->T0, z_out, phy->fsR, phy->meR);

              DdxHeIIdlna_Dxe  = (rec_helium_dxHeIIdlna(xH1s_p, xHeIISaha+Dxe, phy->nH0, phy->T0, phy->fHe, H, z_out, phy->fsR, phy->meR)
                                 -rec_helium_dxHeIIdlna(xH1s_m, xHeIISaha-Dxe, phy->nH0, phy->T0, phy->fHe, H, z_out, phy->fsR, phy->meR))/(2.*Dxe);

              phy->xHeII = xHeIISaha + dxHeIISaha_dlna/DdxHeIIdlna_Dxe;

              /* Check that the post-Saha approximation is still valid. If not, switch it off for future iterations */
              if (fabs(phy->xHeII - xHeIISaha) > DXHEII_MAX){phy->saha_flag = _FALSE_;}
          }
          /* Otherwise integrate ODE */
          else{
            phy->xHeII += phy->dlna * (1.25 * phy->dxHeIIdlna - 0.25 * phy->dxHeIIdlna_prev2);
          }

          /* Update stored derivatives */
          phy->dxHeIIdlna_prev2 = phy->dxHeIIdlna_prev;
          phy->dxHeIIdlna_prev  = phy->dxHeIIdlna;
          //////////////////////////////////

          phy->xH1s = rec_saha_xH1s(phy->xHeII, phy->nH0, phy->T0, z_out, phy->fsR, phy->meR);
          phy->xe_output[iz_out] = (1.-phy->xH1s) + phy->xHeII;
          if(phy->thermohyrec_verbose > 2){printf("[%i] : %.10e \n",iz_out,phy->xe_output[iz_out]);}
        }
        else{
          /* Switch to next stage, and initialize its variables. */
          phy->stage++;
          if(phy->thermohyrec_verbose > 2){printf("Next stage %i \n",phy->stage);}

          phy->dxHIIdlna_prev2 = (phy->xe_output[iz_out-2] - phy->xe_output[iz_out-4])/(2.*phy->dlna) - phy->dxHeIIdlna_prev2;
          phy->dxHIIdlna_prev  = (phy->xe_output[iz_out-1] - phy->xe_output[iz_out-3])/(2.*phy->dlna) - phy->dxHeIIdlna_prev;
          phy->saha_flag       = _TRUE_;
        }

      }
      if(phy->stage == 2){
        /** H II -> I and He II -> I simultaneous recombination (rarely needed but just in case)
              Tm fixed to steady state.
              Integrate H and He simultaneously until xHeII < XHEII_MIN
              Start with post-saha expansion for hydrogen */

        if(phy->xHeII > XHEII_MIN){
          //get_rec_next2_HHe
          //////////////////////////////
          xH1s_in  = phy->xH1s;
          xHeII_in = phy->xHeII;
          xe_in    = xHeII_in + (1.-xH1s_in);

          /* Evolve HeII by solving ODE */
          phy->dxHeIIdlna  = rec_helium_dxHeIIdlna(xH1s_in, xHeII_in, phy->nH0, phy->T0, phy->fHe, H, z_in, phy->fsR, phy->meR);
          phy->xHeII += phy->dlna * (1.25 * phy->dxHeIIdlna - 0.25 * phy->dxHeIIdlna_prev2);

          /* Compute Hydrogen derivative at input time. Even if using the post-Saha expansion, needed for getting the correct radiation field at z_in */
          phy->dxHIIdlna = rec_dxHIIdlna(MODEL, xe_in, (1.-xH1s_in), nH, H, TM, TR, phy->rate_table, phy->twog_params,
                                    phy->Dfminus_hist, phy->Dfminus_Ly_hist, phy->Dfnu_hist, phy->zH0, iz_in-phy->izH0, z_in, phy->fsR, phy->meR);

          /* If Hydrogen is still close to Saha equilibrium do a post-Saha expansion for Hydrogen */
          if(phy->saha_flag == _TRUE_){
            class_call(thermodynamics_hyrec_rec_1Hs_post_saha(phy, iz_out, z_out, phy->xHeII, H, TR, nH, &(phy->xH1s)),
                       phy->error_message,
                       phy->error_message);
          }
          /* Otherwise solve HII ODE */
          else{
            phy->xH1s -= phy->dlna * (1.25 * phy->dxHIIdlna - 0.25 * phy->dxHIIdlna_prev2);
          }

          /* Update derivatives */
          phy->dxHIIdlna_prev2  = phy->dxHIIdlna_prev;
          phy->dxHIIdlna_prev   = phy->dxHIIdlna;
          phy->dxHeIIdlna_prev2 = phy->dxHeIIdlna_prev;
          phy->dxHeIIdlna_prev  = phy->dxHeIIdlna;

          //////////////////////////////
          phy->xe_output[iz_out] = (1.-phy->xH1s) + phy->xHeII;
          if(phy->thermohyrec_verbose > 2){printf("[%i] : %.10e \n",iz_out,phy->xe_output[iz_out]);}
        }
        else{
          /* Switch to next stage, and initialize its variables. */
          phy->stage++;
          if(phy->thermohyrec_verbose > 2){printf("Next stage %i \n",phy->stage);}
        }

      }

      /** Check if pion is low enough if stage is 3 or 4 */
      if(phy->stage == 3 || phy->stage==4){
        Pion_TR_rescaled = TR/phy->fsR/phy->fsR/phy->meR; // TR rescaled for alpha, me
        Pion_RLya        = LYA_FACT(phy->fsR, phy->meR) * H / nH;
        Pion_four_betaB  = SAHA_FACT(phy->fsR, phy->meR) *Pion_TR_rescaled*sqrt(Pion_TR_rescaled) *exp(-0.25*EI/Pion_TR_rescaled) * alphaB_PPB(Pion_TR_rescaled, phy->fsR, phy->meR);
        Pion             = Pion_four_betaB/(3.*Pion_RLya + L2s_rescaled(phy->fsR, phy->meR) + Pion_four_betaB);
        if(phy->thermohyrec_verbose > 3){printf("Calculated pion value %.10e / %.10e ",Pion,PION_MAX);}
      }

      if(phy->stage == 3){
        /** H recombination. Helium assumed entirely neutral.
               Tm fixed to steady-state until its relative difference from Tr is DLNT_MAX */
        if(1.-T_b/T_gamma < DLNT_MAX*0.99){ //1% error margin
          //rec_get_xe_next1_H
          ////////////////////////////


          /* Switch off radiative transfer calculation if needed (param->nzrt and izH0 are evaluated in rec_get_cosmoparam) */
          int model = (Pion < PION_MAX || MODEL != FULL) ? MODEL : EMLA2s2p;
          if(phy->thermohyrec_verbose > 3){printf("[%i,%i] : Model %i (%i)\n",iz_out,phy->stage,model,EMLA2s2p);}

          xe_in = phy->xe_output[iz_in];
          phy->dxHIIdlna = rec_dxHIIdlna(model, xe_in, xe_in, nH, H, TM, TR, phy->rate_table, phy->twog_params,
                                  phy->Dfminus_hist, phy->Dfminus_Ly_hist, phy->Dfnu_hist, phy->zH0, iz_in-phy->izH0, z_in, phy->fsR, phy->meR);

          /* If close to Saha equilibrium (with xHeII = 0), do a post-Saha expansion */
          if (phy->saha_flag == _TRUE_) {
            class_call(thermodynamics_hyrec_rec_1Hs_post_saha(phy, iz_out, z_out, 0., H, TR, nH, &(phy->xH1s)),
                       phy->error_message,
                       phy->error_message);
            phy->xe_output[iz_out] = 1.-phy->xH1s;
          }
          /* Otherwise evolve ODE */
          else{
            phy->xe_output[iz_out] = xe_in + phy->dlna * (1.25 * phy->dxHIIdlna - 0.25 * phy->dxHIIdlna_prev2);
          }
          /* Update previous derivatives */
          phy->dxHIIdlna_prev2 = phy->dxHIIdlna_prev;
          phy->dxHIIdlna_prev  = phy->dxHIIdlna;
          ////////////////////////////
          if(phy->thermohyrec_verbose > 2){printf("[%i] : %.10e \n",iz_out,phy->xe_output[iz_out]);}
        }
        else{
          /* Switch to next stage, and initialize its variables. */
          phy->stage++;
          if(phy->thermohyrec_verbose > 2){printf("Next stage %i \n",phy->stage);}
        }

      }
      if(phy->stage == 4){
        /** Evolve xe and Tm simultaneously until the lower bounds of integration tables are reached.
              Note that the radiative transfer calculation is switched off automatically in the functions
              rec_get_xe_next1_H and rec_get_xe_next2_HTm when it is no longer relevant. */
        if(TR/phy->fsR/phy->fsR/phy->meR > TR_MIN && TM/TR > TM_TR_MIN){

          //rec_get_xe_next2_HTm
          //////////////////////////////
          /* Switch off radiative transfer calculation if needed */
          int model = (Pion > PION_MAX || MODEL != FULL) ? MODEL : EMLA2s2p;
          if(phy->thermohyrec_verbose > 3){printf("[%i,%i] : Model %i (%i)\n",iz_out,phy->stage,model,EMLA2s2p);}

          xe_in = phy->xe_output[iz_in];
          phy->dxHIIdlna = rec_dxHIIdlna(model, xe_in, xe_in, nH, H, TM, TR, phy->rate_table, phy->twog_params,
                             phy->Dfminus_hist, phy->Dfminus_Ly_hist, phy->Dfnu_hist, phy->zH0, iz_in-phy->izH0, z_in, phy->fsR, phy->meR);

          phy->xe_output[iz_out] = xe_in + phy->dlna * (1.25 * phy->dxHIIdlna - 0.25 * phy->dxHIIdlna_prev2);

          phy->dxHIIdlna_prev2 = phy->dxHIIdlna_prev;
          phy->dxHIIdlna_prev  = phy->dxHIIdlna;
          /////////////////////////////////////
          if(phy->thermohyrec_verbose > 2){printf("[%i] : %.10e \n",iz_out,phy->xe_output[iz_out]);}
        }
        else{
          /* Switch to next stage, and initialize its variables. */
          phy->stage++;
          if(phy->thermohyrec_verbose > 2){printf("Next stage %i \n",phy->stage);}
        }

      }
      if(phy->stage == 5){
        /** For low redshifts (z < 20 or so) use Peeble's model (Tm is evolved with xe).
            The precise model does not metter much here as
            1) the free electron fraction is basically zero (~1e-4) in any case and
            2) the universe is going to be reionized around that epoch */
        //rec_get_xe_next2_HTm
        //////////////////////////////
        xe_in = phy->xe_output[iz_in];
        phy->dxHIIdlna = rec_dxHIIdlna(PEEBLES, xe_in, xe_in, nH, H, TM, TR, phy->rate_table, phy->twog_params,
                           phy->Dfminus_hist, phy->Dfminus_Ly_hist, phy->Dfnu_hist, phy->zH0, iz-phy->izH0, z_in, phy->fsR, phy->meR);

        phy->xe_output[iz_out] = xe_in + phy->dlna * (1.25 * phy->dxHIIdlna - 0.25 * phy->dxHIIdlna_prev2);

        phy->dxHIIdlna_prev2 = phy->dxHIIdlna_prev;
        phy->dxHIIdlna_prev  = phy->dxHIIdlna;
        /////////////////////////////////////
        if(phy->thermohyrec_verbose > 2){printf("[%i] : %.10e \n",iz_out,phy->xe_output[iz_out]);}
      }

    }//End for iz not filled

    phy->filled_until_index_z = iz_goal;
    phy->TM_prev = TM;
    phy->TR_prev = TR;
    phy->z_prev = z;
  }
  z_goalm1 = (1.+phy->zstart)*exp(-phy->dlna*(iz_goal-1)) - 1.;
  /**
   * If the table is filled already, or after filling it to the desired value
   *
   * We want to take the values at iz_goal and iz_goal-1 to calculate
   *  the free electron fraction x_e at this position
   * Also obtain logarithmic derivative delta_xe/delta_lna
   *  by numerical derivative
   **/

  frac = ((1.+z)-(1.+z_goalm1))/((1.+z_goal)-(1.+z_goalm1));
  *x_e = frac * phy->xe_output[iz_goal] + (1.-frac)* phy->xe_output[iz_goal-1];
  *dxe_dlna = (phy->xe_output[iz_goal] - phy->xe_output[iz_goal-1])/(phy->dlna);

  if(phy->thermohyrec_verbose > 2){
    //printf("[%i] : %.10e , [%i] : %.10e -> %.10e \n",iz_goal,phy->xe_output[iz_goal],iz_goal-1,phy->xe_output[iz_goal-1],*dxedlna);
    printf("[%i,%.10e] : %.10e , [%i,%.10e] : %.10e -> [%.10e] : %.10e \n",iz_goal,z_goal,phy->xe_output[iz_goal],iz_goal-1,z_goalm1,phy->xe_output[iz_goal-1],z,*x_e);
  }

  return _SUCCESS_;
}

int thermodynamics_hyrec_rec_1Hs_post_saha(struct thermohyrec* phy, int iz_out, double z_out, double xHeII,
                                           double H, double TR, double nH, double* xH1s){

  double xH1sSaha, xHIISaha, dxH1sSaha_dlna, dxH1sdlna_Saha, DdxH1sdlna_DxH1s, Dxe;

  xH1sSaha = rec_saha_xH1s(xHeII, phy->nH0, phy->T0, z_out, phy->fsR, phy->meR);
  xHIISaha = 1.-xH1sSaha;

  dxH1sSaha_dlna = (1.+z_out)*(rec_saha_xH1s(xHeII, phy->nH0, phy->T0, z_out-0.5, phy->fsR, phy->meR)
                              -rec_saha_xH1s(xHeII, phy->nH0, phy->T0, z_out+0.5, phy->fsR, phy->meR));
                                  /* (partial xHII)/(partial lna). Use xH1s = 1-xHII for better accuracy. */
  if (xHeII != 0.){
    //xHe == 0.0 will be passed EXACTLY when we are in the HII -> HI case, where HeII -> HeI has already finished completely
     Dxe = 0.01*(phy->fHe - xHeII);
     dxH1sSaha_dlna += (rec_saha_xH1s(xHeII+Dxe, phy->nH0, phy->T0, z_out, phy->fsR, phy->meR)
                        -rec_saha_xH1s(xHeII-Dxe, phy->nH0, phy->T0, z_out, phy->fsR, phy->meR))/(2.*Dxe)
                       *rec_helium_dxHeIIdlna(xH1sSaha, xHeII, phy->nH0, phy->T0, phy->fHe, H, z_out, phy->fsR, phy->meR);
                          /* (partial xHII)/(partial xHeII).dxHeII/dlna */
  }

  dxH1sdlna_Saha = -rec_dxHIIdlna(MODEL, xHIISaha + xHeII, xHIISaha, nH, H, TR, TR, phy->rate_table, phy->twog_params,
                                  phy->Dfminus_hist, phy->Dfminus_Ly_hist, phy->Dfnu_hist, phy->zH0, iz_out-phy->izH0, z_out, phy->fsR, phy->meR);
  Dxe            = 0.01*xH1sSaha;
  DdxH1sdlna_DxH1s = (rec_dxHIIdlna(MODEL, xHIISaha+Dxe + xHeII, xHIISaha+Dxe, nH, H, TR, TR, phy->rate_table, phy->twog_params,
                                  phy->Dfminus_hist, phy->Dfminus_Ly_hist, phy->Dfnu_hist, phy->zH0, iz_out-phy->izH0, z_out, phy->fsR, phy->meR)
                    -rec_dxHIIdlna(MODEL, xHIISaha-Dxe + xHeII, xHIISaha-Dxe, nH, H, TR, TR, phy->rate_table, phy->twog_params,
                                  phy->Dfminus_hist, phy->Dfminus_Ly_hist, phy->Dfnu_hist, phy->zH0, iz_out-phy->izH0, z_out, phy->fsR, phy->meR))/(2.*Dxe);

  *xH1s = xH1sSaha + (dxH1sSaha_dlna - dxH1sdlna_Saha)/DdxH1sdlna_DxH1s;

  /* Check that we are still close enough to Saha equilibrium. If not, switch post-saha expansion off */
  if (fabs(*xH1s - xH1sSaha) > DXHII_MAX){phy->saha_flag = _FALSE_;}

  return _SUCCESS_;
}

int thermodynamics_hyrec_readrates(struct thermohyrec* phy){
  FILE *fA;
  FILE *fR;

  HRATEEFF *rate_table = phy->rate_table;
  class_open(fA,phy->alpha_file, "r",phy->error_message);
  class_open(fR,phy->rr_file, "r",phy->error_message);

  /* Code below is essentially copied from hyrec (hence the weird indents and all)*/
   unsigned i, j, l;
   int fscanf_result;

   maketab(log(TR_MIN), log(TR_MAX), NTR, rate_table->logTR_tab);
   maketab(TM_TR_MIN, TM_TR_MAX, NTM, rate_table->TM_TR_tab);
   rate_table->DlogTR = rate_table->logTR_tab[1] - rate_table->logTR_tab[0];
   rate_table->DTM_TR = rate_table->TM_TR_tab[1] - rate_table->TM_TR_tab[0];

   for (i = 0; i < NTR; i++) {
      for (j = 0; j < NTM; j++) {
	 for (l = 0; l <= 1; l++) {
           fscanf_result = fscanf(fA, "%le", &(rate_table->logAlpha_tab[l][j][i]));
           if(fscanf_result != 1){printf("Hyrec Warning :: Could not read log Alpha table (Alpha_inf.dat)");}
           rate_table->logAlpha_tab[l][j][i] = log(rate_table->logAlpha_tab[l][j][i]);
        }
      }

      fscanf_result = fscanf(fR, "%le", &(rate_table->logR2p2s_tab[i]));
      if(fscanf_result != 1){printf("Hyrec Warning :: Could not read rate table (R_inf.dat)");}
      rate_table->logR2p2s_tab[i] = log(rate_table->logR2p2s_tab[i]);

   }
   fclose(fA);
   fclose(fR);

  /* End of hyrec code piece*/
  return _SUCCESS_;
}

int thermodynamics_hyrec_readtwogparams(struct thermohyrec* phy){

  FILE *fA = NULL;
  unsigned b;
  double L2s1s_current;
  int fscanf_result;

  TWO_PHOTON_PARAMS *twog = phy->twog_params;
  class_open(fA,phy->twog_file, "r" ,phy->error_message);

  /* Code below is essentially copied from hyrec (hence the weird indents and all)*/
   for (b = 0; b < NVIRT; b++) {
      fscanf_result = 0;
      fscanf_result += fscanf(fA, "%le", &(twog->Eb_tab[b]));
      fscanf_result += fscanf(fA, "%le", &(twog->A1s_tab[b]));
      fscanf_result += fscanf(fA, "%le", &(twog->A2s_tab[b]));
      fscanf_result += fscanf(fA, "%le", &(twog->A3s3d_tab[b]));
      fscanf_result += fscanf(fA, "%le", &(twog->A4s4d_tab[b]));
      if(fscanf_result!=5){printf("Hyrec Warning :: Could not read Two Photon table (two_photon_tables.dat)");}
   }
   fclose(fA);

   /* Normalize 2s--1s differential decay rate to L2s1s (can be set by user in hydrogen.h) */
   L2s1s_current = 0.;
   for (b = 0; b < NSUBLYA; b++) L2s1s_current += twog->A2s_tab[b];
   for (b = 0; b < NSUBLYA; b++) twog->A2s_tab[b] *= L2s1s/L2s1s_current;


  /* Switches for the various effects considered in Hirata (2008) and diffusion:
      Effect A: correct 2s-->1s rate, with stimulated decays and absorptions of non-thermal photons
      Effect B: Sub-Lyman-alpha two-photon decays
      Effect C: Super-Lyman-alpha two-photon decays
      Effect D: Raman scattering */

   #if (EFFECT_A == 0)
     for (b = 0; b < NSUBLYA; b++) twog->A2s_tab[b] = 0;
   #endif
   #if (EFFECT_B == 0)
     for (b = 0; b < NSUBLYA; b++) twog->A3s3d_tab[b] = twog->A4s4d_tab[b] = 0;
   #endif
   #if (EFFECT_C == 0)
      for (b = NSUBLYA; b < NVIRT; b++) twog->A3s3d_tab[b] = twog->A4s4d_tab[b] = 0;
   #endif
   #if (EFFECT_D == 0)
      for (b = NSUBLYA; b < NVIRT; b++) twog->A2s_tab[b] = 0;
      for (b = NSUBLYB; b < NVIRT; b++) twog->A3s3d_tab[b] = 0;
   #endif
   #if (DIFFUSION == 0)
      for (b = 0; b < NVIRT; b++) twog->A1s_tab[b] = 0;
   #endif

  /* End of hyrec code piece*/
  return _SUCCESS_;
}
