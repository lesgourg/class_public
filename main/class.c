#include <math.h>
/** @file class.c
 * Julien Lesgourgues, 17.04.2011
 */

#include "class.h"

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermodynamics th;           /* for thermodynamics */
  struct perturbations pt;         /* for source functions */
  struct primordial pm;       /* for primordial spectra */
  struct fourier fo;        /* for non-linear spectra */
  struct transfer tr;        /* for transfer functions */
  struct harmonic hr;          /* for output spectra */
  struct lensing le;          /* for lensed spectra */
  struct distortions sd;      /* for spectral distortions */
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */

  if (input_init(argc, argv,&pr,&ba,&th,&pt,&tr,&pm,&hr,&fo,&le,&sd,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init \n=>%s\n",errmsg);
    return _FAILURE_;
  }

  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (perturbations_init(&pr,&ba,&th,&pt) == _FAILURE_) {
    printf("\n\nError in perturbations_init \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (primordial_init(&pr,&pt,&pm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (fourier_init(&pr,&ba,&th,&pt,&pm,&fo) == _FAILURE_) {
    printf("\n\nError in fourier_init \n=>%s\n",fo.error_message);
    return _FAILURE_;
  }

  if (transfer_init(&pr,&ba,&th,&pt,&fo,&tr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (harmonic_init(&pr,&ba,&pt,&pm,&fo,&tr,&hr) == _FAILURE_) {
    printf("\n\nError in harmonic_init \n=>%s\n",hr.error_message);
    return _FAILURE_;
  }

  if (lensing_init(&pr,&pt,&hr,&fo,&le) == _FAILURE_) {
    printf("\n\nError in lensing_init \n=>%s\n",le.error_message);
    return _FAILURE_;
  }

  if (distortions_init(&pr,&ba,&th,&pt,&pm,&sd) == _FAILURE_) {
    printf("\n\nError in distortions_init \n=>%s\n",sd.error_message);
    return _FAILURE_;
  }

  if (output_init(&ba,&th,&pt,&pm,&tr,&hr,&fo,&le,&sd,&op) == _FAILURE_) {
    printf("\n\nError in output_init \n=>%s\n",op.error_message);
    return _FAILURE_;
  }

  /****** all calculations done, now free the structures ******/

  if (distortions_free(&sd) == _FAILURE_) {
    printf("\n\nError in distortions_free \n=>%s\n",sd.error_message);
    return _FAILURE_;
  }

  if (lensing_free(&le) == _FAILURE_) {
    printf("\n\nError in lensing_free \n=>%s\n",le.error_message);
    return _FAILURE_;
  }

  if (harmonic_free(&hr) == _FAILURE_) {
    printf("\n\nError in harmonic_free \n=>%s\n",hr.error_message);
    return _FAILURE_;
  }

  if (transfer_free(&tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (fourier_free(&fo) == _FAILURE_) {
    printf("\n\nError in fourier_free \n=>%s\n",fo.error_message);
    return _FAILURE_;
  }

  if (primordial_free(&pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (perturbations_free(&pt) == _FAILURE_) {
    printf("\n\nError in perturbations_free \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  /* ===== HS summary (stdout) ===== */

  if (ba.hs_model == _TRUE_ && ba.hs_R0 > 0.0) {
    double zstar = th.z_rec;     /* recombination redshift from thermo */
    double rs    = th.rs_rec;    /* comoving sound horizon at recombination */
    double DM    = 0.0;

    if (zstar <= 0.0) zstar = 1100.0;

    if (rs > 0.0 && background_D_M_of_z(&ba, zstar, &DM) == _SUCCESS_ && DM > 0.0) {
      double theta_star = rs / DM;                /* radians */
      if (theta_star > 0.0) {
        double ell_star = acos(-1.0) / theta_star;/* π/θ* without M_PI */
        printf("\n[HS] z* = %.3f, rs = %.6f Mpc, D_M(z*) = %.6f Mpc, "
               "theta* = %.9e rad, ell* ≈ %.2f\n",
               zstar, rs, DM, theta_star, ell_star);
      }
    }
  }
  /* =============================== */

    /* ===== HS R0 calibrator (stdout) ===== */
  if (ba.hs_model == _TRUE_ && ba.hs_R0 > 0.0) {
    /* Choose your target l_* (Planck ~ 301). You can tweak this number. */
    const double lstar_target = 301.0;

    /* Ingredients: sound horizon and recombination redshift */
    double zstar = th.z_rec;
    double rs    = th.rs_rec;

    if (zstar > 0.0 && rs > 0.0) {
      /* Compute chi = tau0 - tau(z*) in Mpc using CLASS background */
      double tau0 = 0.0, tauz = 0.0, chi = 0.0;
      if (background_tau_of_z(&ba, 0.0, &tau0) == _SUCCESS_ &&
          background_tau_of_z(&ba, zstar, &tauz) == _SUCCESS_) {

        chi = tau0 - tauz;

        /* Target D_M to hit l_*: D_M = (l_* / pi) * r_s */
        double DM_target = (lstar_target / acos(-1.0)) * rs;

        /* Solve R0 * sin(chi/R0) = DM_target for R0 via simple bisection */
        double Rlo = chi + 1.0;   /* just above chi to avoid sin(pi) issues */
        double Rhi = 1.0e7;       /* big upper bound ~ flat limit */
        double R0  = ba.hs_R0;
        for (int it = 0; it < 80; ++it) {
          double mid = 0.5*(Rlo + Rhi);
          double lhs = mid * sin(chi / mid);
          if (lhs < DM_target) Rlo = mid; else Rhi = mid;
          R0 = mid;
        }

        printf("[HS] calibrate: l*_target=%.2f, rs=%.6f Mpc, chi=%.6f Mpc -> "
               "D_M^target=%.6f Mpc, suggested hs_R0 ≈ %.3f Mpc\n",
               lstar_target, rs, chi, DM_target, R0);
      }
    }
  }
  /* ===================================== */

  /* ===== HS drag-epoch summary (stdout) ===== */
  if (ba.hs_model == _TRUE_ && ba.hs_R0 > 0.0) {
    double zd  = th.z_d;   /* baryon drag redshift */
    double rsd = th.rs_d;  /* comoving sound horizon at drag */
    if (zd > 0.0 && rsd > 0.0) {
      double DMd = 0.0;
      if (background_D_M_of_z(&ba, zd, &DMd) == _SUCCESS_ && DMd > 0.0) {
        double theta_d = rsd / DMd;
        if (theta_d > 0.0) {
          double ell_d = acos(-1.0) / theta_d;
          printf("[HS] z_d = %.3f, r_d = %.6f Mpc, D_M(z_d) = %.6f Mpc, "
                 "theta_d = %.9e rad, ell_d ≈ %.2f\n",
                 zd, rsd, DMd, theta_d, ell_d);
        }
      }
    }
  }
  /* ========================================= */  

  /* ===== HS CMB summary -> hs_cmb_summary.tsv (stdout-compatible) ===== */
  if (ba.hs_model == _TRUE_ && ba.hs_R0 > 0.0) {
    /* Recompute the numbers so this block stands alone */
    double zstar = th.z_rec;       /* recombination */
    double rs    = th.rs_rec;
    double zd    = th.z_d;         /* drag epoch */
    double rsd   = th.rs_d;

    double DMstar = 0.0, DMd = 0.0;
    double theta_star = 0.0, ell_star = 0.0;
    double theta_d    = 0.0, ell_d    = 0.0;

    /* Guard reasonable fallbacks */
    if (zstar <= 0.0) zstar = 1100.0;

    if (background_D_M_of_z(&ba, zstar, &DMstar) == _SUCCESS_ && DMstar > 0.0 && rs > 0.0) {
      theta_star = rs / DMstar;
      if (theta_star > 0.0) ell_star = acos(-1.0) / theta_star; /* π/θ* */
    }
    if (zd > 0.0 && rsd > 0.0 &&
        background_D_M_of_z(&ba, zd, &DMd) == _SUCCESS_ && DMd > 0.0) {
      theta_d = rsd / DMd;
      if (theta_d > 0.0) ell_d = acos(-1.0) / theta_d;
    }

    /* Write a simple TSV next to the executable (working directory) */
    FILE *f = fopen("hs_cmb_summary.tsv", "w");
    if (f != NULL) {
      fprintf(f, "z_star\trs_rec_Mpc\tD_M_zstar_Mpc\ttheta_star_rad\tell_star\t");
      fprintf(f, "z_d\tr_d_Mpc\tD_M_zd_Mpc\ttheta_d_rad\tell_d\n");
      fprintf(f, "%.6f\t%.9f\t%.9f\t%.12e\t%.6f\t", zstar, rs, DMstar, theta_star, ell_star);
      fprintf(f, "%.6f\t%.9f\t%.9f\t%.12e\t%.6f\n", zd, rsd, DMd, theta_d, ell_d);
      fclose(f);
      printf("[HS] wrote hs_cmb_summary.tsv\n");
    } else {
      printf("[HS] warning: could not open hs_cmb_summary.tsv for writing\n");
    }
  }
  /* ===================================================================== */

  /* ===== HS BAO summary -> hs_bao_summary.tsv ===== */
  if (ba.hs_model == _TRUE_ && ba.hs_R0 > 0.0) {
    /* Choose BAO z grid (edit as you like) */
    const double zlist[] = {0.106, 0.15, 0.32, 0.35, 0.57, 0.61, 0.8, 1.0, 1.5, 2.34};
    const int NZ = (int)(sizeof(zlist)/sizeof(zlist[0]));

    FILE *fbao = fopen("hs_bao_summary.tsv","w");
    if (fbao != NULL) {
      fprintf(fbao, "z\tD_M_Mpc\tD_A_Mpc\tH_1perMpc\tD_V_Mpc\n");

      for (int i=0; i<NZ; ++i) {
        double z = zlist[i];
        double DM = 0.0, DA = 0.0, H = 0.0, DV = 0.0;

        /* Our S^3 transverse comoving distance */
        if (background_D_M_of_z(&ba, z, &DM) != _SUCCESS_ || DM <= 0.0) {
          continue;
        }
        DA = DM / (1.0 + z);

        /* H(z) via conformal-time derivative: H = -1 / [(1+z) dτ/dz] */
        double dz = 1e-3;
        double tau1 = 0.0, tau2 = 0.0;
        if (background_tau_of_z(&ba, z,     &tau1) != _SUCCESS_ ||
            background_tau_of_z(&ba, z+dz,  &tau2) != _SUCCESS_) {
          continue;
        }
        double dtau_dz = (tau2 - tau1) / dz;
        if (dtau_dz == 0.0) continue;
        H = -1.0 / ((1.0 + z) * dtau_dz);
        if (H <= 0.0) continue;


        /* Eisenstein+Hu DV definition (CLASS units: c=1, so cz/H -> z/H) */
        double tmp = (1.0+z)*(1.0+z)*DA*DA * (z/H);
        if (tmp > 0.0) DV = pow(tmp, 1.0/3.0);

        fprintf(fbao, "%.6f\t%.9f\t%.9f\t%.12e\t%.9f\n", z, DM, DA, H, DV);
      }
      fclose(fbao);
      printf("[HS] wrote hs_bao_summary.tsv\n");
    } else {
      printf("[HS] warning: could not open hs_bao_summary.tsv for writing\n");
    }
  }
  /* ================================================= */

  /* ===== HS distance grid -> hs_distances.tsv ===== */
  if (ba.hs_model == _TRUE_ && ba.hs_R0 > 0.0) {
    const double zmin = 0.0, zmax = 3.0, dz = 0.01;
    FILE *fdg = fopen("hs_distances.tsv","w");
    if (fdg != NULL) {
      fprintf(fdg, "z\tchi_Mpc\tD_M_Mpc\tD_A_Mpc\tH_1perMpc\n");
      for (double z = zmin; z <= zmax + 1e-12; z += dz) {
        double tau0=0.0, tauz=0.0, chi=0.0;
        double DM=0.0, DA=0.0, H=0.0;

        /* comoving radial distance chi = tau0 - tau(z) */
        if (background_tau_of_z(&ba, 0.0, &tau0) != _SUCCESS_) continue;
        if (background_tau_of_z(&ba, z,   &tauz) != _SUCCESS_) continue;
        chi = tau0 - tauz;

        /* S^3 transverse comoving distance and D_A */
        if (background_D_M_of_z(&ba, z, &DM) != _SUCCESS_ || DM <= 0.0) continue;
        DA = DM / (1.0 + z);

        /* H(z) via finite difference on tau(z): H = -1 / [(1+z) dτ/dz] */
        double dzH = 1e-3, t1=0.0, t2=0.0;
        if (background_tau_of_z(&ba, z,     &t1) != _SUCCESS_) continue;
        if (background_tau_of_z(&ba, z+dzH, &t2) != _SUCCESS_) continue;
        double dtau_dz = (t2 - t1)/dzH;
        if (dtau_dz == 0.0) continue;
        H = -1.0 / ((1.0 + z) * dtau_dz);
        if (H <= 0.0) continue;

        fprintf(fdg, "%.4f\t%.9f\t%.9f\t%.9f\t%.12e\n", z, chi, DM, DA, H);
      }
      fclose(fdg);
      printf("[HS] wrote hs_distances.tsv\n");
    } else {
      printf("[HS] warning: could not open hs_distances.tsv for writing\n");
    }
  }
  /* ================================================ */

if (thermodynamics_free(&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;

}
