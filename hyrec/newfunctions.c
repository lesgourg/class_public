
int thermodynamics_annihilation_coefficients_interpolate(REC_COSMOPARAMS *param,
                                                          double xe,
                                                          double &chi_heat,
                                                          double &chi_lya,
                                                          double &chi_ionH,
                                                          double &chi_ionHe,
                                                          double &chi_lowE
      int thermodynamics_annihilation_coefficients_interpolate(param,
                                                               xe,
                                                               chi_heat,
                                                               chi_lya,
                                                               chi_ionH,
                                                               chi_ionHe,
                                                               chi_lowE
 );
 int thermodynamics_annihilation_coefficients_interpolate(REC_COSMOPARAMS *param,
 double xe,
 double &chi_heat,
 double &chi_lya,
 double &chi_ionH,
 double &chi_ionHe,
 double &chi_lowE
 ) {

 int last_index;
 ErrorMsg error_message;
 array_interpolate_spline(param->annihil_coef_xe,
 param->annihil_coef_num_lines,
 param->annihil_coef_heat,
 param->annihil_coef_dd_heat,
 1,
 xe,
 &last_index,
 &(param->chi_heat),
 1,
 error_message);

 array_interpolate_spline(param->annihil_coef_xe,
 param->annihil_coef_num_lines,
 param->annihil_coef_lya,
 param->annihil_coef_dd_lya,
 1,
 xe,
 &last_index,
 &(param->chi_lya),
 1,
 error_message);

 array_interpolate_spline(param->annihil_coef_xe,
 param->annihil_coef_num_lines,
 param->annihil_coef_ionH,
 param->annihil_coef_dd_ionH,
 1,
 xe,
 &last_index,
 &(param->chi_ionH),
 1,
 error_message);
 array_interpolate_spline(param->annihil_coef_xe,
 param->annihil_coef_num_lines,
 param->annihil_coef_ionHe,
 param->annihil_coef_dd_ionHe,
 1,
 xe,
 &last_index,
 &(param->chi_ionHe),
 1,
 error_message);
 array_interpolate_spline(param->annihil_coef_xe,
 param->annihil_coef_num_lines,
 param->annihil_coef_lowE,
 param->annihil_coef_dd_lowE,
 1,
 xe,
 &last_index,
 &(param->chi_lowE),
 1,
 error_message);

 return _SUCCESS_;

 }
