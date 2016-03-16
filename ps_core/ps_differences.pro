pro ps_differences, power_file1, power_file2, refresh = refresh, $
    savefile_3d = savefile_3d, savefile_2d = savefile_2d, savefiles_1d = savefiles_1d, $
    diff_ratio = diff_ratio, wedge_amp = wedge_amp, wt_cutoffs = wt_cutoffs, wt_measures = wt_measures
    
  test_save = file_test(savefile_3d) *  (1 - file_test(savefile_3d, /zero_length))
  test_save_2d = file_test(savefile_2d) *  (1 - file_test(savefile_2d, /zero_length))
  test_save_1d = file_test(reform(savefiles_1d)) *  (1 - file_test(reform(savefiles_1d) , /zero_length))
  
  if test_save eq 0 or keyword_set(refresh) then begin
  
    if file_test(power_file1) eq 0 then message, 'file not found: ' + power_file1
    if file_test(power_file2) eq 0 then message, 'file not found: ' + power_file2
    
    kx_mpc = getvar_savefile(power_file1, 'kx_mpc')
    kx_mpc2 = getvar_savefile(power_file2, 'kx_mpc')
    if n_elements(kx_mpc) ne n_elements(kx_mpc2) or max(abs(kx_mpc - kx_mpc2)) gt 1.05e-3 then message, 'kx_mpc does not match between cubes'
    ky_mpc = getvar_savefile(power_file1, 'ky_mpc')
    ky_mpc2 = getvar_savefile(power_file2, 'ky_mpc')
    if n_elements(ky_mpc) ne n_elements(ky_mpc2) or max(abs(ky_mpc - ky_mpc2)) gt 1.05e-3 then message, 'ky_mpc does not match between cubes'
    kz_mpc = getvar_savefile(power_file1, 'kz_mpc')
    kz_mpc2 = getvar_savefile(power_file2, 'kz_mpc')
    if n_elements(kz_mpc) ne n_elements(kz_mpc2) or max(abs(kz_mpc - kz_mpc)) gt 1.05e-3 then message, 'kz_mpc does not match between cubes'
    
    ;; if kx, ky, and kz are all the same these should be too
    kperp_lambda_conv = getvar_savefile(power_file1, 'kperp_lambda_conv')
    delay_params = getvar_savefile(power_file1, 'delay_params')
    hubble_param = getvar_savefile(power_file1, 'hubble_param')
    
    power1 = getvar_savefile(power_file1, 'power_3d')
    power2 = getvar_savefile(power_file2, 'power_3d')
    power_diff = power1 - power2
    if max(abs(power_diff)) eq 0 then begin
      print, 'The cubes are identical -- power difference is zero everywhere'
    ;continue
    endif
    undefine, power1, power2
    
    weights1 = getvar_savefile(power_file1, 'weights_3d')
    weights2 = getvar_savefile(power_file2, 'weights_3d')
    
    ;; variance_3d = 1/weights_3d
    var1 = 1./weights1
    wh_wt1_0 = where(weights1 eq 0, count_wt1_0)
    if count_wt1_0 gt 0 then var1[wh_wt1_0] = 0
    var2 = 1./weights2
    wh_wt2_0 = where(weights2 eq 0, count_wt2_0)
    if count_wt2_0 gt 0 then var2[wh_wt2_0] = 0
    undefine, weights1, weights2
    
    var_diff = var1 + var2
    weight_diff = 1/var_diff
    if count_wt1_0 gt 0 then weight_diff[wh_wt1_0] = 0
    if count_wt2_0 gt 0 then weight_diff[wh_wt2_0] = 0
    undefine, var1, var2, var_diff
    
    
    wt_meas_ave1 = getvar_savefile(power_file1, 'wt_meas_ave')
    wt_meas_ave2 = getvar_savefile(power_file2, 'wt_meas_ave')
    wt_meas_ave = wt_meas_ave1 < wt_meas_ave2
    
    wt_meas_min1 = getvar_savefile(power_file1, 'wt_meas_min')
    wt_meas_min2 = getvar_savefile(power_file2, 'wt_meas_min')
    wt_meas_min = wt_meas_min1 < wt_meas_min2
    
    save, file = savefile_3d, power_diff, weight_diff, wt_meas_ave, wt_meas_min, $
      kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param
      
  endif else if min([test_save_1d, test_save_2d]) eq 0 then restore, savefile
  
  if test_save_2d eq 0 or keyword_set(refresh) then begin
  
    if wt_cutoffs gt 0 then begin
      case wt_measures of
        'ave': wt_meas_use = wt_meas_ave
        'min': wt_meas_use = wt_meas_min
      endcase
    endif
    
    power_rebin = kspace_rebinning_2d(power_diff, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, log_kpar = log_kpar, $
      log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, weights = weight_diff, $
      binned_weights = binned_weights, $
      kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoffs)
      
    power = power_rebin
    kperp_edges = kperp_edges_mpc
    kpar_edges = kpar_edges_mpc
    weights = binned_weights
    
    if keyword_set(diff_ratio) then begin
      if n_elements(power1) eq 0 then power1 = getvar_savefile(power_file1, 'power_3d')
      if n_elements(weights1) eq 0 then weights1 = getvar_savefile(power_file1, 'weights_3d')
      
      power_rebin_1 = kspace_rebinning_2d(power1, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, log_kpar = log_kpar, $
        log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, weights = weights1, $
        binned_weights = binned_weights1, $
        kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoffs)
        
      if n_elements(power2) eq 0 then power2 = getvar_savefile(power_file2, 'power_3d')
      if n_elements(weights2) eq 0 then weights2 = getvar_savefile(power_file2, 'weights_3d')
      
      power_rebin_2 = kspace_rebinning_2d(power2, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, log_kpar = log_kpar, $
        log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, weights = weights2, $
        binned_weights = binned_weights2, $
        kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoffs)
        
      power_denom = power_rebin_2
      weights_denom = binned_weights2
      
      power = power/power_denom
      weights = weights/weights_denom
    endif
    
    save, file = savefile_2d, power, weights, kperp_edges, kpar_edges, kperp_bin, kpar_bin, $
      kperp_lambda_conv, delay_params, hubble_param
  endif
  
  if min(test_save_1d) eq 0 or keyword_set(refresh) then begin
  
    for wedge_i=0, n_elements(wedge_amp) do begin
      if wedge_i gt 0 then wedge_amp_use = wedge_amp[wedge_i-1]
      
      if wt_cutoffs gt 0 then begin
        case wt_measures of
          'ave': wt_meas_use = wt_meas_ave
          'min': wt_meas_use = wt_meas_min
        endcase
      endif
      
      power_rebin = kspace_rebinning_1d(power_diff, kx_mpc, ky_mpc, kz_mpc, k_edges_mpc, k_bin = k_bin, log_k = log_k, $
        weights = weight_diff, binned_weights = binned_weights, wedge_amp = wedge_amp_use, $
        kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoffs)
        
      power = power_rebin
      k_edges = k_edges_mpc
      weights = binned_weights
      
      save, file = savefiles_1d[wedge_i], power, weights, k_edges, k_bin, kperp_lambda_conv, delay_params, hubble_param
    endfor
  endif
  
end