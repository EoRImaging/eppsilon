pro ps_differences, power_file1, power_file2, refresh = refresh, $
    savefile_3d = savefile_3d, savefile_2d = savefile_2d, savefiles_1d = savefiles_1d, $
    wedge_amp = wedge_amp, wt_cutoffs = wt_cutoffs, wt_measures = wt_measures, $
    plot_options = plot_options

  test_save = file_valid(savefile_3d)
  test_save_2d = file_valid(savefile_2d)
  test_save_1d = file_valid(reform(savefiles_1d))

  if test_save eq 0 or keyword_set(refresh) then begin

    if file_valid(power_file1) eq 0 then begin
      message, 'file not found: ' + power_file1
    endif
    if file_valid(power_file2) eq 0 then begin
      message, 'file not found: ' + power_file2
    endif

    kx_mpc = getvar_savefile(power_file1, 'kx_mpc')
    kx_mpc2 = getvar_savefile(power_file2, 'kx_mpc')
    if n_elements(kx_mpc) ne n_elements(kx_mpc2) or max(abs(kx_mpc - kx_mpc2)) gt 1.05e-3 then begin
      message, 'kx_mpc does not match between cubes'
    endif
    ky_mpc = getvar_savefile(power_file1, 'ky_mpc')
    ky_mpc2 = getvar_savefile(power_file2, 'ky_mpc')
    if n_elements(ky_mpc) ne n_elements(ky_mpc2) or max(abs(ky_mpc - ky_mpc2)) gt 1.05e-3 then begin
      message, 'ky_mpc does not match between cubes'
    endif
    kz_mpc = getvar_savefile(power_file1, 'kz_mpc')
    kz_mpc2 = getvar_savefile(power_file2, 'kz_mpc')
    if n_elements(kz_mpc) ne n_elements(kz_mpc2) or max(abs(kz_mpc - kz_mpc)) gt 1.05e-3 then begin
      message, 'kz_mpc does not match between cubes'
    endif

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

  endif else if min([test_save_1d, test_save_2d]) eq 0 then restore, savefile_3d

  if test_save_2d eq 0 or keyword_set(refresh) then begin

    if wt_cutoffs gt 0 then begin
      case wt_measures of
        'ave': wt_meas_use = wt_meas_ave
        'min': wt_meas_use = wt_meas_min
      endcase
    endif
    binning_2d_options = create_binning_2d_options()

    power_rebin = kspace_rebinning_2d(power_diff, kx_mpc, ky_mpc, kz_mpc, $
      kperp_edges_mpc, kpar_edges_mpc, $
      weights = weight_diff, binned_weights = binned_weights, $
      kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoffs, $
      binning_2d_options = binning_2d_options)

    power = power_rebin
    kperp_edges = kperp_edges_mpc
    kpar_edges = kpar_edges_mpc
    weights = binned_weights
    kperp_bin = binning_2d_options.kperp_bin
    kpar_bin = binning_2d_options.kpar_bin

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
      binning_1d_options = create_binning_1d_options()

      power_rebin = kspace_rebinning_1d(power_diff, kx_mpc, ky_mpc, kz_mpc, $
        k_edges_mpc, weights = weight_diff, binned_weights = binned_weights, $
        wedge_amp = wedge_amp_use, hubble_param = hubble_param, kperp_lambda_conv, $
        kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoffs, $
        binning_1d_options = binning_1d_options, plot_options = plot_options)

      power = power_rebin
      k_edges = k_edges_mpc
      weights = binned_weights
      k_bin = binning_1d_options.k_bin

      save, file = savefiles_1d[wedge_i], power, weights, k_edges, k_bin, $
        kperp_lambda_conv, delay_params, hubble_param
    endfor
  endif

end
