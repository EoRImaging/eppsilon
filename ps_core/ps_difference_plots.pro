pro ps_difference_plots, info_files, cube_types, pols, $
    plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    note = note, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, plot_wedge_line = plot_wedge_line, $
    quiet = quiet, png = png, eps = eps
    
  if n_elements(info_files) gt 2 then message, 'Only 1 or 2 info_files can be used'
  if n_elements(cube_types) gt 2 then message, 'Only 1 or 2 info_files can be used'
  if n_elements(pols) gt 2 then message, 'Only 1 or 2 info_files can be used'
  
  if max([n_elements(info_files), n_elements(cube_types), n_elements(pols)]) eq 1 then $
    message, 'at least one of info_files, cube_types, pols must be a 2 element vector'
    
  if n_elements(cube_types) eq 1 then cube_types = [cube_types, cube_types]
  if n_elements(pols) eq 1 then pols = [pols, pols]
  type_pol_str = cube_types + '_' + pols
  
  ;; default to including baseline axis & delay axis
  if n_elements(baseline_axis) eq 0 then baseline_axis = 1
  if n_elements(delay_axis) eq 0 then delay_axis = 1
  
  ;; default to blackman-harris spectral window
  if not keyword_set(no_spec_window) then begin
    if n_elements(spec_window_type) eq 0 then spec_window_type = 'Blackman-Harris'
  endif else undefine, spec_window_type
  
  ;; default to plot wedge line
  if n_elements(plot_wedge_line) eq 0 then plot_wedge_line=1
  
  ;; default to hinv
  if n_elements(hinv) eq 0 then hinv = 1
  
  file_struct_arr1 = fhd_file_setup(info_files[0], pol_inc, spec_window_type = spec_window_type)
  
  pol_ind1 = where(pol_inc eq pols[0], count_pol)
  if count_pol eq 0 then message, 'pol ' + pols[0] + ' not included in info_file: ' + info_files[0]
  
  cube_ind1 = where(file_struct_arr1.type_pol_str eq type_pol_str[0], count_typepol)
  if count_typepol eq 0 then message, 'type ' + cube_types[0] + ' not included in info_file: ' + info_files[0]
  
  power2dfile_1 = file_struct_arr1[cube_ind1].savefile_froot + file_struct_arr1[cube_ind1].savefilebase + file_struct_arr1[cube_ind1].power_tag + '_2dkpower.idlsave'
  
  power_file1 = file_struct_arr1[cube_ind1].power_savefile
  if file_test(power_file1) eq 0 then message, 'No power file for ' + type_pol_str[0] + ' and info_file: ' + info_files[0]
  
  if n_elements(info_files) eq 2 then file_struct_arr2 = fhd_file_setup(info_files[1], pol_inc, spec_window_type = spec_window_type) $
  else file_struct_arr2 = file_struct_arr1
  
  pol_ind2 = where(pol_inc eq pols[1], count_pol)
  if count_pol eq 0 then message, 'pol ' + pols[1] + ' not included in info_file: ' + info_files[n_elements(info_files)-1]
  
  cube_ind2 = where(file_struct_arr2.type_pol_str eq type_pol_str[1], count_typepol)
  if count_typepol eq 0 then message, 'type ' + cube_types[1] + ' not included in info_file: ' + info_files[n_elements(info_files)-1]
  
  power2dfile_2 = file_struct_arr2[cube_ind2].savefile_froot + file_struct_arr2[cube_ind2].savefilebase + file_struct_arr2[cube_ind2].power_tag + '_2dkpower.idlsave'
  
  power_file2 = file_struct_arr2[cube_ind2].power_savefile
  if file_test(power_file2) eq 0 then message, 'No power file for ' + type_pol_str[1] + ' and info_file: ' + info_files[n_elements(info_files)-1]
  
  kx1 = getvar_savefile(power_file1, 'kx_mpc')
  kx2 = getvar_savefile(power_file2, 'kx_mpc')
  if n_elements(kx1) ne n_elements(kx2) or total(abs(kx1 - kx2)) ne 0 then message, 'kx_mpc does not match between cubes'
  ky1 = getvar_savefile(power_file1, 'ky_mpc')
  ky2 = getvar_savefile(power_file2, 'ky_mpc')
  if n_elements(ky1) ne n_elements(ky2) or total(abs(ky1 - ky2)) ne 0 then message, 'ky_mpc does not match between cubes'
  kz1 = getvar_savefile(power_file1, 'kz_mpc')
  kz2 = getvar_savefile(power_file2, 'kz_mpc')
  if n_elements(kz1) ne n_elements(kz2) or total(abs(kz1 - kz2)) ne 0 then message, 'kz_mpc does not match between cubes'
  
  ;; if kx, ky, and kz are all the same these should be too
  kperp_lambda_conv = getvar_savefile(power_file1, 'kperp_lambda_conv')
  delay_params = getvar_savefile(power_file1, 'delay_params')
  hubble_param = getvar_savefile(power_file1, 'hubble_param')
  
  
  power1 = getvar_savefile(power_file1, 'power_3d')
  power2 = getvar_savefile(power_file2, 'power_3d')
  power_diff = power1 - power2
  if max(abs(power_diff)) eq 0 then message, 'The cubes are identical -- power difference is zero everywhere'
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
  
  
  power_rebin = kspace_rebinning_2d(power_diff, kx1, ky1, kz1, kperp_edges_mpc, kpar_edges_mpc, log_kpar = log_kpar, $
    log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, weights = weight_diff, $
    binned_weights = binned_weights)
    
  power = power_rebin
  kperp_edges = kperp_edges_mpc
  kpar_edges = kpar_edges_mpc
  weights = binned_weights
  
  wh_good_kperp = where(total(weights, 2) gt 0, count)
  if count eq 0 then stop
  kperp_plot_range = [min(kperp_edges[wh_good_kperp]), max(kperp_edges[wh_good_kperp+1])]
  
  
  ;; if not save_path specified, save in folder of first info file
  if n_elements(save_path) eq 0 then save_path = file_dirname(info_files[0], /mark_directory)
  if n_elements(plot_path) eq 0 then plot_path = save_path
  
  if n_elements(info_files) eq 2 then begin
    fileparts_1 = strsplit(file_struct_arr1[0].general_filebase, '_', /extract)
    fileparts_2 = strsplit(file_struct_arr2[0].general_filebase, '_', /extract)
    match_test = strcmp(fileparts_1, fileparts_2)
    wh_diff = where(match_test eq 0, count_diff, complement = wh_same, ncomplement = count_same)
  endif
  
  if n_elements(savefilebase) eq 0 then begin
    if n_elements(info_files) eq 1 then savefilebase = file_struct_arr1[0].general_filebase + $
      '_' + type_pol_str[0] + '_minus_' + type_pol_str[1] else begin
      if count_diff eq 0 then savefilebase = file_struct_arr1[0].general_filebase + '_diff' else begin
        if count_same gt 0 then savefilebase = strjoin(fileparts_1[wh_same], '_') + '__' + $
          strjoin(fileparts_1[wh_diff]) + '_minus_' + strjoin(fileparts_2[wh_diff]) $
        else $
          savefilebase = file_struct_arr1[0].general_filebase + '_' + type_pol_str[0] + $
          '_minus_' + file_struct_arr2[0].general_filebase + '_' + type_pol_str[1]
      endelse
    endelse
  endif
  if n_elements(plot_filebase) eq 0 then plot_filebase = savefilebase
  
  savefile = save_path + savefilebase + '.idlsave'
  if keyword_set(png) or keyword_set(eps) then pub = 1 else pub = 0
  if pub then begin
    plotfile = plot_path + plot_filebase
    
    plotfile_2d_1 = plot_path + file_struct_arr1[cube_ind1].savefilebase + '_2dkpower'
    plotfile_2d_2 = plot_path + file_struct_arr2[cube_ind2].savefilebase + '_2dkpower'
  endif
  
  title = type_pol_str[0] + '-' + type_pol_str[1]
  
  save, file = savefile, power, weights, kperp_edges, kpar_edges, kperp_bin, kpar_bin, $
    kperp_lambda_conv, delay_params, hubble_param
    
  kperp_plot_range = [5./kperp_lambda_conv, min([file_struct_arr1.kspan/2.,file_struct_arr1.max_baseline_lambda])/kperp_lambda_conv]
  
  if keyword_set(hinv) then kperp_plot_range = kperp_plot_range / hubble_param
  
  if keyword_set(plot_wedge_line) then begin
    z0_freq = 1420.40 ;; MHz
    redshifts1 = z0_freq/file_struct_arr1[0].frequencies - 1
    redshifts2 = z0_freq/file_struct_arr2[0].frequencies - 1
    mean_redshift = mean([redshifts1, redshifts2])
    
    cosmology_measures, mean_redshift, wedge_factor = wedge_factor
    ;; assume 20 degrees from pointing center to first null
    source_dist = 20d * !dpi / 180d
    fov_amp = wedge_factor * source_dist
    
    ;; calculate angular distance to horizon
    max_theta = max([file_struct_arr1[0].max_theta, file_struct_arr2[0].max_theta])
    horizon_amp = wedge_factor * ((max_theta+90d) * !dpi / 180d)
    
    wedge_amp = [fov_amp, horizon_amp]
  endif else wedge_amp = 0d
  
  if keyword_set(pub) then font = 1 else font = -1
  
  if not keyword_set(quiet) then begin
    kpower_2d_plots, power2dfile_1, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
      data_range = data_range, png = png, eps = eps, plotfile = plotfile_2d_1, full_title=title, window_num = 1, $
      kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, baseline_axis = baseline_axis, delay_axis = delay_axis, $
      wedge_amp = wedge_amp, plot_wedge_line = plot_wedge_line, hinv = hinv
      
    kpower_2d_plots, power2dfile_2, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
      data_range = data_range, png = png, eps = eps, plotfile = plotfile_2d_2, full_title=title, window_num = 2, $
      kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, baseline_axis = baseline_axis, delay_axis = delay_axis, $
      wedge_amp = wedge_amp, plot_wedge_line = plot_wedge_line, hinv = hinv
      
    kpower_2d_plots, savefile, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, note = note, $
      data_range = diff_range, png = png, eps = eps, plotfile = plotfile, full_title=title, window_num = 3, color_profile = 'sym_log', $
      kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, baseline_axis = baseline_axis, delay_axis = delay_axis, $
      wedge_amp = wedge_amp, plot_wedge_line = plot_wedge_line, hinv = hinv
  endif
  
end