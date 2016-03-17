pro ps_ratio_plots, folder_names, obs_info, cube_types, pols, all_pol_diff_ratio = all_pol_diff_ratio, $
    freq_ch_range = freq_ch_range, plot_path = plot_path, plot_filebase = plot_filebase, $
    note = note, spec_window_types = spec_window_types, ave_removal = ave_removal, data_range = data_range, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, $
    diff_ratio = diff_ratio, diff_range = diff_range, diff_min_abs = diff_min_abs, $
    wt_cutoffs = wt_cutoffs, wt_measures = wt_measures, invert_colorbar = invert_colorbar, $
    plot_wedge_line = plot_wedge_line, quiet = quiet, png = png, eps = eps, pdf = pdf, $
    window_num = window_num
    
  if keyword_set(png) or keyword_set(eps) or keyword_set(pdf) then pub = 1 else pub = 0
  if pub then begin
    if keyword_set(png) and keyword_set(eps) and keyword_set(pdf) then begin
      print, 'only one of eps, pdf and png can be set, using png'
      eps = 0
    endif
    
    if keyword_set(png) then begin
      plot_exten = '.png'
      delete_ps = 1
    endif else if keyword_set(pdf) then begin
      plot_exten = '.pdf'
      delete_ps = 1
    endif else if keyword_set(eps) then begin
      plot_exten = '.eps'
      delete_ps = 0
    endif
  endif
  
  fadd_2dbin = ''
  ;;if keyword_set(fill_holes) then fadd_2dbin = fadd_2dbin + '_nohole'
  if keyword_set(no_kzero) then fadd_2dbin = fadd_2dbin + '_nok0'
  if keyword_set(log_kpar) then fadd_2dbin = fadd_2dbin + '_logkpar'
  if keyword_set(log_kperp) then fadd_2dbin = fadd_2dbin + '_logkperp'
  
  
  compare_plot_prep, folder_names, obs_info, cube_types, pols, 'ratio', compare_files, $
    plot_slices = plot_slices, slice_type = slice_type, fadd_2dbin = fadd_2dbin, $
    spec_window_types = spec_window_types, freq_ch_range = freq_ch_range, ave_removal = ave_removal, $
    baseline_axis = baseline_axis, delay_axis = delay_axis, plot_wedge_line = plot_wedge_line, hinv = hinv, $
    plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    axis_type_1d = axis_type_1d, note = note, wt_cutoffs = wt_cutoffs, wt_measures = wt_measures, $
    pub = pub, plot_exten = plot_exten, full_compare = full_compare
    
  test_save1 = file_test(compare_files.input_savefile1) *  (1 - file_test(compare_files.input_savefile1, /zero_length))
  test_save2 = file_test(compare_files.input_savefile2) *  (1 - file_test(compare_files.input_savefile2, /zero_length))
  
  if min(test_save1) eq 0 then message, '2D savefile not found: ' + strjoin(compare_files.input_savefile1[where(test_save1 eq 0)], ', ')
  if min(test_save2) eq 0 then message, '2D savefile not found: ' + strjoin(compare_files.input_savefile2[where(test_save2 eq 0)], ', ')
  
  
  if keyword_set(diff_ratio) or keyword_set(all_pol_diff_ratio) then begin
  
    for i=0, n_sets/2-1 do begin
stop
      kperp_edges = getvar_savefile(savefiles_2d[2*i,0], 'kperp_edges')
      kperp_edges2 = getvar_savefile(savefiles_2d[2*i+1,0], 'kperp_edges')
      if n_elements(kperp_edges) ne n_elements(kperp_edges2) or max(abs(kperp_edges - kperp_edges2)) gt 1.05e-3 then $
        message, 'kperp_edges does not match in savefiles'
      kpar_edges = getvar_savefile(savefiles_2d[2*i,0], 'kpar_edges')
      kpar_edges2 = getvar_savefile(savefiles_2d[2*i+1,0], 'kpar_edges')
      if n_elements(kpar_edges) ne n_elements(kpar_edges2) or max(abs(kpar_edges - kpar_edges2)) gt 1.05e-3 then $
        message, 'kpar_edges does not match in savefiles'
        
      kperp_bin = getvar_savefile(savefiles_2d[2*i,0], 'kperp_bin')
      kperp_bin2 = getvar_savefile(savefiles_2d[2*i+1,0], 'kperp_bin')
      if abs(kperp_bin - kperp_bin2) gt 1.05e-3 then $
        message, 'kperp_bin does not match in savefiles'
      kpar_bin = getvar_savefile(savefiles_2d[2*i,0], 'kpar_bin')
      kpar_bin2 = getvar_savefile(savefiles_2d[2*i+1,0], 'kpar_bin')
      if abs(kpar_bin - kpar_bin2) gt 1.05e-3 then $
        message, 'kpar_edges does not match in savefiles'
        
      kperp_lambda_conv = getvar_savefile(savefiles_2d[2*i,0], 'kperp_lambda_conv')
      kperp_lambda_conv2 = getvar_savefile(savefiles_2d[2*i+1,0], 'kperp_lambda_conv')
      if abs(kperp_lambda_conv - kperp_lambda_conv2)/kperp_lambda_conv gt 1.05e-3  then message, 'kperp_lambda_conv do not match in savefiles'
      delay_params = getvar_savefile(savefiles_2d[2*i,0], 'delay_params')
      delay_params2 = getvar_savefile(savefiles_2d[2*i+1,0], 'delay_params')
      if max(abs(delay_params - delay_params2)) gt 1.05e-3  then message, 'delay_params do not match in savefiles'
      hubble_param = getvar_savefile(savefiles_2d[2*i,0], 'hubble_param')
      if total(abs(hubble_param - getvar_savefile(savefiles_2d[2*i+1,0], 'hubble_param'))) ne 0 then message, 'hubble_param do not match in savefiles'
      
      power1 = getvar_savefile(savefiles_2d[2*i,0], 'power')
      power2 = getvar_savefile(savefiles_2d[2*i,1], 'power')
      
      power_ratio1 = power1 / power2
      wh0 = where(power2 eq 0, count0)
      if count0 gt 0 then power_ratio1[wh0] = 0
      
      power1 = getvar_savefile(savefiles_2d[2*i+1,0], 'power')
      power2 = getvar_savefile(savefiles_2d[2*i+1,1], 'power')
      
      power_ratio2 = power1 / power2
      wh0 = where(power2 eq 0, count0)
      if count0 gt 0 then power_ratio2[wh0] = 0
      
      power_diff_ratio = power_ratio1-power_ratio2
      
      if i eq 0 then begin
        kperp_plot_range = [5./kperp_lambda_conv, min([file_struct_arr1.kspan/2.,file_struct_arr1.max_baseline_lambda])/kperp_lambda_conv]
        if keyword_set(hinv) then kperp_plot_range = kperp_plot_range / hubble_param
        
        
        ncol=3
        nrow=n_sets/2
        start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
      endif else pos_use = positions[*,3*i]
      
      kpower_2d_plots, power = power_ratio1, multi_pos = pos_use, start_multi_params = start_multi_params, $
        kperp_edges = kperp_edges, kpar_edges = kpar_edges, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
        kperp_lambda_conv = kperp_lambda_conv, delay_params = delay_params, hubble_param = hubble_param, $
        png = png, eps = eps, pdf = pdf, plotfile = plotfile_2d, window_num = window_num, $
        kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
        data_range = data_range, full_title = titles[i,0], $
        plot_wedge_line = plot_wedge_line, hinv = hinv, /pwr_ratio, $
        wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, $
        kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
        
      if i eq 0 then begin
        positions = pos_use
        undefine, start_multi_params
      endif
      pos_use = positions[*,3*i+1]
      
      kpower_2d_plots, power = power_ratio2, multi_pos = pos_use, start_multi_params = start_multi_params, $
        kperp_edges = kperp_edges, kpar_edges = kpar_edges, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
        kperp_lambda_conv = kperp_lambda_conv, delay_params = delay_params, hubble_param = hubble_param, $
        png = png, eps = eps, pdf = pdf, plotfile = plotfile_2d, window_num = window_num, $
        kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
        data_range = data_range, full_title = titles[i,1], $
        plot_wedge_line = plot_wedge_line, hinv = hinv, /pwr_ratio, $
        wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, $
        kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
        
      pos_use = positions[*,3*i+2]
      
      kpower_2d_plots, power = power_diff_ratio, multi_pos = pos_use, start_multi_params = start_multi_params, $
        kperp_edges = kperp_edges, kpar_edges = kpar_edges, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
        kperp_lambda_conv = kperp_lambda_conv, delay_params = delay_params, hubble_param = hubble_param, $
        png = png, eps = eps, pdf = pdf, plotfile = plotfile_2d, window_num = window_num, color_profile = 'sym_log', $
        kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
        data_range = diff_range, data_min_abs = diff_min_abs, full_title = titles[i,2], note = note + ', ' + kperp_density_names, $
        plot_wedge_line = plot_wedge_line, hinv = hinv, /pwr_ratio, $
        wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, $
        kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
        
    endfor
    
    if keyword_set(pub) then begin
      cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density=600
    endif
    
    
  endif else begin
  
    kperp_edges = getvar_savefile(savefiles_2d[0], 'kperp_edges')
    kperp_edges2 = getvar_savefile(savefiles_2d[1], 'kperp_edges')
    if n_elements(kperp_edges) ne n_elements(kperp_edges2) or max(abs(kperp_edges - kperp_edges2)) gt 1.05e-3 then $
      message, 'kperp_edges does not match in savefiles'
    kpar_edges = getvar_savefile(savefiles_2d[0], 'kpar_edges')
    kpar_edges2 = getvar_savefile(savefiles_2d[1], 'kpar_edges')
    if n_elements(kpar_edges) ne n_elements(kpar_edges2) or max(abs(kpar_edges - kpar_edges2)) gt 1.05e-3 then $
      message, 'kpar_edges does not match in savefiles'
      
    kperp_bin = getvar_savefile(savefiles_2d[0], 'kperp_bin')
    kperp_bin2 = getvar_savefile(savefiles_2d[1], 'kperp_bin')
    if abs(kperp_bin - kperp_bin2) gt 1.05e-3 then $
      message, 'kperp_bin does not match in savefiles'
    kpar_bin = getvar_savefile(savefiles_2d[0], 'kpar_bin')
    kpar_bin2 = getvar_savefile(savefiles_2d[1], 'kpar_bin')
    if abs(kpar_bin - kpar_bin2) gt 1.05e-3 then $
      message, 'kpar_edges does not match in savefiles'
      
    kperp_lambda_conv = getvar_savefile(savefiles_2d[0], 'kperp_lambda_conv')
    kperp_lambda_conv2 = getvar_savefile(savefiles_2d[1], 'kperp_lambda_conv')
    if abs(kperp_lambda_conv - kperp_lambda_conv2)/kperp_lambda_conv gt 1.05e-3  then message, 'kperp_lambda_conv do not match in savefiles'
    delay_params = getvar_savefile(savefiles_2d[0], 'delay_params')
    delay_params2 = getvar_savefile(savefiles_2d[1], 'delay_params')
    if max(abs(delay_params - delay_params2)) gt 1.05e-3  then message, 'delay_params do not match in savefiles'
    hubble_param = getvar_savefile(savefiles_2d[0], 'hubble_param')
    if total(abs(hubble_param - getvar_savefile(savefiles_2d[1], 'hubble_param'))) ne 0 then message, 'hubble_param do not match in savefiles'
    
    power1 = getvar_savefile(savefiles_2d[0], 'power')
    power2 = getvar_savefile(savefiles_2d[1], 'power')
    
    power_ratio = power1 / power2
    wh0 = where(power2 eq 0, count0)
    if count0 gt 0 then power_ratio[wh0] = 0
    
    kperp_plot_range = [5./kperp_lambda_conv, min([file_struct_arr1.kspan/2.,file_struct_arr1.max_baseline_lambda])/kperp_lambda_conv]
    if keyword_set(hinv) then kperp_plot_range = kperp_plot_range / hubble_param
    
    
    kpower_2d_plots, power = power_ratio, kperp_edges = kperp_edges, kpar_edges = kpar_edges, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
      kperp_lambda_conv = kperp_lambda_conv, delay_params = delay_params, hubble_param = hubble_param, $
      png = png, eps = eps, pdf = pdf, plotfile = plotfile_2d, window_num = window_num, $
      kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
      data_range = data_range, title_prefix = title, note = note + ', ' + kperp_density_names, $
      plot_wedge_line = plot_wedge_line, hinv = hinv, /pwr_ratio, $
      wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, $
      kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
      
  endelse
  
end