pro ps_ratio_wrapper, folder_names, obs_names_in, exact_obsnames = exact_obsnames, cube_types = cube_types, $
    pols = pols, all_pol_diff_ratio = all_pol_diff_ratio, freq_ch_range = freq_ch_range, $
    spec_window_types = spec_window_types, delta_uv_lambda = delta_uv_lambda, ave_removal = ave_removal, $
    diff_ratio = diff_ratio, diff_range = diff_range, $
    png = png, eps = eps, pdf = pdf, data_range = data_range, invert_colorbar = invert_colorbar, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, sim = sim, window_num = window_num, $
    uvf_input = uvf_input, diff_save_path = diff_save_path, plot_path = diff_plot_path, ps_foldernames=ps_foldernames
    
  if n_elements(folder_names) gt 2 then message, 'only 1 or 2 folder_names allowed'
  if n_elements(folder_names) eq 0 then message, 'at least 1 folder name must be specified'
  if n_elements(obs_names_in) gt 2 then message, 'only 1 or 2 obs_names_in allowed'
  if n_elements(spec_window_types) gt 2 then message, 'only 1 or 2 spec_window_types allowed'
  if n_elements(delta_uv_lambda) gt 1 then message, 'only 1 delta_uv_lambda allowed'
  if n_elements(ave_removal) gt 2 then message, 'only 1 or 2 ave_removal values allowed'
  
  spawn, 'hostname', hostname
  if stregex(hostname, 'mit.edu', /boolean) eq 1 then loc_name = 'mit'
  if stregex(hostname, 'enterprise', /boolean) eq 1 then loc_name = 'enterprise'
  case loc_name of
    'mit':  folder_names = mit_folder_locs(folder_names, rts = rts)
    'enterprise': folder_names = enterprise_folder_locs(folder_names, rts = rts)
    else: begin
      folder_test = file_test(folder_names, /directory)
      if min(folder_test) eq 0 then message, 'machine not recognized and folder not found'
    endelse
  endcase
  
  obs_info = ps_filenames(folder_names, obs_names_in, dirty_folder = dirty_folder, exact_obsnames = exact_obsnames, rts = rts, sim = sim, $
    uvf_input = uvf_input, casa = casa, data_subdirs = data_subdirs, $
    ps_foldernames = ps_foldernames, save_paths = save_paths, plot_paths = plot_paths, refresh_info = refresh_info, no_wtvar_rts = no_wtvar_rts)
    
  if n_elements(diff_plot_path) eq 0 then begin
    if n_elements(diff_save_path) gt 0 then diff_plot_path = diff_save_path + path_sep() + 'plots' + path_sep()
  endif
  
  wh_noinfo = where(obs_info.info_files eq '', count_noinfo)
  if count_noinfo gt 0 then message, 'Info files are not all present'
  
  if n_elements(data_range) eq 0 then data_range = [1e-3, 1e1]
  
  ps_ratio_plots, folder_names, obs_info, cube_types, pols, all_pol_diff_ratio = all_pol_diff_ratio, $
    freq_ch_range = freq_ch_range,plot_path = diff_plot_path, plot_filebase = plot_filebase, $
    note = note, spec_window_types = spec_window_types, delta_uv_lambda = delta_uv_lambda, $
    ave_removal = ave_removal, data_range = data_range, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, $
    diff_ratio = diff_ratio, diff_range = diff_range, invert_colorbar = invert_colorbar, $
    plot_wedge_line = plot_wedge_line, quiet = quiet, png = png, eps = eps, pdf = pdf, $
    window_num = window_num
end
