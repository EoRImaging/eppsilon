pro hellebore_ratio_wrapper, folder_names, obs_names_in, cube_types = cube_types, pols = pols, $
    spec_window_types = spec_window_types, diff_ratio = diff_ratio, diff_range = diff_range, $
    png = png, eps = eps, pdf = pdf, data_range = data_range, data_min_abs = data_min_abs, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, sim = sim, window_num = window_num
    
  if n_elements(folder_names) gt 2 then message, 'only 1 or 2 folder_names allowed'
  if n_elements(folder_names) eq 0 then message, 'at least 1 folder name must be specified'
  if n_elements(obs_names_in) gt 2 then message, 'only 1 or 2 obs_names_in allowed'
  if n_elements(spec_window_types) gt 2 then message, 'only 1 or 2 spec_window_types allowed'
  
  
  obs_info = hellebore_filenames(folder_names, obs_names_in, sim = sim)
  wh_noinfo = where(obs_info.info_files eq '', count_noinfo)
  if count_noinfo gt 0 then message, 'Info files are not all present'
  
  if n_elements(data_range) eq 0 then data_range = [1e-3, 1e1]
  
  ps_ratio_plots, folder_names, obs_info, cube_types, pols, $
    plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    note = note, spec_window_types = spec_window_types, data_range = data_range, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, diff_ratio = diff_ratio, diff_range = diff_range, $
    plot_wedge_line = plot_wedge_line, quiet = quiet, png = png, eps = eps, pdf = pdf, window_num = window_num
end
