pro mit_cube_images, folder_names, obs_names_in, exact_obsnames = exact_obsnames, data_subdirs=data_subdirs, cube_types = cube_types, $
    pols = pols, evenodd = evenodd, $
    rts = rts, sim = sim, casa = casa, png = png, eps = eps, pdf = pdf, slice_range = slice_range, sr2 = sr2, $
    nvis_norm = nvis_norm, ratio = ratio, diff_ratio = diff_ratio, diff_frac = diff_frac, $
    log = log, data_range = data_range, color_profile = color_profile, sym_color = sym_color, window_num = window_num, plot_as_map = plot_as_map
    
  message, 'Error: This wrapper is depricated and will not continue to be supported. Please call ps_cube_images instead. ' + $
    'This line can be commented out to allow the deprecacated code to run.'
    
  if n_elements(folder_names) gt 2 then message, 'No more than 2 folder_names can be supplied'
  if n_elements(evenodd) eq 0 then evenodd = 'even'
  if n_elements(evenodd) gt 2 then message, 'No more than 2 evenodd values can be supplied'
  if n_elements(obs_names_in) gt 2 then message, 'No more than 2 obs_names can be supplied'
  
  folder_names = mit_folder_locs(folder_names, rts = rts)
  
  obs_info = ps_filenames(folder_names, obs_names_in, exact_obsnames = exact_obsnames, rts = rts, sim = sim, casa = casa, $
    data_subdirs = data_subdirs, save_paths = save_paths, plot_paths = plot_paths)
    
  cube_images, folder_names, obs_info, nvis_norm = nvis_norm, pols = pols, cube_types = cube_types, evenodd = evenodd, $
    png = png, eps = eps, pdf = pdf, slice_range = slice_range, sr2 = sr2, $
    ratio = ratio, diff_ratio = diff_ratio, diff_frac = diff_frac, $
    log = log, data_range = data_range, color_profile = color_profile, sym_color = sym_color, $
    window_num = window_num, plot_as_map = plot_as_map
    
end