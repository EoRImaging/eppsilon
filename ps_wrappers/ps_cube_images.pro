pro ps_cube_images, folder_names, obs_names_in, exact_obsnames = exact_obsnames, data_subdirs=data_subdirs, cube_types = cube_types, $
    pols = pols, evenodd = evenodd, $
    rts = rts, sim = sim, casa = casa, png = png, eps = eps, pdf = pdf, slice_range = slice_range, sr2 = sr2, $
    nvis_norm = nvis_norm, ratio = ratio, diff_ratio = diff_ratio, diff_frac = diff_frac, $
    log = log, data_range = data_range, data_min_abs = data_min_abs, color_profile = color_profile, sym_color = sym_color, $
    window_num = window_num, plot_as_map = plot_as_map, plot_path = plot_path
    
  ;; default to normalizing by nvis
  if n_elements(nvis_norm) eq 0 then nvis_norm = 1
    
  if n_elements(folder_names) gt 2 then message, 'No more than 2 folder_names can be supplied'
  if n_elements(evenodd) eq 0 then evenodd = 'even'
  if n_elements(evenodd) gt 2 then message, 'No more than 2 evenodd values can be supplied'
  if n_elements(obs_names_in) gt 2 then message, 'No more than 2 obs_names can be supplied'
  
  spawn, 'hostname', hostname
  if stregex(hostname, 'mit.edu', /boolean) eq 1 then loc_name = 'mit'
  if stregex(hostname, 'enterprise', /boolean) eq 1 then loc_name = 'enterprise'
  if stregex(hostname, 'constellation', /boolean) eq 1 then loc_name = 'enterprise'
  if stregex(hostname, 'defiant', /boolean) eq 1 then loc_name = 'enterprise'
  case loc_name of
    'mit':  folder_names = mit_folder_locs(folder_names, rts = rts)
    'enterprise': folder_names = enterprise_folder_locs(folder_names, rts = rts)
    else: begin
      folder_test = file_test(folder_names, /directory)
      if min(folder_test) eq 0 then message, 'machine not recognized and folder not found'
    endelse
  endcase
  
  obs_info = ps_filenames(folder_names, obs_names_in, dirty_folder = dirty_folder, exact_obsnames = exact_obsnames, rts = rts, $
    uvf_input = uvf_input, casa = casa, data_subdirs = data_subdirs, $
    save_paths = save_paths, plot_paths = plot_paths, refresh_info = refresh_info, no_wtvar_rts = no_wtvar_rts)
    
  if n_elements(plot_path) eq 0 then if tag_exist(obs_info, 'diff_plot_path') then plot_path = obs_info.diff_plot_path $
  else plot_path = obs_info.plot_paths[0]
  
  if keyword_set(rts) then begin
  
    if obs_info.info_files[0] ne '' then datafile = obs_info.info_files[0] else $
      if obs_info.cube_files.(0)[0] ne '' then datafile = obs_info.cube_files.(0) else $
      datafile = rts_fits2idlcube(obs_info.datafiles.(0), obs_info.weightfiles.(0), obs_info.variancefiles.(0), $
      pol_inc, save_path = obs_info.folder_names[0]+path_sep(), refresh = refresh_dft, no_wtvar = no_wtvar_rts)
      
    if keyword_set(refresh_rtscube) then datafile = rts_fits2idlcube(obs_info.datafiles.(0), obs_info.weightfiles.(0), obs_info.variancefiles.(0), $
      pol_inc, save_path = obs_info.folder_names[0]+path_sep(), /refresh, no_wtvar = no_wtvar_rts)
      
    if keyword_set(no_wtvar_rts) then stop
    
  endif
  
  cube_images, folder_names, obs_info, nvis_norm = nvis_norm, pols = pols, cube_types = cube_types, evenodd = evenodd, rts = rts, $
    png = png, eps = eps, pdf = pdf, slice_range = slice_range, sr2 = sr2, $
    ratio = ratio, diff_ratio = diff_ratio, diff_frac = diff_frac, $
    log = log, data_range = data_range, data_min_abs = data_min_abs, color_profile = color_profile, sym_color = sym_color, $
    window_num = window_num, plot_as_map = plot_as_map, plot_path = plot_pat
    
end