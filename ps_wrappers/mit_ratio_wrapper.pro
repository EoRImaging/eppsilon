pro mit_ratio_wrapper, folder_names, obs_names_in, exact_obsnames = exact_obsnames, cube_types = cube_types, pols = pols, $
    spec_window_types = spec_window_types, diff_ratio = diff_ratio, diff_range = diff_range, $
    png = png, eps = eps, pdf = pdf, data_range = data_range, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, sim = sim, window_num = window_num
    
    
  if n_elements(folder_names) gt 2 then message, 'only 1 or 2 folder_names allowed'
  if n_elements(folder_names) eq 0 then message, 'at least 1 folder name must be specified'
  if n_elements(obs_names_in) gt 2 then message, 'only 1 or 2 obs_names_in allowed'
  if n_elements(spec_window_types) gt 2 then message, 'only 1 or 2 spec_window_types allowed'
  
  
  for i=0, n_elements(folder_names)-1 do begin
    ;; check for folder existence, otherwise look for common folder names to figure out full path. If none found, try '/nfs/mwa-09/r1/djc/EoR2013/Aug23/'
    start_path = '/nfs/mwa-09/r1/djc/'
    folder_test = file_test(folder_names[i], /directory)
    if folder_test eq 0 then begin
      pos_eor2013 = strpos(folder_names[i], 'EoR2013')
      if pos_eor2013 gt -1 then begin
        test_name = start_path + strmid(folder_names[i], pos_eor2013)
        folder_test = file_test(test_name, /directory)
        if folder_test eq 1 then folder_names[i] = test_name
      endif
    endif
    if folder_test eq 0 then begin
      pos_aug23 = strpos(folder_names[i], 'Aug23')
      if pos_aug23 gt -1 then begin
        test_name = start_path + 'EoR2013/' + strmid(folder_names[i], pos_aug23)
        folder_test = file_test(test_name, /directory)
        if folder_test eq 1 then folder_names[i] = test_name
      endif
    endif
    if folder_test eq 0 then begin
      pos_aug26 = strpos(folder_names[i], 'Aug26')
      if pos_aug26 gt -1 then begin
        test_name = start_path + 'EoR2013/' + strmid(folder_names[i], pos_aug26)
        folder_test = file_test(test_name, /directory)
        if folder_test eq 1 then folder_names[i] = test_name
      endif
    endif
    if folder_test eq 0 then begin
      pos_week1 = strpos(folder_names[i], 'week1')
      if pos_week1 gt -1 then begin
        test_name = start_path + 'EoR2013/' + strmid(folder_names[i], pos_week1)
        folder_test = file_test(test_name, /directory)
        if folder_test eq 1 then folder_names[i] = test_name
      endif
    endif
    if folder_test eq 0 then begin
      test_name = start_path + 'EoR2013/Aug23/' + folder_names[i]
      folder_test = file_test(test_name, /directory)
      if folder_test eq 1 then folder_names[i] = test_name
    endif
    if folder_test eq 0 then begin
      test_name = start_path + 'EoR2013/Aug26/' + folder_names[i]
      folder_test = file_test(test_name, /directory)
      if folder_test eq 1 then folder_names[i] = test_name
    endif
    if folder_test eq 0 then begin
      test_name = start_path + 'EoR2013/week1/' + folder_names[i]
      folder_test = file_test(test_name, /directory)
      if folder_test eq 1 then folder_names[i] = test_name
    endif
    
    if folder_test eq 0 then message, 'folder not found'
  endfor
  
  save_paths = folder_names + '/ps/'
  obs_info = ps_filenames(folder_names, obs_names_in, exact_obsnames = exact_obsnames, rts = rts, sim = sim, casa = casa, $
    data_subdirs = 'Healpix/', save_paths = save_paths, plot_paths = save_path)
    
  wh_noinfo = where(obs_info.info_files eq '', count_noinfo)
  if count_noinfo gt 0 then message, 'Info files are not all present'
  
  if n_elements(data_range) eq 0 then data_range = [1e-3, 1e1]
  
  ps_ratio_plots, folder_names, obs_info, cube_types, pols, $
    plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    note = note, spec_window_types = spec_window_types, data_range = data_range, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, diff_ratio = diff_ratio, diff_range = diff_range, $
    plot_wedge_line = plot_wedge_line, quiet = quiet, png = png, eps = eps, pdf = pdf, window_num = window_num
end
