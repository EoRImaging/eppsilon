pro mit_cube_images, folder_names, obs_names_in, exact_obsnames = exact_obsnames, data_subdirs=data_subdirs, cube_types = cube_types, $
    pols = pols, evenodd = evenodd, $
    rts = rts, sim = sim, casa = casa, png = png, eps = eps, pdf = pdf, slice_range = slice_range, sr2 = sr2, $
    nvis_norm = nvis_norm, ratio = ratio, diff_ratio = diff_ratio, diff_frac = diff_frac, $
    log = log, data_range = data_range, color_profile = color_profile, sym_color = sym_color, window_num = window_num, plot_as_map = plot_as_map
    
  if n_elements(folder_names) gt 2 then message, 'No more than 2 folder_names can be supplied'
  if n_elements(evenodd) eq 0 then evenodd = 'even'
  if n_elements(evenodd) gt 2 then message, 'No more than 2 evenodd values can be supplied'
  if n_elements(obs_names_in) gt 2 then message, 'No more than 2 obs_names can be supplied'
  
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
    
    start_path2 = '/nfs/mwa-03/r1/'
    if folder_test eq 0 then begin
      pos_eor2013 = strpos(folder_names[i], 'EoR2013')
      if pos_eor2013 gt -1 then begin
        test_name = start_path2 + strmid(folder_names[i], pos_eor2013)
        folder_test = file_test(test_name, /directory)
        if folder_test eq 1 then folder_names[i] = test_name
      endif
    endif
    if folder_test eq 0 then begin
      test_name = start_path2 + 'EoR2013/' + folder_names[i]
      folder_test = file_test(test_name, /directory)
      if folder_test eq 1 then folder_names[i] = test_name
    endif
    
    if folder_test eq 0 then message, 'folder not found'
  endfor
  
  save_paths = folder_names + '/ps/'
  plot_paths = save_paths + 'plots/'
  if n_elements(data_subdirs) eq 0 then data_subdirs = 'Healpix/' else if n_elements(data_subdirs) gt 2 then message, 'No more than 2 data_subdirs can be supplied.'
  obs_info = ps_filenames(folder_names, obs_names_in, exact_obsnames = exact_obsnames, rts = rts, sim = sim, casa = casa, $
    data_subdirs = data_subdirs, save_paths = save_paths, plot_paths = plot_paths)
    
  cube_images, folder_names, obs_info, nvis_norm = nvis_norm, pols = pols, cube_types = cube_types, evenodd = evenodd, $
    png = png, eps = eps, pdf = pdf, slice_range = slice_range, sr2 = sr2, $
    ratio = ratio, diff_ratio = diff_ratio, diff_frac = diff_frac, $
    log = log, data_range = data_range, color_profile = color_profile, sym_color = sym_color, $
    window_num = window_num, plot_as_map = plot_as_map
    
end