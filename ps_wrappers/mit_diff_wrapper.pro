pro mit_diff_wrapper, folder_names, obs_names_in, cube_types = cube_types, pols = pols, refresh_diff = refresh_diff, $
    spec_window_types = spec_window_types, all_type_pol = all_type_pol, $
    png = png, eps = eps, pdf = pdf, data_range = data_range, data_min_abs = data_min_abs, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, sim = sim, axis_type_1d=axis_type_1d
    
    
  if n_elements(folder_names) gt 2 then message, 'only 1 or 2 folder_names allowed'
  if n_elements(obs_names_in) gt 2 then message, 'only 1 or 2 obs_names_in allowed'
  if n_elements(spec_window_types) gt 2 then message, 'only 1 or 2 spec_window_types allowed'
  
  if n_elements(spec_window_types) eq 2 then $
    if spec_window_types[0] eq spec_window_types[1] then spec_window_types = spec_window_types[0]
    
  if (n_elements(folder_names) eq 2 or n_elements(obs_names_in) eq 2 or n_elements(spec_window_types) eq 2) $
    and n_elements(cube_types) eq 0 and n_elements(pols) eq 0 and n_elements(all_type_pol) eq 0 then all_type_pol = 1
    
  if keyword_set(all_type_pol) and n_elements(folder_names) lt 2 and n_elements(obs_names_in) lt 2 and n_elements(spec_window_types) lt 2 then $
    message, 'all_type_pol can only be set with 2 folder names and/or 2 obs names and/or 2 spec windows'
    
  if not keyword_set(all_type_pol) then begin
    if n_elements(cube_types) eq 0 then if n_elements(folder_names) eq 1 then cube_types = ['dirty', 'res'] else cube_types = 'res'
    if n_elements(pols) eq 0 then pols = 'xx'
  endif else begin
    undefine, cube_types, pols
  endelse
  
  if n_elements(obs_names_in) eq 2 then diff_str = '_' + obs_names_in
  if n_elements(spec_window_types) eq 2 then begin
    type_list = ['Hann', 'Hamming', 'Blackman', 'Nutall', 'Blackman-Nutall', 'Blackman-Harris', 'None']
    sw_tag_list = ['hann', 'ham', 'blm', 'ntl', 'bn', 'bh', '']
    sw_tags = strarr(2)
    for i=0, 1 do begin
      wh_type = where(strlowcase(type_list) eq strlowcase(spec_window_types[i]), count_type)
      if count_type eq 0 then wh_type = where(strlowcase(sw_tag_list) eq strlowcase(spec_window_types[i]), count_type)
      if count_type eq 0 then message, 'Spectral window type not recognized.' else begin
        spec_window_types[i] = type_list[wh_type[0]]
        if spec_window_types[i] eq 'None' then sw_tags[i] = '_nosw' else sw_tags[i] = '_' + sw_tag_list[wh_type[0]]
      endelse
    endfor
    
    if n_elements(diff_str) gt 0 then diff_str = diff_str + sw_tags else diff_str = sw_tags
  endif
  if n_elements(diff_str) eq 0 then diff_str = strarr(2)
  
  if n_elements(folder_names) eq 0 then folder_names = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_apb_pipeline_paper_deep_1/Healpix/'
  
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
  obs_info = ps_filenames(folder_names, obs_names_in, rts = rts, sim = sim, casa = casa, data_subdirs = 'Healpix/', save_paths = save_paths, plot_paths = save_path)
  
  wh_noinfo = where(obs_info.info_files eq '', count_noinfo)
  if count_noinfo gt 0 then message, 'Info files are not all present'
  
  if tag_exist(obs_info, 'diff_note') then obs_info = create_struct(obs_info, 'diff_plot_path', obs_info.diff_save_path)
  
  if n_elements(obs_info.folder_names) eq 2 then begin
    save_path = obs_info.diff_save_path
    note = obs_info.diff_note
    plot_path = obs_info.diff_plot_path
  endif else begin
    save_path = obs_info.save_paths[0]
    note = obs_info.fhd_types[0]
    plot_path = obs_info.plot_paths[0]
  endelse
  
  if n_elements(spec_window_types) eq 2 then note = note + ' ' + spec_window_types[0] + ' minus ' + spec_window_types[1]
  
  if keyword_set(all_type_pol) then begin
  
    if n_elements(folder_names) eq 2 then plot_filebase = obs_info.name_same_parts + '__' + obs_info.name_diff_parts[0] + diff_str[0] + $
      '_minus_' + obs_info.name_diff_parts[1] + diff_str[1] $
    else plot_filebase = obs_info.fhd_types[0] + '__' + diff_str[0] + '_minus' + diff_str[1]
    
  endif else begin
  
    if n_elements(folder_names) eq 1 then begin
      if n_elements(obs_info.obs_names) gt 1 then begin
        plot_filebase = obs_info.fhd_types[0] + '_' + obs_info.obs_names[0] + '_' + cube_types[0] + '_' + pols[0] + diff_str[0] + $
          '_minus_' + obs_info.obs_names[0] + '_' + cube_types[max_type] + '_' + pols[max_pol] + diff_str[1]
      endif else begin
        if obs_info.integrated[0] eq 0 then plot_start = obs_info.fhd_types[0] + '_' + obs_info.obs_names[0] else plot_start = obs_info.fhd_types[0]
        
        plot_filebase = plot_start + '_' + cube_types[0] + '_' + pols[0] + diff_str[0] + $
          '_minus_' + cube_types[max_type] + '_' + pols[max_pol] + diff_str[1]
      endelse
    endif else begin
      plot_filebase = obs_info.name_same_parts + '__' + strjoin([obs_info.name_diff_parts[0], cube_types[0], pols[0]], '_') + diff_str[0]  + $
        '_minus_' + strjoin([obs_info.name_diff_parts[1], cube_types[n_elements(cube_types)-1], pols[n_elements(pols)-1]], '_') + diff_str[1]
    endelse
  endelse
  
  if not file_test(plot_path, /directory) then file_mkdir, plot_path
  if not file_test(save_path, /directory) then file_mkdir, save_path
  
  ps_difference_plots, obs_info.info_files, cube_types, pols, spec_window_types = spec_window_types, all_type_pol = all_type_pol, refresh_diff = refresh_diff, $
    plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    note = note, kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, axis_type_1d=axis_type_1d, $
    data_range = data_range, data_min_abs = data_min_abs, $
    quiet = quiet, png = png, eps = eps, pdf = pdf
    
end
