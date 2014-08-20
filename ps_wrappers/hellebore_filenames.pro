function hellebore_filenames, folder_names, obs_names_in, rts = rts, sim = sim, casa = casa

  if keyword_set(rts) then begin
    std_savepath = base_path('data') + 'rts_data/'
    std_plotpath = base_path('plots') + 'power_spectrum/rts_data/'
  endif else if keyword_set(casa) then begin
    std_savepath = base_path('data') + 'mit_data/'
    std_plotpath = base_path('plots') + 'power_spectrum/mit_data/'
  endif else if keyword_set(sim) then begin
    std_savepath = base_path('data') + 'fhd_sim_data/'
    std_plotpath = base_path('plots') + 'power_spectrum/fhd_sim/'
  endif else begin
    std_savepath = base_path('data') + 'fhd_ps_data/'
    std_plotpath = base_path('plots') + 'power_spectrum/fhd_data/'
  endelse
  
  if keyword_set(rts) then begin
    plot_paths = strarr(n_elements(folder_names))
    
    for i=0, n_elements(folder_names)-1 do begin
      ;; check for folder existence, otherwise look for common folder names to figure out full path. If none found, try base_path('data') + 'rts_data/'
      folder_test = file_test(folder_names[i], /directory)
      if folder_test eq 0 then begin
        pos_rts_data = strpos(folder_names[i], 'rts_data')
        if pos_rts_data gt -1 then begin
          test_name = base_path('data') + strmid(folder_names[i], pos_rts_data)
          folder_test = file_test(test_name, /directory)
          if folder_test eq 1 then folder_names[i] = test_name
        endif
      endif
      if folder_test eq 0 then begin
        test_name = base_path('data') + 'rts_data/' + folder_names[i]
        folder_test = file_test(test_name, /directory)
        if folder_test eq 1 then folder_names[i] = test_name
      endif
      
      if folder_test eq 0 then message, 'folder not found'
      
      pos = strpos(folder_names[i], std_savepath)
      if pos ne -1 then save_path_ext = strmid(folder_names[i], pos + strlen(std_savepath)) $
      else save_path_ext = ''
      plot_paths[i] = std_plotpath + save_path_ext + path_sep()
    endfor
    
  endif else if keyword_set(casa) then begin
    plot_paths = strarr(n_elements(folder_names))
    
    for i=0, n_elements(folder_names)-1 do begin
      ;; check for folder existence, otherwise look for common folder names to figure out full path. If none found, try base_path('data') + 'mit_data/'
      folder_test = file_test(folder_names[i], /directory)
      if folder_test eq 0 then begin
        pos_rts_data = strpos(folder_names[i], 'mit_data')
        if pos_rts_data gt -1 then begin
          test_name = base_path('data') + strmid(folder_names[i], pos_rts_data)
          folder_test = file_test(test_name, /directory)
          if folder_test eq 1 then folder_names[i] = test_name
        endif
      endif
      if folder_test eq 0 then begin
        test_name = base_path('data') + 'mit_data/' + folder_names[i]
        folder_test = file_test(test_name, /directory)
        if folder_test eq 1 then folder_names[i] = test_name
      endif
      
      if folder_test eq 0 then message, 'folder not found'
      
      pos = strpos(folder_names[i], std_savepath)
      if pos ne -1 then save_path_ext = strmid(folder_names[i], pos + strlen(std_savepath)) $
      else save_path_ext = ''
      plot_paths[i] = std_plotpath + save_path_ext + path_sep()
    endfor
    
  endif else begin
    plot_paths = strarr(n_elements(folder_names))
    
    for i=0, n_elements(folder_names)-1 do begin
      ;; check for folder existence, otherwise look for common folder names to figure out full path. If none found, try base_path('data') + 'fhd_ps_data/128T_cubes/'
      if keyword_set(sim) then begin
        folder_test = file_test(folder_names[i], /directory)
        if folder_test eq 0 then begin
          pos_sim_data = strpos(folder_names[i], 'fhd_sim_data')
          if pos_sim_data gt -1 then begin
            test_name = base_path('data') + strmid(folder_names[i], pos_sim_data)
            folder_test = file_test(test_name, /directory)
            if folder_test eq 1 then folder_names[i] = test_name
          endif
        endif
        if folder_test eq 0 then begin
          test_name = base_path('data') + 'fhd_sim_data/' + folder_names[i]
          folder_test = file_test(test_name, /directory)
          if folder_test eq 1 then folder_names[i] = test_name
        endif
      endif else begin
        folder_test = file_test(folder_names[i], /directory)
        if folder_test eq 0 then begin
          pos_fhd_data = strpos(folder_names[i], 'fhd_ps_data')
          if pos_fhd_data gt -1 then begin
            test_name = base_path('data') + strmid(folder_names[i], pos_fhd_data)
            folder_test = file_test(test_name, /directory)
            if folder_test eq 1 then folder_names[i] = test_name
          endif
        endif
        if folder_test eq 0 then begin
          pos_fhd_128 = strpos(folder_names[i], '128T_cubes')
          if pos_fhd_128 gt -1 then begin
            test_name = base_path('data') + 'fhd_ps_data/' + strmid(folder_names[i], pos_fhd_128)
            folder_test = file_test(test_name, /directory)
            if folder_test eq 1 then folder_names[i] = test_name
          endif
        endif
        if folder_test eq 0 then begin
          test_name = base_path('data') + 'fhd_ps_data/128T_cubes/' + folder_names[i]
          folder_test = file_test(test_name, /directory)
          if folder_test eq 1 then folder_names[i] = test_name
        endif
      endelse
      
      if folder_test eq 0 then message, 'folder not found'
      
      pos = strpos(folder_names[i], std_savepath)
      if pos ne -1 then save_path_ext = strmid(folder_names[i], pos + strlen(std_savepath)) $
      else save_path_ext = ''
      plot_paths[i] = std_plotpath + save_path_ext + path_sep()
    endfor
    
  endelse
  
  obs_info = ps_filenames(folder_names, obs_names_in, rts = rts, sim = sim, casa = casa, plot_paths = plot_paths)
  
  if tag_exist(obs_info, 'diff_note') then begin
    pos = strpos(obs_info.diff_save_path, std_savepath)
    if pos ne -1 then diff_save_path_ext = strmid(obs_info.diff_save_path, pos + strlen(std_savepath)) $
    else diff_save_path_ext = ''
    diff_plot_path = base_path('plots') + 'power_spectrum/fhd_data/' + diff_save_path_ext
    
    obs_info = create_struct(obs_info, 'diff_plot_path', diff_plot_path)
  endif
  
  return, obs_info
  
end
