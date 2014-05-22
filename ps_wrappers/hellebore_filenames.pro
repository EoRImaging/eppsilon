function hellebore_filenames, folder_names, obs_names_in, rts = rts, sim = sim, casa = casa

  n_filesets = max([n_elements(folder_names), n_elements(obs_names_in)])
  
  if n_filesets gt 1 then begin
    if n_elements(folder_names) eq 1 then folder_names = replicate(folder_names, n_filesets)
    if n_elements(folder_names) ne n_filesets then message, 'If both folder_names and obs_names_in are arrays, the number of elements must match'
    
    if n_elements(obs_names_in) eq 1 then obs_names_in = replicate(obs_names_in, n_filesets)
    if n_elements(obs_names_in) gt 0 and n_elements(obs_names_in) ne n_filesets then message, 'If both folder_names and obs_names_in are arrays, the number of elements must match'
  endif
  
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
    obs_names = strarr(n_filesets)
    rts_types = strarr(n_filesets)
    info_files = strarr(n_filesets)
    ;integrated = intarr(n_filesets)
    plot_paths = strarr(n_filesets)
    
    for i=0, n_filesets-1 do begin
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
      
      rts_types[i] = file_basename(folder_names[i])
      
      pos = strpos(folder_names[i], std_savepath)
      if pos ne -1 then save_path_ext = strmid(folder_names[i], pos + strlen(std_savepath)) $
      else save_path_ext = ''
      plot_paths[i] = std_plotpath + save_path_ext + path_sep()
      
      if n_elements(obs_names_in) gt 0 then begin
        if size(obs_names_in,/type) eq 7 then begin
          obs_names[i] = obs_names_in[i]
          obs_name_single = obs_names[i]
        endif else begin
          obs_names[i] = number_formatter(obs_names_in[i])
          obs_name_single = obs_names[i]
        endelse
      endif else begin
        obs_names[i] = ''
        obs_name_single = ''
      endelse
      
      ;; first look for info files
      info_file = file_search(folder_names[i] + '/' + obs_names[i] + '*info*', count = n_infofile)
      if n_infofile gt 0 then begin
        if obs_names[i] eq '' then begin
          if n_infofile gt 1 then print, 'More than 1 info files found, using first one'
          info_files[i] = info_file[0]
          obs_names[i] = stregex(info_files[i], '[0-9]+.[0-9]+_', /extract)
        endif else begin
          if n_infofile gt 1 then message, 'More than one info file found with given obs_range'
          info_files[i] = info_file[0]
        endelse
        
      endif
      
      ;; then look for cube files
      cube_file_list = file_search(folder_names[i] + '/' + obs_names[i] + '*_image*.fits', count = n_cubefiles)
      if n_cubefiles gt 0 then begin
        if obs_names[i] eq '' then begin
          obs_name_arr = stregex(cube_file_list, '[0-9]+.[0-9]+_', /extract)
          wh_first = where(obs_name_arr eq obs_name_arr[0], count_first)
          if count_first lt n_cubefiles then $
            print, 'More than one obs_range found, using first range (' + obs_name_arr[0] + ', ' + number_formatter(count_first) + ' files)'
          datafiles = cube_file_list[wh_first]
          obs_names[i] = obs_name_arr[0]
        endif else begin
          datafiles = cube_file_list
        endelse
        
      endif
      
      ;; now get weights & variance files
      weightfile_list = file_search(folder_names[i] + '/' + obs_names[i] + '*_weights*.fits', count = n_wtfiles)
      if n_wtfiles ne n_elements(datafiles) and info_files[i] eq '' then message, 'number of weight files does not match number of datafiles'
      
      variancefile_list = file_search(folder_names[i] + '/' + obs_names[i] + '*_weights*.fits', count = n_varfiles)
      if n_varfiles ne n_elements(datafiles) and info_files[i] eq '' then message, 'number of variance files does not match number of datafiles'
      
      if n_elements(datafiles) eq 0 and info_files[i] eq '' then message, 'No cube or info files found in folder ' + folder_name
      
      if n_elements(datafiles) eq 0 then begin
        datafiles = ''
        weightfile_list = ''
        variancefile_list = ''
      endif
      
      if i eq 0 then begin
        tag = 'fs' + number_formatter(i)
        cube_files = create_struct(tag, datafiles)
        weightfiles = create_struct(tag, weightfile_list)
        variancefiles = create_struct(tag, variancefile_list)
      endif else begin
        tag = 'fs' + number_formatter(i)
        cube_files = create_struct(cube_files, tag, datafile)
        weightfiles = create_struct(weightfiles, tag, weightfile_list)
        variancefiles = create_struct(variancefiles, tag, variancefile_list)
      endelse
      undefine, datafiles, weightfile_list, variancefile_list
      
    endfor
    
    obs_info = {folder_names:folder_names, obs_names:obs_names, info_files:info_files, cube_files:cube_files, $
      weightfiles:weightfiles, variancefiles:variancefiles, rts_types:rts_types, plot_paths:plot_paths}
      
  endif else if keyword_set(casa) then begin
    obs_names = strarr(n_filesets)
    casa_types = strarr(n_filesets)
    info_files = strarr(n_filesets)
    ;integrated = intarr(n_filesets)
    plot_paths = strarr(n_filesets)
    
    for i=0, n_filesets-1 do begin
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
      
      casa_types[i] = file_basename(folder_names[i])
      
      pos = strpos(folder_names[i], std_savepath)
      if pos ne -1 then save_path_ext = strmid(folder_names[i], pos + strlen(std_savepath)) $
      else save_path_ext = ''
      plot_paths[i] = std_plotpath + save_path_ext + path_sep()
      
      if n_elements(obs_names_in) gt 0 then begin
        if size(obs_names_in,/type) eq 7 then begin
          obs_names[i] = obs_names_in[i]
          obs_name_single = obs_names[i]
        endif else begin
          obs_names[i] = number_formatter(obs_names_in[i])
          obs_name_single = obs_names[i]
        endelse
      endif else begin
        obs_names[i] = ''
        obs_name_single = ''
      endelse
      
      ;; first look for info files
      info_file = file_search(folder_names[i] + '/' + obs_names[i] + '*info*', count = n_infofile)
      if n_infofile gt 0 then begin
        if obs_names[i] eq '' then begin
          if n_infofile gt 1 then print, 'More than 1 info files found, using first one'
          info_files[i] = info_file[0]
        ;obs_names[i] = stregex(info_files[i], '[0-9]+.[0-9]+_', /extract)
        endif else begin
          if n_infofile gt 1 then message, 'More than one info file found with given obs_range'
          info_files[i] = info_file[0]
        endelse
        
      endif
      
      ;; then look for cube files
      cube_file_list = file_search(folder_names[i] + '/' + obs_names[i] + '*_holo_[xy]*.fits', count = n_cubefiles)
      if n_cubefiles gt 0 then begin
        ;        if obs_names[i] eq '' then begin
        ;          obs_name_arr = stregex(cube_file_list, '[0-9]+.[0-9]+_', /extract)
        ;          wh_first = where(obs_name_arr eq obs_name_arr[0], count_first)
        ;          if count_first lt n_cubefiles then $
        ;            print, 'More than one obs_range found, using first range (' + obs_name_arr[0] + ', ' + number_formatter(count_first) + ' files)'
        ;          datafiles = cube_file_list[wh_first]
        ;          obs_names[i] = obs_name_arr[0]
        ;        endif else begin
        datafiles = cube_file_list
      ;        endelse
        
      endif
      
      ;; now get weights & variance files
      weightfile_list = file_search(folder_names[i] + '/' + obs_names[i] + '*_psf_*.fits', count = n_wtfiles)
      if n_wtfiles ne n_elements(datafiles) and info_files[i] eq '' then message, 'number of weight files does not match number of datafiles'
      
      variancefile_list = file_search(folder_names[i] + '/' + obs_names[i] + '*psfbeamsquare*.fits', count = n_varfiles)
      if n_varfiles ne n_elements(datafiles) and info_files[i] eq '' then message, 'number of variance files does not match number of datafiles'
      
      if n_elements(datafiles) eq 0 and info_files[i] eq '' then message, 'No cube or info files found in folder ' + folder_names[i]
      
      if n_elements(datafiles) eq 0 then begin
        datafiles = ''
        weightfile_list = ''
        variancefile_list = ''
      endif
      
      if i eq 0 then begin
        tag = 'fs' + number_formatter(i)
        cube_files = create_struct(tag, datafiles)
        weightfiles = create_struct(tag, weightfile_list)
        variancefiles = create_struct(tag, variancefile_list)
      endif else begin
        tag = 'fs' + number_formatter(i)
        cube_files = create_struct(cube_files, tag, datafiles)
        weightfiles = create_struct(weightfiles, tag, weightfile_list)
        variancefiles = create_struct(variancefiles, tag, variancefile_list)
      endelse
      undefine, datafiles, weightfile_list, variancefile_list
      
    endfor
    
    obs_info = {folder_names:folder_names, obs_names:obs_names, info_files:info_files, cube_files:cube_files, $
      weightfiles:weightfiles, variancefiles:variancefiles, casa_types:casa_types, plot_paths:plot_paths}
      
  endif else begin
    obs_names = strarr(n_filesets)
    fhd_types = strarr(n_filesets)
    info_files = strarr(n_filesets)
    integrated = intarr(n_filesets)
    plot_paths = strarr(n_filesets)
    
    for i=0, n_filesets-1 do begin
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
      
      fhd_types[i] = file_basename(folder_names[i])
      
      pos = strpos(folder_names[i], std_savepath)
      if pos ne -1 then save_path_ext = strmid(folder_names[i], pos + strlen(std_savepath)) $
      else save_path_ext = ''
      plot_paths[i] = std_plotpath + save_path_ext + path_sep()
      
      if n_elements(obs_names_in) gt 0 then begin
        if size(obs_names_in,/type) eq 7 then begin
          obs_names[i] = obs_names_in[i]
          obs_name_single = obs_names[i]
        endif else begin
          obs_names[i] = number_formatter(obs_names_in[i])
          obs_name_single = obs_names[i]
        endelse
      endif else begin
        obs_names[i] = ''
        obs_name_single = ''
      endelse
      
      ;; first look for integrated info files with names like Combined_obs_...
      info_file = file_search(folder_names[i] + '/Combined_obs_' + obs_names[i] + '*info*', count = n_infofile)
      if n_infofile gt 0 then begin
        if obs_names[i] eq '' then begin
          if n_infofile gt 1 then print, 'More than 1 info files found, using first one'
          info_files[i] = info_file[0]
          if stregex(info_files[i], '[0-9]+-[0-9]+', /boolean) then obs_names[i] = stregex(info_files[i], '[0-9]+-[0-9]+', /extract) else begin
            start_pos = strpos(info_file, 'Combined_obs_') + strlen('Combined_obs_')
            end_pos = strpos(strmid(info_file, start_pos), '_cube')
            obs_names[i] = strmid(info_file, start_pos, end_pos)
          endelse
        endif else begin
          if n_infofile gt 1 then message, 'More than one info file found with given obs_range'
          info_files[i] = info_file
          if stregex(info_files[i], '[0-9]+-[0-9]+', /boolean) then obs_names[i] = stregex(info_files[i], '[0-9]+-[0-9]+', /extract) else begin
            start_pos = strpos(info_file, 'Combined_obs_') + strlen('Combined_obs_')
            end_pos = strpos(strmid(info_file, start_pos), '_cube')
            obs_names[i] = strmid(info_file, start_pos, end_pos)
          endelse
        endelse
        integrated[i]=1
        
      endif else if n_elements(obs_range) lt 2 then begin
        ;; then look for single obs info files
        info_file = file_search(folder_names[i] + '/' + obs_name_single + '*info*', count = n_infofile)
        if n_infofile gt 0 then begin
          info_basename = file_basename(info_file)
          if obs_names[i] eq '' then begin
            if n_infofile gt 1 then print, 'More than 1 info files found, using first one'
            info_files[i] = info_file[0]
            if stregex(info_basename[0], '[0-9]+', /boolean) then obs_names[i] = stregex(info_basename[0], '[0-9]+', /extract) else begin
              end_pos = strpos(info_basename[0], '_cube')
              obs_names[i] = strmid(info_basename[0], 0, end_pos)
            endelse
          endif else begin
            if n_infofile gt 1 then message, 'More than one info file found with given obs_range'
            info_files[i] = info_file
            if stregex(info_basename[0], '[0-9]+', /boolean) then obs_names[i] = stregex(info_basename[0], '[0-9]+', /extract) else begin
              end_pos = strpos(info_basename[0], '_cube')
              obs_names[i] = strmid(info_basename[0], 0, end_pos)
            endelse
          endelse
          integrated[i]=0
          
        endif
      endif
      
      ;; first look for integrated cube files with names like Combined_obs_...
      cube_file_list = file_search(folder_names[i] + '/Combined_obs_' + obs_names[i] + '*_cube.sav', count = n_cubefiles)
      if n_cubefiles gt 0 then begin
        if obs_names[i] eq '' then begin
          obs_name_arr = strarr(n_cubefiles)
          for j=0, n_cubefiles-1 do begin
            start_pos = strpos(cube_file_list[j], 'Combined_obs_') + strlen('Combined_obs_')
            end_pos_even = strpos(strmid(cube_file_list[j], start_pos), '_even')
            end_pos_odd = strpos(strmid(cube_file_list[j], start_pos), '_odd')
            end_pos_cube = strpos(strmid(cube_file_list[j], start_pos), '_cube') ;; always > -1
            end_pos = end_pos_even > end_pos_odd
            wh_noend = where(end_pos eq -1, count_noend)
            if count_noend gt 0 then end_pos[wh_noend] = end_pos_cube[wh_noend]
            
            ;obs_name_arr = stregex(cube_file_list, '[0-9]+-[0-9]+', /extract)
            obs_name_arr[j] = strmid(cube_file_list[j], start_pos, end_pos)
          endfor
          
          wh_first = where(obs_name_arr eq obs_name_arr[0], count_first)
          if count_first lt n_elements(cube_file_list) then $
            print, 'More than one obs_range found, using first range (' + obs_name_arr[0] + ', ' + number_formatter(count_first) + ' files)'
          if count_first gt 2 then message, 'More than two cubes found with first obs_range'
          datafile = cube_file_list[wh_first]
          obs_names[i] = obs_name_arr[0]
          integrated[i]=1
        endif else begin
          if n_elements(cube_file_list) gt 2 then message, 'More than two cubes found with given obs_name'
          datafile = cube_file_list
        endelse
        
      endif else begin
        ;; then look for single obs cube files
        cube_file_list = file_search(folder_names[i] + '/' + obs_names[i] + '*_cube.sav', count = n_cubefiles)
        if n_cubefiles gt 0 then begin
          cube_basename = file_basename(cube_file_list)
          if obs_names[i] eq '' then begin
            obs_name_arr = strarr(n_cubefiles)
            for j=0, n_cubefiles-1 do begin
              end_pos_even = strpos(strmid(cube_basename[j], 0), '_even')
              end_pos_odd = strpos(strmid(cube_basename[j], 0), '_odd')
              end_pos_cube = strpos(strmid(cube_basename[j], 0), '_cube') ;; always > -1
              end_pos = end_pos_even > end_pos_odd
              wh_noend = where(end_pos eq -1, count_noend)
              if count_noend gt 0 then end_pos[wh_noend] = end_pos_cube[wh_noend]
              
              obs_name_arr[j] = strmid(cube_basename[j], 0, end_pos)
            endfor
            ;obs_name_arr = stregex(cube_basename, '[0-9]+', /extract)
            
            wh_first = where(obs_name_arr eq obs_name_arr[0], count_first)
            if count_first lt n_elements(cube_file_list) then $
              print, 'More than one obs_range found, using first range (' + obs_name_arr[0] + ', ' + number_formatter(count_first) + ' files)'
            if count_first gt 2 then message, 'More than two cubes found with first obs_range'
            datafile = cube_file_list[wh_first]
            obs_names[i] = obs_name_arr[0]
            integrated[i]=0
          endif else begin
            if n_elements(cube_file_list) gt 2 then message, 'More than two cubes found with given obs_range'
            datafile = cube_file_list
          endelse
        endif
        
      endelse
      
      if n_elements(datafile) eq 0 and info_files[i] eq '' then message, 'No cube or info files found in folder ' + folder_name
      
      if n_elements(datafile) eq 0 then datafile = ''
      if i eq 0 then cube_files = create_struct('fs0', datafile) else begin
        tag = 'fs' + number_formatter(i)
        cube_files = create_struct(cube_files, tag, datafile)
      endelse
      undefine, datafile
      
    endfor
    
    if n_elements(folder_names) eq 2 then begin
      folderparts_1 = strsplit(folder_names[0], path_sep(), /extract)
      folderparts_2 = strsplit(folder_names[1], path_sep(), /extract)
      match_test = strcmp(folderparts_1, folderparts_2)
      wh_diff = where(match_test eq 0, count_diff, complement = wh_same, ncomplement = count_same)
      
      if count_diff eq 0 then begin
        ;; folders are the same
        if obs_info.obs_names[0] eq obs_info.obs_names[1] then diff_note = obs_info.fhd_types[0] $
        else diff_note = obs_info.fhd_types[0] + ' ' + obs_info.obs_names[0] + '-' + obs_info.obs_names[1]
        diff_save_path = folder_names[0] + path_sep()
      endif else begin
        joint_path = strjoin(folderparts_1[wh_same], path_sep())
        if strmid(folder_names[0], 0,1) eq path_sep() then joint_path = path_sep() + joint_path
        
        
        fnameparts_1 = strsplit(file_basename(folder_names[0]), '_', /extract, count = nfileparts_1)
        fnameparts_2 = strsplit(file_basename(folder_names[1]), '_', /extract, count = nfileparts_2)
        if nfileparts_1 ne nfileparts_2 then begin
          if nfileparts_1 gt nfileparts_2 then fnameparts_2 = [fnameparts_2, strarr(nfileparts_1-nfileparts_2)] $
          else fnameparts_1 = [fnameparts_1, strarr(nfileparts_2-nfileparts_1)]
        endif
        match_name_test = strcmp(fnameparts_1, fnameparts_2)
        wh_name_diff = where(match_name_test eq 0, count_name_diff, complement = wh_name_same, ncomplement = count_name_same)
        
        if count_name_diff eq 0 then begin
          diff_dir = file_basename(folder_names[0]) + '_diff'
          diff_note = strjoin(folderparts_1[wh_diff], path_sep()) + ' - ' + strjoin(folderparts_2[wh_diff], path_sep()) + ' ' + file_basename(folder_names[0])
        endif else begin
          if min(wh_name_diff) ge nfileparts_1 or min(wh_name_diff) ge nfileparts_2 then begin
            wh_name_diff = [max(wh_name_same), wh_name_diff]
            count_name_diff = count_name_diff + 1
            if count_name_same gt 1 then begin
              wh_name_same = wh_name_same[0:count_name_same-2]
              count_name_same = count_name_same-1
            endif else count_name_same = 0
          endif
          
          str1_diff = strjoin(fnameparts_1[wh_name_diff[where((wh_name_diff lt nfileparts_1) gt 0)]], '_')
          str2_diff = strjoin(fnameparts_2[wh_name_diff[where((wh_name_diff lt nfileparts_2) gt 0)]], '_')
          
          if count_name_same gt 0 then begin
            str_same = strjoin(fnameparts_1[wh_name_same], '_')
            diff_dir = str_same + '__' + str1_diff + '_minus_' + str2_diff
            diff_note = str_same + ' ' + str1_diff + ' - ' + str2_diff
          endif else begin
            diff_dir = str1_diff + '_minus_' + str2_diff
            diff_note = str1_diff + ' - ' + str2_diff
          endelse
        endelse
        
        diff_save_path = joint_path + path_sep() + diff_dir + path_sep()
        
        if count_name_same gt 0 then name_same_parts = strjoin(fnameparts_1[wh_name_same], '_') else name_same_parts = ''
        name_diff_parts = [strjoin(fnameparts_1[wh_name_diff], '_'), strjoin(fnameparts_2[wh_name_diff], '_')]
      endelse
      
      pos = strpos(diff_save_path, std_savepath)
      if pos ne -1 then diff_save_path_ext = strmid(diff_save_path, pos + strlen(std_savepath)) $
      else diff_save_path_ext = ''
      diff_plot_path = base_path('plots') + 'power_spectrum/fhd_data/' + diff_save_path_ext
    endif
    
    obs_info = {folder_names:folder_names, obs_names:obs_names, info_files:info_files, cube_files:cube_files, $
      fhd_types:fhd_types, integrated:integrated, plot_paths:plot_paths}
      
    if n_elements(diff_note) gt 0 then obs_info = create_struct(obs_info, 'diff_note', diff_note, 'diff_save_path', $
      diff_save_path, 'diff_plot_path', diff_plot_path)
      
    if n_elements(name_same_parts) gt 0 then obs_info = create_struct(obs_info, 'name_same_parts', name_same_parts, 'name_diff_parts',name_diff_parts)
    
  endelse
  
  return, obs_info
  
end
