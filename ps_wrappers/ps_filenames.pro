function ps_filenames, folder_names, obs_names_in, rts = rts, sim = sim, casa = casa, data_subdirs = data_subdirs, plot_paths = plot_paths, save_paths = save_paths

  n_filesets = max([n_elements(folder_names), n_elements(obs_names_in)])
  
  if n_elements(data_subdirs) eq 0 then data_subdirs = '' else begin
    ;; make sure there is no path separator at beginning and there is one at the end
    for i=0, n_elements(data_subdirs)-1 do if data_subdirs[i] ne '' then data_subdirs[i] = strjoin(strsplit(data_subdirs[i], path_sep(), /extract), path_sep()) + path_sep()
  endelse
  
  if n_filesets gt 1 then begin
    if n_elements(folder_names) eq 1 then folder_names = replicate(folder_names, n_filesets)
    if n_elements(folder_names) ne n_filesets then message, 'If both folder_names and obs_names_in are arrays, the number of elements must match'
    
    if n_elements(obs_names_in) eq 1 then obs_names_in = replicate(obs_names_in, n_filesets)
    if n_elements(obs_names_in) gt 0 and n_elements(obs_names_in) ne n_filesets then message, 'If both folder_names and obs_names_in are arrays, the number of elements must match'
    
    if n_elements(data_subdirs) eq 1 then data_subdirs = replicate(data_subdirs, n_filesets)
    if n_elements(data_subdirs) ne n_filesets then message, 'If data_subdirs is an array, the number of elements must match the max of folder_names & obs_names_in'
    
    if n_elements(save_paths) gt 0 then begin
      if n_elements(save_paths) eq 1 then save_paths = replicate(plosave_pathst_paths, n_filesets)
      if n_elements(save_paths) ne n_filesets then message, 'If save_paths is an array, the number of elements must match the max of folder_names & obs_names_in'
    endif else save_paths = folder_names + data_subdirs
    
    if n_elements(plot_paths) gt 0 then begin
      if n_elements(plot_paths) eq 1 then plot_paths = replicate(plot_paths, n_filesets)
      if n_elements(plot_paths) ne n_filesets then message, 'If plot_paths is an array, the number of elements must match the max of folder_names & obs_names_in'
    endif else plot_paths = folder_names + data_subdirs
  endif else begin
    if n_elements(save_paths) eq 0 then save_paths = folder_names + data_subdirs
    if n_elements(plot_paths) eq 0 then plot_paths = folder_names + data_subdirs
  endelse
  
  ;; make sure save & plot paths end in path_sep()
  pos = strpos(save_paths, path_sep(), /reverse_search)
  wh_nosep = where(pos+1-strlen(save_paths) lt 0, count_nosep)
  if count_nosep gt 0 then save_paths[wh_nosep] = save_paths[wh_nosep] + path_sep()
  
  pos = strpos(plot_paths, path_sep(), /reverse_search)
  wh_nosep = where(pos+1-strlen(plot_paths) lt 0, count_nosep)
  if count_nosep gt 0 then plot_paths[wh_nosep] = plot_paths[wh_nosep] + path_sep()
  
  if keyword_set(rts) then begin
    obs_names = strarr(n_filesets)
    rts_types = strarr(n_filesets)
    folder_basenames = strarr(n_filesets)
    info_files = strarr(n_filesets)
    ;integrated = intarr(n_filesets)
    
    for i=0, n_filesets-1 do begin
    
      folder_basenames[i] = file_basename(folder_names[i])
      rts_types[i] = folder_basenames[i]
      
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
      
      ;; first look for info files in save_path
      info_file = file_search(save_paths[i] + obs_names[i] + '*info*', count = n_infofile)
      if n_infofile gt 0 then begin
        if obs_names[i] eq '' then begin
          info_files[i] = info_file[0]
          obs_names[i] = stregex(file_basename(info_files[i]), '[0-9]+.[0-9]+_', /extract)
          if n_infofile gt 1 then begin
            print, 'More than 1 info files found, using first one'
            rts_types[i] = rts_types[i] + '_' + obs_names[i]
          endif
        endif else begin
          if n_infofile gt 1 then message, 'More than one info file found with given obs_name'
          info_files[i] = info_file[0]
          test_other_obsnames = file_search(save_paths[i] + '*info*', count = n_all_infofile)
          if n_all_infofile gt n_infofile then rts_types[i] = rts_types[i] + '_' + obs_names[i]
        endelse
        
      endif
      
      ;; then look for combined cube files in folder + data_subdir
      cube_file_list = file_search(folder_names[i] + '/' + data_subdirs[i] + obs_names[i] + '*_cube.idlsave', count = n_cubefiles)
      if n_cubefiles gt 0 then begin
        if obs_names[i] eq '' then begin
          obs_name_arr = stregex(file_basename(cube_file_list), '[0-9]+.[0-9]+_', /extract)
          wh_first = where(obs_name_arr eq obs_name_arr[0], count_first)
          cube_files = cube_file_list[wh_first]
          obs_names[i] = obs_name_arr[0]
          if count_first lt n_cubefiles then begin
            print, 'More than one obs_name found, using first obs_name (' + obs_name_arr[0] + ', ' + number_formatter(count_first) + ' files)'
            if n_elements(info_files) eq 0 then rts_types[i] = rts_types[i] + '_' + obs_names[i]
          endif
        endif else begin
          cube_files = cube_file_list
          
          if n_elements(info_files) eq 0 then begin
            test_other_obsnames = file_search(folder_names[i] + '/' + data_subdirs[i] + '*_cube.idlsave', count = n_all_cubefiles)
            if n_all_cubefiles gt n_cubefiles then rts_types[i] = rts_types[i] + '_' + obs_names[i]
          endif
        endelse
        
      endif
      
      ;; then look for original fits files
      fits_file_list = file_search(folder_names[i] + '/' + data_subdirs[i] + obs_names[i] + '*_image*.fits', count = n_fitsfiles)
      if n_fitsfiles gt 0 then begin
        if obs_names[i] eq '' then begin
          obs_name_arr = stregex(file_basename(fits_file_list), '[0-9]+.[0-9]+_', /extract)
          wh_first = where(obs_name_arr eq obs_name_arr[0], count_first)
          if count_first lt n_fitsfiles then $
            print, 'More than one obs_name found, using first obs_name (' + obs_name_arr[0] + ', ' + number_formatter(count_first) + ' files)'
            
          fits_files = fits_file_list[wh_first]
          obs_names[i] = obs_name_arr[0]
        endif else begin
          fits_files = fits_file_list
        endelse
        
      endif
      
      if n_elements(cube_files) eq 0 and n_elements(fits_files) and info_files[i] eq '' then message, 'No cube or info files found in folder ' + folder_names[i]
      
      if n_elements(fits_files) eq 0 then begin
        fits_file_list = ''
        weightfile_list = ''
        variancefile_list = ''
      endif else begin
        ;; now get weights & variance files
        weightfile_list = file_search(folder_names[i] + '/' + data_subdirs[i] + obs_names[i] + '*_weights*.fits', count = n_wtfiles)
        if n_wtfiles ne n_elements(fits_files) and info_files[i] eq '' then message, 'number of weight files does not match number of datafiles'
        
        variancefile_list = file_search(folder_names[i] + '/' + data_subdirs[i] + obs_names[i] + '*_weights*.fits', count = n_varfiles)
        if n_varfiles ne n_elements(fits_files) and info_files[i] eq '' then message, 'number of variance files does not match number of datafiles'
      endelse
      
      if i eq 0 then begin
        tag = 'fs' + number_formatter(i)
        if n_cubefiles gt 0 then cubefiles = create_struct(tag, cube_files) else cubefiles = create_struct(tag, '')
        if n_fitsfiles gt 0 then datafiles = create_struct(tag, fits_files) else datafiles = create_struct(tag, '')
        weightfiles = create_struct(tag, weightfile_list)
        variancefiles = create_struct(tag, variancefile_list)
      endif else begin
        tag = 'fs' + number_formatter(i)
        if n_cubefiles gt 0 then cubefiles = create_struct(cubefiles, tag, cube_files) else cubefiles = create_struct(cubefiles, tag, '')
        if n_fitsfiles gt 0 then datafiles = create_struct(datafiles, tag, fits_files) else datafiles = create_struct(datafiles, tag, '')
        weightfiles = create_struct(weightfiles, tag, weightfile_list)
        variancefiles = create_struct(variancefiles, tag, variancefile_list)
      endelse
      undefine, fits_files, weightfile_list, variancefile_list, cube_files
      
    endfor
    
    obs_info = {folder_names:folder_names, folder_basenames:folder_basenames, obs_names:obs_names, info_files:info_files, cube_files:cubefiles, $
      datafiles:datafiles, weightfiles:weightfiles, variancefiles:variancefiles, rts_types:rts_types, plot_paths:plot_paths, save_paths:save_paths}
      
  endif else if keyword_set(casa) then begin
    obs_names = strarr(n_filesets)
    casa_types = strarr(n_filesets)
    info_files = strarr(n_filesets)
    ;integrated = intarr(n_filesets)
    
    for i=0, n_filesets-1 do begin
      ;; check for folder existence, otherwise look for common folder names to figure out full path. If none found, try base_path('data') + 'mit_data/'
    
      casa_types[i] = file_basename(folder_names[i])
      
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
      
      ;; first look for info files in save_paths
      info_file = file_search(save_paths[i] + obs_names[i] + '*info*', count = n_infofile)
      if n_infofile gt 0 then begin
        if obs_names[i] eq '' then begin
          if n_infofile gt 1 then print, 'More than 1 info files found, using first one'
          info_files[i] = info_file[0]
        ;obs_names[i] = stregex(info_files[i], '[0-9]+.[0-9]+_', /extract)
        endif else begin
          if n_infofile gt 1 then message, 'More than one info file found with given obs_name'
          info_files[i] = info_file[0]
        endelse
        
      endif
      
      ;; then look for cube files in folder + data_subdir
      cube_file_list = file_search(folder_names[i] + '/' + data_subdirs[i] + obs_names[i] + '*_holo_[xy]*.fits', count = n_cubefiles)
      if n_cubefiles gt 0 then begin
        ;        if obs_names[i] eq '' then begin
        ;          obs_name_arr = stregex(cube_file_list, '[0-9]+.[0-9]+_', /extract)
        ;          wh_first = where(obs_name_arr eq obs_name_arr[0], count_first)
        ;          if count_first lt n_cubefiles then $
        ;            print, 'More than one obs_name found, using first obs_name (' + obs_name_arr[0] + ', ' + number_formatter(count_first) + ' files)'
        ;          datafiles = cube_file_list[wh_first]
        ;          obs_names[i] = obs_name_arr[0]
        ;        endif else begin
        datafiles = cube_file_list
      ;        endelse
        
      endif
      
      ;; now get weights & variance files
      weightfile_list = file_search(folder_names[i] + '/' + data_subdirs[i] + obs_names[i] + '*_psf_*.fits', count = n_wtfiles)
      if n_wtfiles ne n_elements(datafiles) and info_files[i] eq '' then message, 'number of weight files does not match number of datafiles'
      
      variancefile_list = file_search(folder_names[i] + '/' + data_subdirs[i] + obs_names[i] + '*psfbeamsquare*.fits', count = n_varfiles)
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
      weightfiles:weightfiles, variancefiles:variancefiles, casa_types:casa_types, plot_paths:plot_paths, save_paths:save_paths}
      
  endif else begin
    obs_names = strarr(n_filesets)
    fhd_types = strarr(n_filesets)
    folder_basenames = strarr(n_filesets)
    info_files = strarr(n_filesets)
    integrated = intarr(n_filesets)
    
    for i=0, n_filesets-1 do begin
    
      folder_basenames[i] = file_basename(folder_names[i])
      fhd_types[i] = folder_basenames[i]
      
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
      
      ;; first look for integrated info files in save_paths with names like Combined_obs_...
      info_file = file_search(save_paths[i] +  'Combined_obs_' + obs_names[i] + '*info*', count = n_infofile)
      if n_infofile gt 0 then begin
        if obs_names[i] eq '' then begin
          info_files[i] = info_file[0]
          if stregex(info_files[i], '[0-9]+-[0-9]+', /boolean) then obs_names[i] = stregex(info_files[i], '[0-9]+-[0-9]+', /extract) else begin
            start_pos = strpos(info_files[i], 'Combined_obs_') + strlen('Combined_obs_')
            end_pos = strpos(strmid(info_files[i], start_pos), '_cube')
            obs_names[i] = strmid(info_files[i], start_pos, end_pos)
          endelse
          if n_infofile gt 1 then begin
            print, 'More than 1 info files found, using first one'
            fhd_types[i] = fhd_types[i] + '_' + obs_names[i]
          endif
        endif else begin
          if n_infofile gt 1 then message, 'More than one info file found with given obs_name'
          info_files[i] = info_file
          if stregex(info_files[i], '[0-9]+-[0-9]+', /boolean) then obs_names[i] = stregex(info_files[i], '[0-9]+-[0-9]+', /extract) else begin
            start_pos = strpos(info_files[i], 'Combined_obs_') + strlen('Combined_obs_')
            end_pos = strpos(strmid(info_files[i], start_pos), '_cube')
            obs_names[i] = strmid(info_files[i], start_pos, end_pos)
          endelse
          test_other_obsnames = file_search(save_paths[i] +  'Combined_obs_' + '*info*', count = n_all_infofile)
          if n_all_infofile gt n_infofile then fhd_types[i] = fhd_types[i] + '_' + obs_names[i]
        endelse
        integrated[i]=1
        
      endif else begin
        ;; then look for single obs info files
        info_file = file_search(save_paths[i] + obs_name_single + '*info*', count = n_infofile)
        if n_infofile gt 0 then begin
          info_basename = file_basename(info_file)
          if obs_names[i] eq '' then begin
            info_files[i] = info_file[0]
            if stregex(info_basename[0], '[0-9]+', /boolean) then obs_names[i] = stregex(info_basename[0], '[0-9]+', /extract) else begin
              end_pos = strpos(info_basename[0], '_cube')
              obs_names[i] = strmid(info_basename[0], 0, end_pos)
            endelse
            if n_infofile gt 1 then begin
              print, 'More than 1 info files found, using first one'
              fhd_types[i] = fhd_types[i] + '_' + obs_names[i]
            endif
          endif else begin
            if n_infofile gt 1 then message, 'More than one info file found with given obs_name'
            info_files[i] = info_file
            if stregex(info_basename[0], '[0-9]+', /boolean) then obs_names[i] = stregex(info_basename[0], '[0-9]+', /extract) else begin
              end_pos = strpos(info_basename[0], '_cube')
              obs_names[i] = strmid(info_basename[0], 0, end_pos)
            endelse
            test_other_obsnames = file_search(save_paths[i] + '*info*', count = n_all_infofile)
            if n_all_infofile gt n_infofile then fhd_types[i] = fhd_types[i] + '_' + obs_names[i]
          endelse
          integrated[i]=0
          
        endif
      endelse
      
      ;; first look for integrated cube files in folder + data_dir with names like Combined_obs_...
      cube_file_list = file_search(folder_names[i] + '/' + data_subdirs[i] + 'Combined_obs_' + obs_names[i] + '*_cube*.sav', count = n_cubefiles)
      if n_cubefiles gt 0 then begin
        cube_basename = file_basename(cube_file_list)
        pol_exist = stregex(cube_basename, '[xy][xy]', /boolean, /fold_case)
        
        if obs_names[i] eq '' then begin
          obs_name_arr = strarr(n_cubefiles)
          for j=0, n_cubefiles-1 do begin
            start_pos = strpos(cube_basename[j], 'Combined_obs_') + strlen('Combined_obs_')
            end_pos_even = strpos(strmid(cube_basename[j], start_pos), '_even')
            end_pos_odd = strpos(strmid(cube_basename[j], start_pos), '_odd')
            end_pos_cube = strpos(strmid(cube_basename[j], start_pos), '_cube') ;; always > -1
            end_pos = end_pos_even > end_pos_odd
            wh_noend = where(end_pos eq -1, count_noend)
            if count_noend gt 0 then end_pos[wh_noend] = end_pos_cube[wh_noend]
            
            ;obs_name_arr = stregex(cube_basename, '[0-9]+-[0-9]+', /extract)
            obs_name_arr[j] = strmid(cube_basename[j], start_pos, end_pos)
          endfor
          
          wh_first = where(obs_name_arr eq obs_name_arr[0], count_first)
          if max(pol_exist[wh_first]) gt 0 then begin
            if min(pol_exist[wh_first]) eq 0 then message, 'some files with first obs_name have pol identifiers and some do not'
            pols = stregex(cube_basename[wh_first], '[xy][xy]', /extract, /fold_case)
            
            pols_inc = pols[0]
            pol_num = intarr(count_first)
            for pol_i=0, count_first-1 do begin
              wh_pol = where(pols_inc eq pols[pol_i], count_pol)
              if count_pol eq 1 then pol_num[pol_i] = wh_pol[0] else begin
                pols_inc = [pols_inc, pols[pol_i]]
                pol_num[pol_i] = n_elements(pols_inc)-1
              endelse
            endfor
            
            for pol_i=0, n_elements(pols_inc) do begin
              wh_pol = where(pol_num eq pol_i, count_pol)
              if count_pol gt 2 then message, 'More than two cubes found with first obs_name and the same polarization'
            endfor
          endif else if count_first gt 2 then message, 'More than two cubes found with first obs_name'
          
          datafile = cube_file_list[wh_first]
          obs_names[i] = obs_name_arr[0]
          if count_first lt n_elements(cube_basename) then begin
            print, 'More than one obs_name found, using first obs_name (' + obs_name_arr[0] + ', ' + number_formatter(count_first) + ' files)'
            fhd_types[i] = fhd_types[i] + '_' + obs_names[i]
          endif
        endif else begin
          if max(pol_exist) gt 0 then begin
            if min(pol_exist) eq 0 then message, 'some files with given obs_name have pol identifiers and some do not'
            pols = stregex(cube_basename, '[xy][xy]', /extract, /fold_case)
            
            pols_inc = pols[0]
            pol_num = intarr(n_cubefiles)
            for pol_i=0, n_cubefiles-1 do begin
              wh_pol = where(pols_inc eq pols[pol_i], count_pol)
              if count_pol eq 1 then pol_num[pol_i] = wh_pol[0] else begin
                pols_inc = [pols_inc, pols[pol_i]]
                pol_num[pol_i] = n_elements(pols_inc)-1
              endelse
            endfor
            
            for pol_i=0, n_elements(pols_inc) do begin
              wh_pol = where(pol_num eq pol_i, count_pol)
              if count_pol gt 2 then message, 'More than two cubes found with given obs_name and the same polarization'
            endfor
          endif else if n_cubefiles gt 2 then message, 'More than two cubes found with given obs_name'
          
          datafile = cube_file_list
          test_other_obsnames = file_search(folder_names[i] + '/' + data_subdirs[i] + 'Combined_obs_' + '*_cube*.sav', count = n_all_infofile)
          if n_all_infofile gt n_infofile then fhd_types[i] = fhd_types[i] + '_' + obs_names[i]
        endelse
        integrated[i]=1
        
      endif else begin
        ;; then look for single obs cube files
        cube_file_list = file_search(folder_names[i] + '/' + data_subdirs[i] + obs_names[i] + '*_cube*.sav', count = n_cubefiles)
        if n_cubefiles gt 0 then begin
          cube_basename = file_basename(cube_file_list)
          pol_exist = stregex(cube_basename, '[xy][xy]', /boolean, /fold_case)
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
            if max(pol_exist[wh_first]) gt 0 then begin
              if min(pol_exist[wh_first]) eq 0 then message, 'some files with first obs_name have pol identifiers and some do not'
              pols = stregex(cube_basename[wh_first], '[xy][xy]', /extract, /fold_case)
              
              pols_inc = pols[0]
              pol_num = intarr(count_first)
              for pol_i=0, count_first-1 do begin
                wh_pol = where(pols_inc eq pols[pol_i], count_pol)
                if count_pol eq 1 then pol_num[pol_i] = wh_pol[0] else begin
                  pols_inc = [pols_inc, pols[pol_i]]
                  pol_num[pol_i] = n_elements(pols_inc)-1
                endelse
              endfor
              
              for pol_i=0, n_elements(pols_inc) do begin
                wh_pol = where(pol_num eq pol_i, count_pol)
                if count_pol gt 2 then message, 'More than two cubes found with first obs_name and the same polarization'
              endfor
            endif else if count_first gt 2 then message, 'More than two cubes found with first obs_name'
            datafile = cube_file_list[wh_first]
            obs_names[i] = obs_name_arr[0]
            integrated[i]=0
            
            if count_first lt n_elements(cube_file_list) then begin
              print, 'More than one obs_name found, using first obs_name (' + obs_name_arr[0] + ', ' + number_formatter(count_first) + ' files)'
              fhd_types[i] = fhd_types[i] + '_' + obs_names[i]
            endif
            
          endif else begin
            if max(pol_exist) gt 0 then begin
              if min(pol_exist) eq 0 then message, 'some files with given obs_name have pol identifiers and some do not'
              pols = stregex(cube_basename, '[xy][xy]', /extract, /fold_case)
              
              pols_inc = pols[0]
              pol_num = intarr(n_cubefiles)
              for pol_i=0, n_cubefiles-1 do begin
                wh_pol = where(pols_inc eq pols[pol_i], count_pol)
                if count_pol eq 1 then pol_num[pol_i] = wh_pol[0] else begin
                  pols_inc = [pols_inc, pols[pol_i]]
                  pol_num[pol_i] = n_elements(pols_inc)-1
                endelse
              endfor
              
              for pol_i=0, n_elements(pols_inc) do begin
                wh_pol = where(pol_num eq pol_i, count_pol)
                if count_pol gt 2 then message, 'More than two cubes found with given obs_name and the same polarization'
              endfor
            endif else if n_cubefiles gt 2 then message, 'More than two cubes found with given obs_name'
            datafile = cube_file_list
            
            test_other_obsnames = file_search(folder_names[i] + '/' + data_subdirs[i] + '*_cube*.sav', count = n_all_infofile)
            if n_all_infofile gt n_infofile then fhd_types[i] = fhd_types[i] + '_' + obs_names[i]
          endelse
        endif
        
      endelse
      
      if n_elements(datafile) eq 0 then begin
        ;; finally look for uniformly gridded uvf files
        ;; Combined
        cube_file_list = file_search(folder_names[i] + '/' + data_subdirs[i] + 'Combined_obs_' + obs_names[i] + '*_uvf.sav', count = n_cubefiles)
        if n_cubefiles gt 0 then begin
          cube_basename = file_basename(cube_file_list)
          pol_exist = stregex(cube_basename, '[xy][xy]', /boolean, /fold_case)
          
          if obs_names[i] eq '' then begin
            obs_name_arr = strarr(n_cubefiles)
            for j=0, n_cubefiles-1 do begin
              start_pos = strpos(cube_basename[j], 'Combined_obs_') + strlen('Combined_obs_')
              end_pos_even = strpos(strmid(cube_basename[j], start_pos), '_even')
              end_pos_odd = strpos(strmid(cube_basename[j], start_pos), '_odd')
              end_pos_cube = strpos(strmid(cube_basename[j], start_pos), '_cube') ;; always > -1
              end_pos = end_pos_even > end_pos_odd
              wh_noend = where(end_pos eq -1, count_noend)
              if count_noend gt 0 then end_pos[wh_noend] = end_pos_cube[wh_noend]
              
              ;obs_name_arr = stregex(cube_basename, '[0-9]+-[0-9]+', /extract)
              obs_name_arr[j] = strmid(cube_basename[j], start_pos, end_pos)
            endfor
            
            wh_first = where(obs_name_arr eq obs_name_arr[0], count_first)
            if max(pol_exist[wh_first]) gt 0 then begin
              if min(pol_exist[wh_first]) eq 0 then message, 'some files with first obs_name have pol identifiers and some do not'
              pols = stregex(cube_basename[wh_first], '[xy][xy]', /extract, /fold_case)
              
              pols_inc = pols[0]
              pol_num = intarr(count_first)
              for pol_i=0, count_first-1 do begin
                wh_pol = where(pols_inc eq pols[pol_i], count_pol)
                if count_pol eq 1 then pol_num[pol_i] = wh_pol[0] else begin
                  pols_inc = [pols_inc, pols[pol_i]]
                  pol_num[pol_i] = n_elements(pols_inc)-1
                endelse
              endfor
              
              for pol_i=0, n_elements(pols_inc) do begin
                wh_pol = where(pol_num eq pol_i, count_pol)
                if count_pol gt 2 then message, 'More than two cubes found with first obs_name and the same polarization'
              endfor
            endif else if count_first gt 2 then message, 'More than two cubes found with first obs_name'
            
            if count_first lt n_elements(cube_basename) then $
              print, 'More than one obs_name found, using first obs_name (' + obs_name_arr[0] + ', ' + number_formatter(count_first) + ' files)'
              
            datafile = cube_file_list[wh_first]
            obs_names[i] = obs_name_arr[0]
          endif else begin
            if max(pol_exist) gt 0 then begin
              if min(pol_exist) eq 0 then message, 'some files with given obs_name have pol identifiers and some do not'
              pols = stregex(cube_basename, '[xy][xy]', /extract, /fold_case)
              
              pols_inc = pols[0]
              pol_num = intarr(n_cubefiles)
              for pol_i=0, n_cubefiles-1 do begin
                wh_pol = where(pols_inc eq pols[pol_i], count_pol)
                if count_pol eq 1 then pol_num[pol_i] = wh_pol[0] else begin
                  pols_inc = [pols_inc, pols[pol_i]]
                  pol_num[pol_i] = n_elements(pols_inc)-1
                endelse
              endfor
              
              for pol_i=0, n_elements(pols_inc) do begin
                wh_pol = where(pol_num eq pol_i, count_pol)
                if count_pol gt 2 then message, 'More than two cubes found with given obs_name and the same polarization'
              endfor
            endif else if n_cubefiles gt 2 then message, 'More than two cubes found with given obs_name'
            
            datafile = cube_file_list
          endelse
          integrated[i]=1
          uvf_input=1
          
          ;; look for beam files
          beam_file_list = file_search(folder_names[i] + '/' + data_subdirs[i] + 'Combined_obs_' + obs_names[i] + '*gridded*beam2*image*.sav', count = n_beamfiles)
          if n_beamfiles gt 0 then begin
            if n_elements(datafile) gt 1 then begin
              if n_beamfiles eq 1 then beamfiles = strarr(n_elements(datafile)) + beam_file_list $
              else if n_beamfiles ne n_elements(datafile) then stop
            endif else if n_beamfiles ne n_elements(datafile) then stop
          endif
          
        endif else begin
          ;; single
          cube_file_list = file_search(folder_names[i] + '/' + data_subdirs[i] + obs_names[i] + '*_uvf.sav', count = n_cubefiles)
          if n_cubefiles gt 0 then begin
            cube_basename = file_basename(cube_file_list)
            pol_exist = stregex(cube_basename, '[xy][xy]', /boolean, /fold_case)
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
              if max(pol_exist[wh_first]) gt 0 then begin
                if min(pol_exist[wh_first]) eq 0 then message, 'some files with first obs_name have pol identifiers and some do not'
                pols = stregex(cube_basename[wh_first], '[xy][xy]', /extract, /fold_case)
                
                pols_inc = pols[0]
                pol_num = intarr(count_first)
                for pol_i=0, count_first-1 do begin
                  wh_pol = where(pols_inc eq pols[pol_i], count_pol)
                  if count_pol eq 1 then pol_num[pol_i] = wh_pol[0] else begin
                    pols_inc = [pols_inc, pols[pol_i]]
                    pol_num[pol_i] = n_elements(pols_inc)-1
                  endelse
                endfor
                
                for pol_i=0, n_elements(pols_inc) do begin
                  wh_pol = where(pol_num eq pol_i, count_pol)
                  if count_pol gt 2 then message, 'More than two cubes found with first obs_name and the same polarization'
                endfor
              endif else if count_first gt 2 then message, 'More than two cubes found with first obs_name'
              
              if count_first lt n_elements(cube_file_list) then $
                print, 'More than one obs_name found, using first obs_name (' + obs_name_arr[0] + ', ' + number_formatter(count_first) + ' files)'
                
              datafile = cube_file_list[wh_first]
              obs_names[i] = obs_name_arr[0]
              integrated[i]=0
            endif else begin
              if max(pol_exist) gt 0 then begin
                if min(pol_exist) eq 0 then message, 'some files with given obs_name have pol identifiers and some do not'
                pols = stregex(cube_basename, '[xy][xy]', /extract, /fold_case)
                
                pols_inc = pols[0]
                pol_num = intarr(n_cubefiles)
                for pol_i=0, n_cubefiles-1 do begin
                  wh_pol = where(pols_inc eq pols[pol_i], count_pol)
                  if count_pol eq 1 then pol_num[pol_i] = wh_pol[0] else begin
                    pols_inc = [pols_inc, pols[pol_i]]
                    pol_num[pol_i] = n_elements(pols_inc)-1
                  endelse
                endfor
                
                for pol_i=0, n_elements(pols_inc) do begin
                  wh_pol = where(pol_num eq pol_i, count_pol)
                  if count_pol gt 2 then message, 'More than two cubes found with given obs_name and the same polarization'
                endfor
              endif else if n_cubefiles gt 2 then message, 'More than two cubes found with given obs_name'

              datafile = cube_file_list
            endelse
          endif
          uvf_input=1
          
          ;; look for beam files
          beam_file_list = file_search(folder_names[i] + '/' + data_subdirs[i] + obs_names[i] + '*gridded*beam2*image*.sav', count = n_beamfiles)
          if n_beamfiles gt 0 then begin
            if n_elements(datafile) gt 1 then begin
              if n_beamfiles eq 1 then beamfiles = strarr(n_elements(datafile)) + beam_file_list[0] $
              else if n_beamfiles ne n_elements(datafile) then stop else beamfiles = beam_file_list
            endif else if n_beamfiles ne n_elements(datafile) then stop else beamfiles = beam_file_list
          endif
          
        endelse
      endif
      
      if n_elements(datafile) eq 0 and info_files[i] eq '' then message, 'No cube or info files found in folder ' + folder_names[i]
      
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
        else diff_note = obs_info.folder_names[0] + ' ' + obs_info.obs_names[0] + '-' + obs_info.obs_names[1]
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
        
        fhdtypeparts_1 = strsplit(file_basename(fhd_types[0]), '_', /extract, count = nfhdtypeparts_1)
        fhdtypeparts_2 = strsplit(file_basename(fhd_types[1]), '_', /extract, count = nfhdtypeparts_2)
        if nfhdtypeparts_1 ne nfhdtypeparts_2 then begin
          if nfhdtypeparts_1 gt nfhdtypeparts_2 then fhdtypeparts_2 = [fhdtypeparts_2, strarr(nfhdtypeparts_1-nfhdtypeparts_2)] $
          else fhdtypeparts_1 = [fhdtypeparts_1, strarr(nfhdtypeparts_2-nfhdtypeparts_1)]
        endif
        match_fhdtype_test = strcmp(fhdtypeparts_1, fhdtypeparts_2)
        wh_fhdtype_diff = where(match_fhdtype_test eq 0, count_fhdtype_diff, complement = wh_fhdtype_same, ncomplement = count_fhdtype_same)
        
        if count_name_diff eq 0 then begin
          ;; same folder name, different directories
          diff_dir = file_basename(folder_names[0]) + '_diff'
          if count_fhdtype_diff eq 0 then diff_note = strjoin(folderparts_1[wh_diff], path_sep()) + ' - ' + strjoin(folderparts_2[wh_diff], path_sep()) + ' ' + file_basename(folder_names[0]) $
          else $
            diff_note = strjoin(folderparts_1[wh_diff], path_sep()) + strjoin(fhdtypeparts_1[wh_fhdtype_diff], '_') + ' - ' + $
            strjoin(folderparts_2[wh_diff], path_sep()) + strjoin(fhdtypeparts_2[wh_fhdtype_diff], '_') + ' ' + file_basename(folder_names[0])
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
        if count_name_diff eq 0 then name_diff_parts = strarr(2) else name_diff_parts = [str1_diff, str2_diff]
        
        if min(wh_fhdtype_diff) ge nfhdtypeparts_1 or min(wh_fhdtype_diff) ge nfhdtypeparts_2 then begin
          wh_fhdtype_diff = [max(wh_fhdtype_same), wh_fhdtype_diff]
          count_fhdtype_diff = count_fhdtype_diff + 1
          if count_fhdtype_same gt 1 then begin
            wh_fhdtype_same = wh_fhdtype_same[0:count_fhdtype_same-2]
            count_fhdtype_same = count_fhdtype_same-1
          endif else count_fhdtype_same = 0
        endif
        
        str1_diff_fhdtype = strjoin(fhdtypeparts_1[wh_fhdtype_diff[where((wh_fhdtype_diff lt nfhdtypeparts_1) gt 0)]], '_')
        str2_diff_fhdtype = strjoin(fhdtypeparts_2[wh_fhdtype_diff[where((wh_fhdtype_diff lt nfhdtypeparts_2) gt 0)]], '_')
        
        if count_fhdtype_same gt 0 then fhdtype_same_parts = strjoin(fhdtypeparts_1[wh_fhdtype_same], '_') else fhdtype_same_parts = ''
        if count_fhdtype_diff eq 0 then fhdtype_diff_parts = strarr(2) else fhdtype_diff_parts = [str1_diff_fhdtype, str2_diff_fhdtype]
        
      endelse
      
    endif
    
    obs_info = {folder_names:folder_names, folder_basenames:folder_basenames, obs_names:obs_names, info_files:info_files, cube_files:cube_files, $
      fhd_types:fhd_types, integrated:integrated, plot_paths:plot_paths, save_paths:save_paths}
      
    if n_elements(diff_note) gt 0 then obs_info = create_struct(obs_info, 'diff_note', diff_note, 'diff_save_path', $
      diff_save_path)
      
    if n_elements(name_same_parts) gt 0 then obs_info = create_struct(obs_info, 'name_same_parts', name_same_parts, 'name_diff_parts',name_diff_parts)
    
    if n_elements(fhdtype_same_parts) gt 0 then obs_info = create_struct(obs_info, 'fhdtype_same_parts', fhdtype_same_parts, 'fhdtype_diff_parts',fhdtype_diff_parts)
    
    
    if n_elements(uvf_input) gt 0 then obs_info = create_struct(obs_info, 'uvf_input', uvf_input)
    
    if n_elements(beamfiles) gt 0 then obs_info = create_struct(obs_info, 'beam_files', beamfiles)
    
  endelse
  
  return, obs_info
  
end
