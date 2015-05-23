; Wrapper to make movies for FHD output image cubes.
; Movies can run over frequencies, obsids or a set of supplied files.
; This wrapper also supports making movies over differenced or ratioed images
;    (anything supported by ps_core/cube_images)
; To make a movie over frequency, simply specify 1-2 folder_names and
;    optionally 1-2 obs_names_in, cube_types, pols & even/odd designations
;    (specifying 2 sets results in differencing by default but ratios and difference ratios are optional)
;   by default, movies over frequency will make one frame per frequency,
;   for averaging use the n_freq_frames keyword (number of channels averaged per frame = n_freq/n_freq_frames)
; To make a movie over all obsids in a given folder, specify the folder and set keyword movie_axis = 'obsid'
;   As above, cube_types, pols & even/odd designations are available
; To make a movie over a set of supplied files, supply a list of folder_names &/or obs_names_in
;
; The individual frame images are deleted by default, to prevent this set delete_images = 0
;
; The movie file name can be defined using the movie_file keyword,
;   otherwise a default choice is made and the filename used is reported
;
;
pro mit_cube_movie, folder_names, obs_names_in, exact_obsnames = exact_obsnames, data_subdirs=data_subdirs, $
    movie_axis = movie_axis, cube_types = cube_types, pols = pols, evenodd = evenodd, $
    rts = rts, sim = sim, casa = casa, png = png, slice_range = slice_range, sr2 = sr2, $
    nvis_norm = nvis_norm, ratio = ratio, diff_ratio = diff_ratio, diff_frac = diff_frac, $
    log = log, data_range = data_range, color_profile = color_profile, sym_color = sym_color, $
    window_num = window_num, plot_as_map = plot_as_map, $
    delete_images = delete_images, n_freq_frames = n_freq_frames, movie_file = movie_file
    
  if n_elements(evenodd) eq 0 then evenodd = 'even'
  
  if n_elements(movie_axis) eq 0 then begin
    if max([n_elements(folder_names), n_elements(obs_names_in)]) gt 2 then movie_axis = 'file' $
    else movie_axis = 'frequency'
  endif
  
  if n_elements(delete_images) eq 0 then delete_images = 1
  
  movie_axis_enum = ['file', 'frequency', 'obsid']
  where_axis = where(movie_axis_enum eq movie_axis, count)
  if count eq 0 then message, 'movie_axis not recognized'
  
  if movie_axis eq 'frequency' then begin
    if n_elements(folder_names) gt 2 then message, 'No more than 2 folder_names can be supplied'
    if n_elements(evenodd) gt 2 then message, 'No more than 2 evenodd values can be supplied'
    if n_elements(obs_names_in) gt 2 then message, 'No more than 2 obs_names can be supplied'
  endif
  
  if movie_axis eq 'obsid' then begin
    if n_elements(folder_names) gt 2 then message, 'No more than 2 folder_names can be supplied'
    if n_elements(evenodd) gt 2 then message, 'No more than 2 evenodd values can be supplied'
  endif
  
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
  
  if n_elements(data_subdirs) eq 0 then data_subdirs = 'Healpix/'
  
  if movie_axis eq 'obsid' and n_elements(obs_names_in) eq 0 then begin
    for i=0, n_elements(folder_names)-1 do begin
      ;; look for single obs cube files to get all single obs names
      cube_file_list = file_search(folder_names[i] + '/' + data_subdirs + '*_cube*.sav', count = n_cubefiles)
      if n_cubefiles gt 0 then begin
        ;;first find single obs
        cube_basename = file_basename(cube_file_list)
        single_obs_test = stregex(cube_basename, '^[0-9]+[^-]*', /boolean)
        if total(single_obs_test) gt 0 then begin
          cube_file_list = cube_file_list[where(single_obs_test)]
          cube_basename = cube_basename[where(single_obs_test)]
          n_cubefiles = n_elements(cube_file_list)
        endif else message, 'no single_obs cube files found'
        
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
        obs_names_use = obs_name_arr[uniq(obs_name_arr, sort(obs_name_arr))]
        
      endif else message, 'no single_obs cube files found'
      
      if i eq 0 then obs_names = obs_names_use else stop
    endfor
  endif
  
  if movie_axis eq 'obsid' or movie_axis eq 'file' then begin
    for i=0, n_elements(obs_names)-1 do begin
      if movie_axis eq 'obsid' then begin
        folder_use = folder_names
        obs_use = obs_names[i]
        exact_obsnames = 1
      endif else begin
        if n_elements(folder_names) gt 1 then folder_use = folder_names[i] else folder_use = folder_names
        if n_elements(obs_names_in) gt 1 then obs_use = obs_names_in[i] else $
          if n_elements(obs_names_in) gt 0 then obs_use = obs_names_in
      endelse
      
      save_paths = folder_use + '/ps/'
      plot_paths = save_paths + 'plots/'
      this_obs_info = ps_filenames(folder_use, obs_use, exact_obsnames = exact_obsnames, $
        rts = rts, sim = sim, casa = casa, $
        data_subdirs = data_subdirs, save_paths = save_paths, plot_paths = plot_paths)
        
      if i eq 0 then obs_info = [this_obs_info] else obs_info = [obs_info, this_obs_info]
    endfor
    
  endif else begin
    save_paths = folder_names + '/ps/'
    plot_paths = save_paths + 'plots/'
    
    obs_info = ps_filenames(folder_names, obs_names_in, exact_obsnames = exact_obsnames, $
      rts = rts, sim = sim, casa = casa, $
      data_subdirs = data_subdirs, save_paths = save_paths, plot_paths = plot_paths)
  endelse
  
  if movie_axis eq 'frequency'then begin
  
    if obs_info.integrated[0] eq 1 then obs_varname = 'obs_arr' else obs_varname = 'obs'
    obs_arr = getvar_savefile(obs_info.cube_files.(0)[0], obs_varname)
    n_freq_orig = obs_arr[0].n_freq
    n_avg = getvar_savefile(obs_info.cube_files.(0)[0], 'n_avg')
    n_freq = n_freq_orig / n_avg
    undefine_fhd, obs_arr
    
    if n_elements(n_freq_frames) gt 0 then begin
      n_freq_per_frame = floor(n_freq / n_freq_frames)
      n_frames = n_freq_frames
    endif else begin
      n_frames = n_freq
      n_freq_per_frame = 1
    endelse
    
    plotfiles = strarr(n_frames)
    
    for i=0, n_frames-1 do begin
      if n_freq_per_frame eq 1 then slice_range = i else $
        slice_range = n_freq_per_frame*i + [0,n_freq_per_frame-1]
        
      cube_images, folder_names, obs_info, nvis_norm = nvis_norm, pols = pols, cube_types = cube_types, evenodd = evenodd, $
        png = png, slice_range = slice_range, sr2 = sr2, $
        ratio = ratio, diff_ratio = diff_ratio, diff_frac = diff_frac, $
        log = log, data_range = data_range, color_profile = color_profile, sym_color = sym_color, $
        window_num = window_num, plot_as_map = plot_as_map, plotfile = temp
        
      plotfiles[i] = temp
      undefine, temp
    endfor
    
  endif else begin
  
    n_frames = n_elements(obs_info)
    plotfiles = strarr(n_frames)
    if n_elements(slice_range) gt 0 then slice_range_use = slice_range
    if n_elements(sr2) gt 0 then sr2_use = sr2
    
    for i=0, n_frames-1 do begin
      if movie_axis eq 'file' and n_elements(folder_names) gt 1 then $
        folder_use = folder_names[i] else folder_use = folder_names
        
      cube_images, folder_use, obs_info[i], nvis_norm = nvis_norm, pols = pols, cube_types = cube_types, evenodd = evenodd, $
        png = png, slice_range = slice_range_use, sr2 = sr2_use, $
        ratio = ratio, diff_ratio = diff_ratio, diff_frac = diff_frac, $
        log = log, data_range = data_range, color_profile = color_profile, sym_color = sym_color, $
        window_num = window_num, plot_as_map = plot_as_map, plotfile = temp
        
      plotfiles[i] = temp
      undefine, temp
      
      if n_elements(slice_range) eq 0 then undefine, slice_range_use
    endfor
  endelse
  
  if keyword_set(png) then begin
    have_exten = stregex(plotfiles, '.png' ,/boolean)
    if min(have_exten) eq 0 then plotfiles[where(have_exten eq 0)] += '.png'
    
    if n_elements(movie_file) eq 0 then begin
      case movie_axis of
        'frequency':  begin
          ch_loc = stregex(plotfiles[0], '_ch[0-9]+')
          movie_file = strmid(plotfiles[0], 0, ch_loc)
          if n_freq_per_frame gt 1 then movie_file = movie_file + '_' + number_formatter(n_freq_per_frame) + 'chave'
          movie_file = movie_file + '_freq_movie.gif'
        end
        'obsid': begin
          temp = stregex(plotfiles[0], obs_info[0].obs_names, length = temp2)
          movie_file = strmid(plotfiles[0], 0, temp) + obs_info[0].obs_names + '-' + $
            obs_info[n_frames-1].obs_names + strmid(plotfiles[0], temp + temp2)
            
          movie_file = movie_file + '_obsid_movie.gif'
        end
        'file': begin
          movie_file = plotfiles[0] + '_file_movie.gif'
        end
      endcase
    endif
    
    cmd = 'convert -delay 20 -loop 0 ' + strjoin(plotfiles, ' ') + ' ' + movie_file
    
    Spawn, cmd, result, error_result
    
    if error_result[0] ne "" then begin
      for k=0,N_Elements(error_result)-1 DO Print, error_result[k]
    endif else Print, 'Output file located here: ' + movie_file
    
    ; Have you been asked to delete the image files?
    if Keyword_Set(delete_images) then begin
      file_delete, plotfiles
    endif
  endif
  
end