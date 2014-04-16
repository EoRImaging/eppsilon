pro hellebore_fhd_cube_images, folder_names, cube_types = cube_types, pols = pols, evenodd = evenodd, $
    png = png, eps = eps, slice_range = slice_range
    
    
  if n_elements(folder_names) eq 0 then folder_names = base_path('data') + 'fhd_ps_data/128T_cubes/aug23_3hr_first/'
  if n_elements(evenodd) eq 0 then evenodd = 'even'
  
  filenames = strarr(max([n_elements(folder_names), n_elements(evenodd)]))
  
  if n_elements(cube_types) eq 0 then cube_types = 'res'
  if n_elements(pols) eq 0 then pols = 'xx'
  
  n_cubes = max([n_elements(filenames), n_elements(cube_types), n_elements(pols)])
  
  info_files = strarr(n_elements(folder_names))
  obs_names = strarr(n_elements(folder_names))
  fhd_types = strarr(n_elements(folder_names))
  
  for i=0, n_elements(folder_names)-1 do begin
    ;; check for folder existence, otherwise look for common folder names to figure out full path. If none found, try base_path('data') + 'fhd_ps_data/128T_cubes/'
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
    
    if folder_test eq 0 then message, 'folder not found'
    
    fhd_types[i] = file_basename(folder_names[i])
    
    ;; first look for integrated cube files with names like Combined_obs_...
    cube_files = file_search(folder_names[i] + '/Combined_obs_*_cube.sav', count = n_cubefiles)
    if n_cubefiles gt 0 then begin
      obs_name_arr = stregex(cube_files, '[0-9]+-[0-9]+', /extract)
      wh_first = where(obs_name_arr eq obs_name_arr[0], count_first)
      if count_first lt n_elements(cube_files) then $
        print, 'More than one obs_range found, using first range (' + obs_name_arr[0] + ', ' + number_formatter(count_first) + ' files)'
      if count_first gt 2 then message, 'More than two cubes found with first obs_range'
      datafile = cube_files[wh_first]
      obs_names[i] = obs_name_arr[0]
      
    endif else if n_elements(obs_range) lt 2 then begin
      ;; then look for single obs cube files
      cube_files = file_search(folder_names[i] + '/*_cube.sav', count = n_cubefiles)
      if n_cubefiles gt 0 then begin
        cube_basename = file_basename(cube_files)
        obs_name_arr = stregex(cube_basename, '[0-9]+', /extract)
        wh_first = where(obs_name_arr eq obs_name_arr[0], count_first)
        if count_first lt n_elements(cube_files) then $
          print, 'More than one obs_range found, using first range (' + obs_name_arr[0] + ', ' + number_formatter(count_first) + ' files)'
        if count_first gt 2 then message, 'More than two cubes found with first obs_range'
        datafile = cube_files[wh_first]
        obs_names[i] = obs_name_arr[0]
        
      endif
    endif
    
    if n_elements(filenames) eq 1 then begin
      ;; only 1 folder name & 1 evenodd
      evenodd_mask = stregex(datafile, evenodd, /boolean)
      if total(evenodd_mask) gt 0 then filenames = datafile[(where(evenodd_mask eq 1))[0]] else message, 'requested file does not exist'
    endif else begin
      ;; 2 of folder name and/or evenodd
      if n_elements(evenodd) eq 1 then begin
        ;; 2 folder names, 1 evenodd
        evenodd_mask = stregex(datafile, evenodd, /boolean)
        if total(evenodd_mask) gt 0 then filenames[i] = datafile[(where(evenodd_mask eq 1))[0]] else message, 'requested file does not exist'
      endif else begin
        if n_elements(folder_names) gt 1 then begin
          ;; 2 of each folder name & evenodd
          evenodd_mask = stregex(datafile, evenodd[i], /boolean)
          if total(evenodd_mask) gt 0 then filenames[i] = datafile[(where(evenodd_mask eq 1))[0]] else message, 'requested file does not exist'
        endif else begin
          ;; 1 foldername, 2 evenodd
          for j=0, n_elements(evenodd)-1 do begin
            evenodd_mask = stregex(datafile, evenodd[j], /boolean)
            if total(evenodd_mask) gt 0 then filenames[j] = datafile[(where(evenodd_mask eq 1))[0]] else message, 'requested file does not exist'
          endfor
        endelse
      endelse
    endelse
    undefine, datafile
  endfor
  
  ;; save_path specifies a location to save the power spectrum files.
  ;; This is also where the code looks for intermediate save files to avoid re-running code.
  
  std_savepath = base_path('data') + 'fhd_ps_data/'
  
  if n_elements(folder_names) eq 2 then begin
    folderparts_1 = strsplit(folder_names[0], path_sep(), /extract)
    folderparts_2 = strsplit(folder_names[1], path_sep(), /extract)
    match_test = strcmp(folderparts_1, folderparts_2)
    wh_diff = where(match_test eq 0, count_diff, complement = wh_same, ncomplement = count_same)
    
    if count_diff eq 0 then begin
      ;; folders are the same
      folder_names = folder_names[0]
      note = file_basename(info_file)
      save_path = file_dirname(info_file, /mark_directory)
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
        note = strjoin(folderparts_1[wh_diff], path_sep()) + ' - ' + strjoin(folderparts_2[wh_diff], path_sep()) + ' ' + file_basename(folder_names[0])
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
          note = str_same + ' ' + str1_diff + ' - ' + str2_diff
        endif else begin
          diff_dir = str1_diff + '_minus_' + str2_diff
          note = str1_diff + ' - ' + str2_diff
        endelse
      endelse
      
      save_path = joint_path + path_sep() + diff_dir + path_sep()
    endelse
    if file_test(save_path) eq 0 then file_mkdir, save_path
  endif else begin
    save_path = folder_names[0] + path_sep()
    note = fhd_types[0]
  endelse
  
  pos = strpos(save_path, std_savepath)
  if pos ne -1 then save_path_ext = strmid(save_path, pos + strlen(std_savepath)) $
  else save_path_ext = ''
  
  max_file = n_elements(filenames)-1
  max_type = n_elements(cube_types)-1
  max_pol = n_elements(pols)-1
  max_eo = n_elements(evenodd)-1
  
  
  ;; title to use:
  if n_cubes gt 1 then begin
    if n_elements(folder_names) eq 1 then diff_title = evenodd[0] + '_' + cube_types[0] + '_' + pols[0] + $
      ' - ' + evenodd[max_eo] + '_' + cube_types[max_type] + '_' + pols[max_pol] $
    else $
      diff_title = evenodd[0] + '_' + cube_types[0] + '_' + pols[0] + $
      ' - ' + evenodd[max_eo] + '_' + cube_types[max_type] + '_' + pols[max_pol]
  endif else diff_title = evenodd[0] + '_' + cube_types[0] + '_' + pols[0]
  
  hpx_inds1 = getvar_savefile(filenames[0], 'hpx_inds')
  if n_elements(filenames) gt 1 then begin
    hpx_inds2 = getvar_savefile(filenames[1], 'hpx_inds')
    if total(abs(hpx_inds2-hpx_inds1)) gt 0 then message, 'healpix pixels do not match between the 2 files'
  endif
  
  nside1 = getvar_savefile(filenames[0], 'nside')
  if n_elements(filenames) gt 1 then begin
    nside2 = getvar_savefile(filenames[1], 'nside')
    if total(abs(hpx_inds2-hpx_inds1)) gt 0 then message, 'nsides do not match between the 2 files'
  endif
    
  cube1 = getvar_savefile(filenames[0], cube_types[0] + '_' + pols[0] + '_cube')
  n_freq1 = (size(cube1,/dimension))[1]
  if n_cubes gt 1 then begin
    cube2 = getvar_savefile(filenames[max_file], cube_types[max_type] + '_' + pols[max_pol] + '_cube')
    n_freq2 = (size(cube2,/dimension))[1]
    if n_freq1 ne n_freq2 then message, 'number of frequencies do not match between the 2 files'
  endif
  
  print, 'nside, n pixels: ' + number_formatter(nside1) + ', ' + number_formatter(n_elements(hpx_inds1))
  
  case n_elements(slice_range) of
    0: begin
      slice_range = [0, n_freq1-1]
      title_range = 'freq. added'
    end
    1: begin
       title_range = 'slice ' + number_formatter(slice_range)
    end
    2: begin
      if min(slice_range) lt 0 then message, 'slice_range cannot be less than zero'
      if max(slice_range) ge n_freq1 then message, 'slice_range cannot be more than ' + number_formatter(n_freq1-1)
      if slice_range[1] lt slice_range[0] then message, 'slice_range[1] cannot be less than slice_range[0]'
      
      title_range = 'slices [' + number_formatter(slice_range[0]) + ':' + number_formatter(slice_range[1]) + ']'
    end
    else: begin
      message, 'slice_range must be a 1 or 2 element vector'
    end
  endcase
  
  
  if keyword_set(png) or keyword_set(eps) then pub = 1 else pub = 0
  if pub then begin
    ;; plot_path specifies a location to save plot files.
    plot_path = base_path('plots') + 'power_spectrum/fhd_data/' + save_path_ext
    
    if not file_test(plot_path, /directory) then file_mkdir, plot_path
    
    ;; plot_filebase specifies a base name to use for the plot files
    if n_cubes gt 1 then begin
      if n_elements(folder_names) eq 1 then plot_filebase = fhd_types[0] + '_' + evenodd[0] + '_' + cube_types[0] + '_' + pols[0] + $
        '_minus_' + evenodd[max_eo] + '_' + cube_types[max_type] + '_' + pols[max_pol] $
      else $
        plot_filebase = strjoin(fnameparts_1[wh_name_same], '_') + '__' + strjoin([fnameparts_1[wh_name_diff], evenodd[0], cube_types[0], pols[0]], '_')  + $
        '_minus_' + strjoin([fnameparts_2[wh_name_diff], evenodd[max_eo], cube_types[max_type], pols[max_pol]], '_')
    endif else plot_filebase = fhd_types[0] + '_' + evenodd[0] + '_' + cube_types[0] + '_' + pols[0]
    
    plotfile = plot_path + plot_filebase + '_image'
  endif
  
  if n_cubes gt 1 then begin
    temp = cube1-cube2
    if max(abs(temp)) eq 0 then message, 'cubes are identical.'
    if n_elements(slice_range) eq 1 then temp = temp[*,slice_range] else temp = total(temp[*, slice_range[0]:slice_range[1]],2)
  endif else if n_elements(slice_range) eq 1 then temp = cube1[*,slice_range] else temp = total(cube1[*, slice_range[0]:slice_range[1]],2)
  
  title = diff_title + ', ' + title_range
  
  healpix_quickimage, temp, hpx_inds1, nside1, title = title, savefile = plotfile, note=note, slice_ind = slice_ind
  
end
