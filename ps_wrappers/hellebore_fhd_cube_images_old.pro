pro hellebore_fhd_cube_images_old, folder_name, obs_range, png = png, eps = eps

  if n_elements(folder_name) eq 0 then folder_name = base_path('data') + 'fhd_ps_data/128T_cubes/aug23_3hr_first'
  
  ;; check for folder existence, otherwise look for common folder names to figure out full path. If none found, try base_path('data') + 'fhd_ps_data/128T_cubes/'
  folder_test = file_test(folder_name, /directory)
  if folder_test eq 0 then begin
    pos_fhd_data = strpos(folder_name, 'fhd_ps_data')
    if pos_fhd_data gt -1 then begin
      test_name = base_path('data') + strmid(folder_name, pos_fhd_data)
      folder_test = file_test(test_name, /directory)
      if folder_test eq 1 then folder_name = test_name
    endif
  endif
  if folder_test eq 0 then begin
    pos_fhd_128 = strpos(folder_name, '128T_cubes')
    if pos_fhd_128 gt -1 then begin
      test_name = base_path('data') + 'fhd_ps_data/' + strmid(folder_name, pos_fhd_128)
      folder_test = file_test(test_name, /directory)
      if folder_test eq 1 then folder_name = test_name
    endif
  endif
  if folder_test eq 0 then begin
    test_name = base_path('data') + 'fhd_ps_data/128T_cubes/' + folder_name
    folder_test = file_test(test_name, /directory)
    if folder_test eq 1 then folder_name = test_name
  endif
  
  if folder_test eq 0 then message, 'folder not found'
  
  fhd_type = file_basename(folder_name)
  
  if n_elements(obs_range) gt 0 then begin
    if size(obs_range,/type) eq 7 then begin
      if n_elements(obs_range) gt 1 then $
        message, 'obs_range can be specified as a single string to use as the name or as a 2 element obsid range'
      obs_name = obs_range
      obs_range = long(strsplit(obs_name, '-', /extract))
      if n_elements(obs_range) eq 1 then obs_name_single = obs_name
    endif else begin
      if n_elements(obs_range) gt 2 then message, 'obs_range can be specified as a single string to use as the name or as a 2 element obsid range'
      if n_elements(obs_range) eq 2 then obs_name = number_formatter(obs_range[0]) + '-' + number_formatter(obs_range[1]) else begin
        obs_name = number_formatter(obs_range[0]) + '-' + number_formatter(obs_range[0])
        obs_name_single = number_formatter(obs_range[0])
      endelse
    endelse
  endif else begin
    obs_name = ''
    obs_name_single = ''
  endelse
  
  ;; first look for integrated cube files with names like Combined_obs_...
  cube_files = file_search(folder_name + '/Combined_obs_' + obs_name + '*_cube.sav', count = n_cubefiles)
  if n_cubefiles gt 0 then begin
    if obs_name eq '' then begin
      obs_name_arr = stregex(cube_files, '[0-9]+-[0-9]+', /extract)
      wh_first = where(obs_name_arr eq obs_name_arr[0], count_first)
      if count_first lt n_elements(cube_files) then $
        print, 'More than one obs_range found, using first range (' + obs_name_arr[0] + ', ' + number_formatter(count_first) + ' files)'
      if count_first gt 2 then message, 'More than two cubes found with first obs_range'
      datafile = cube_files[wh_first]
      obs_name = obs_name_arr[0]
      obs_range = long(strsplit(obs_name, '-', /extract))
    endif else begin
      if n_elements(cube_files) gt 2 then message, 'More than two cubes found with given obs_range'
      
      obs_name_arr = stregex(cube_files, '[0-9]+-[0-9]+', /extract)
      if obs_name_arr[1] ne obs_name_arr[0] then message, 'Cube files do not have the same obs ranges.'
      obs_name = obs_name_arr[0]
      obs_range = long(strsplit(obs_name, '-', /extract))
    endelse
    
  endif else if n_elements(obs_range) lt 2 then begin
    ;; then look for single obs cube files
    cube_files = file_search(folder_name + '/' + obs_name + '*_cube.sav', count = n_cubefiles)
    if n_cubefiles gt 0 then begin
      cube_basename = file_basename(cube_files)
      if obs_name eq '' then begin
        obs_name_arr = stregex(cube_basename, '[0-9]+', /extract)
        wh_first = where(obs_name_arr eq obs_name_arr[0], count_first)
        if count_first lt n_elements(cube_files) then $
          print, 'More than one obs_range found, using first range (' + obs_name_arr[0] + ', ' + number_formatter(count_first) + ' files)'
        if count_first gt 2 then message, 'More than two cubes found with first obs_range'
        datafile = cube_files[wh_first]
        obs_name = obs_name_arr[0]
        obs_range = long(obs_name)
      endif else begin
        if n_elements(cube_files) gt 2 then message, 'More than two cubes found with given obs_range'
        
        datafile = cube_files
        obs_name_arr = stregex(cube_basename, '[0-9]+', /extract)
        if obs_name_arr[1] ne obs_name_arr[0] then message, 'Cube files do not have the same obs ranges.'
        obs_name = obs_name_arr[0]
        obs_range = long(obs_name)
      endelse
      
    endif
  endif
  
  if n_elements(datafile) eq 0 then message, 'No cube or info files found in folder ' + folder_name
  
  note = fhd_type
  
  if keyword_set(png) or keyword_set(eps) then pub = 1 else pub = 0
  
  
  std_savepath = base_path('data') + 'fhd_ps_data/'
  pos = strpos(file_dirname(datafile[0], /mark_directory), std_savepath)
  if pos ne -1 then save_path_ext = strmid(file_dirname(datafile[0], /mark_directory), pos + strlen(std_savepath)) $
  else save_path_ext = ''
  
  if pub then begin
    plot_path = base_path('plots') + 'power_spectrum/fhd_data/' + save_path_ext
    if not file_test(plot_path, /directory) then file_mkdir, plot_path
    
    plot_filebase = fhd_type + '_' + obs_name
  endif
  
  even_mask = stregex(datafile, 'even', /boolean)
  odd_mask = stregex(datafile, 'odd', /boolean)
  
  if total(even_mask ) gt 0 then even_filename = datafile[(where(even_mask eq 1))[0]]
  if total(odd_mask ) gt 0 then odd_filename = datafile[(where(odd_mask eq 1))[0]]
  
  if n_elements(even_filename) gt 0 then even = 1 else even = 0
  if n_elements(odd_filename) gt 0 then odd = 1 else even = 0
  if even and odd then both = 1
  if even eq 0 and odd eq 0 then begin
    print, 'no even or odd filename'
    return
  endif
  
  if even then begin
    hpx_inds = getvar_savefile(even_filename, 'hpx_inds')
    nside = getvar_savefile(even_filename, 'nside')
  endif else begin
    hpx_inds = getvar_savefile(odd_filename, 'hpx_inds')
    nside = getvar_savefile(odd_filename, 'nside')
  endelse
  
  if even then begin
    even_res = getvar_savefile(even_filename, 'res_xx_cube')
    even_model = getvar_savefile(even_filename, 'model_xx_cube')
    even_weights = getvar_savefile(even_filename, 'weights_xx_cube')
    even_dirty = getvar_savefile(even_filename, 'dirty_xx_cube')
  endif
  if odd then begin
    odd_res = getvar_savefile(odd_filename, 'res_xx_cube')
    odd_model = getvar_savefile(odd_filename, 'model_xx_cube')
    odd_weights = getvar_savefile(odd_filename, 'weights_xx_cube')
    odd_dirty = getvar_savefile(odd_filename, 'dirty_xx_cube')
  endif
  if both then begin
    diff_res = even_res - odd_res
    diff_weights = even_weights - odd_weights
    diff_model = even_model - odd_model
    diff_dirty = even_dirty - odd_dirty
  endif
  
  if both then begin
    if pub then begin
      weights_plotfile = plot_path + plot_filebase + '_weights_diff_image.png'
      model_plotfile = plot_path + plot_filebase + '_model_diff_image.png'
      res_plotfile = plot_path + plot_filebase + '_res_diff_image.png'
      dirty_plotfile = plot_path + plot_filebase + '_diff_dirty_image.png'
    endif
    healpix_quickimage, total(diff_res,2), hpx_inds, nside, title = 'res even-odd, freq. averaged', savefile = res_plotfile
    healpix_quickimage, total(diff_model,2), hpx_inds, nside, title = 'model even-odd, freq. averaged', savefile = model_plotfile
    healpix_quickimage, total(diff_weights,2), hpx_inds, nside, title = 'weights even-odd, freq. averaged', savefile = weights_plotfile
    healpix_quickimage, total(diff_dirty,2), hpx_inds, nside, title = 'dirty even-odd, freq. averaged', savefile = dirty_plotfile
  endif
  
  if even then begin
    if pub then begin
      even_res_plotfile = plot_path + plot_filebase + '_even_res_image.png'
      even_model_plotfile = plot_path + plot_filebase + '_even_model_image.png'
      even_dirty_plotfile = plot_path + plot_filebase + '_even_dirty_image.png'
    endif
    healpix_quickimage, total(even_res,2), hpx_inds, nside, title = 'res even, freq. averaged', savefile = even_res_plotfile
    healpix_quickimage, total(even_model,2), hpx_inds, nside, title = 'model even, freq. averaged', savefile = even_model_plotfile
    healpix_quickimage, total(even_dirty,2), hpx_inds, nside, title = 'dirty even, freq. averaged', savefile = even_dirty_plotfile
  endif else begin
    if pub then begin
      odd_res_plotfile = plot_path + plot_filebase + '_odd_res_image.png'
      odd_model_plotfile = plot_path + plot_filebase + '_odd_model_image.png'
      odd_dirty_plotfile = plot_path + plot_filebase + '_odd_dirty_image.png'
    endif
    healpix_quickimage, total(odd_res,2), hpx_inds, nside, title = 'res odd, freq. averaged', savefile = odd_res_plotfile
    healpix_quickimage, total(odd_model,2), hpx_inds, nside, title = 'model odd, freq. averaged', savefile = odd_model_plotfile
    healpix_quickimage, total(odd_dirty,2), hpx_inds, nside, title = 'dirty odd, freq. averaged', savefile = odd_dirty_plotfile
  endelse
  
  
end