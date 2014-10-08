pro hellebore_fhd_cube_images, folder_names, obs_names_in, cube_types = cube_types, pols = pols, evenodd = evenodd, $
    rts = rts, sim = sim, casa = casa, png = png, eps = eps, slice_range = slice_range, ratio = ratio, diff_ratio = diff_ratio, $
    log = log, data_range = data_range, color_profile = color_profile, sym_color = sym_color, window_num = window_num, plot_as_map = plot_as_map
    
  if n_elements(folder_names) gt 2 then message, 'No more than 2 folder_names can be supplied'
  if n_elements(evenodd) eq 0 then evenodd = 'even'
  if n_elements(evenodd) gt 2 then message, 'No more than 2 evenodd values can be supplied'
  if n_elements(obs_names_in) gt 2 then message, 'No more than 2 obs_names can be supplied'
  
  obs_info = hellebore_filenames(folder_names, obs_names_in, rts = rts, sim = sim, casa = casa)
  
  filenames = strarr(max([n_elements(obs_info.obs_names), n_elements(evenodd)]))
  
  if n_elements(cube_types) eq 0 then cube_types = 'res'
  if n_elements(cube_types) gt 2 then message, 'No more than 2 cube_types can be supplied'
  if n_elements(pols) eq 0 then pols = 'xx'
  if n_elements(pols) gt 2 then message, 'No more than 2 pols can be supplied'
  
  n_cubes = max([n_elements(filenames), n_elements(cube_types), n_elements(pols)])
  
  if keyword_set(ratio) and keyword_set(diff_ratio) then message, 'only one of ratio & diff_ratio keywords can be set'
  
  if keyword_set(diff_ratio) and n_cubes eq 1 then begin
    print, 'diff_ratio keyword only applies when 2 cubes are specified.'
    undefine, diff_ratio
  endif
  
  if keyword_set(ratio) and n_cubes eq 1 then begin
    print, 'ratio keyword only applies when 2 cubes are specified.'
    undefine, ratio
  endif
  
  if n_cubes eq 2 and n_elements(data_range) eq 0 and n_elements(sym_color) eq 0 and not keyword_set(ratio) then sym_color=1
  if keyword_set(sym_color) and keyword_set(log) then color_profile = 'sym_log'
  
  for i=0, n_elements(folder_names)-1 do begin
  
    if n_elements(filenames) eq 1 then begin
      ;; only 1 folder name & 1 evenodd
      evenodd_mask = stregex(obs_info.cube_files.(i), evenodd, /boolean)
      if total(evenodd_mask) gt 0 then filenames = obs_info.cube_files.(i)[(where(evenodd_mask eq 1))[0]] else message, 'requested file does not exist'
    endif else begin
      ;; 2 of folder name and/or evenodd
      if n_elements(evenodd) eq 1 then begin
        ;; 2 folder names, 1 evenodd
        evenodd_mask = stregex(obs_info.cube_files.(i), evenodd, /boolean)
        if total(evenodd_mask) gt 0 then filenames[i] = obs_info.cube_files.(i)[(where(evenodd_mask eq 1))[0]] else message, 'requested file does not exist'
      endif else begin
        if n_elements(folder_names) gt 1 then begin
          ;; 2 of each folder name & evenodd
          evenodd_mask = stregex(obs_info.cube_files.(i), evenodd[i], /boolean)
          if total(evenodd_mask) gt 0 then filenames[i] = obs_info.cube_files.(i)[(where(evenodd_mask eq 1))[0]] else message, 'requested file does not exist'
        endif else begin
          ;; 1 foldername, 2 evenodd
          for j=0, n_elements(evenodd)-1 do begin
            evenodd_mask = stregex(obs_info.cube_files.(i), evenodd[j], /boolean)
            if total(evenodd_mask) gt 0 then filenames[j] = obs_info.cube_files.(i)[(where(evenodd_mask eq 1))[0]] else message, 'requested file does not exist'
          endfor
        endelse
      endelse
    endelse
    
  endfor
  
  
  if n_elements(obs_info.folder_names) eq 2 then begin
    save_path = obs_info.diff_save_path
    note = obs_info.diff_note
    if keyword_set(ratio) then note = strjoin(strsplit(note, '-', /extract), '/')
    plot_path = obs_info.diff_plot_path
  endif else begin
    save_path = obs_info.folder_names[0] + path_sep()
    if keyword_set(rts) then note = obs_info.rts_types[0] else note = obs_info.fhd_types[0]
    plot_path = obs_info.plot_paths[0]
  endelse
  
  if file_test(save_path) eq 0 then file_mkdir, save_path
  
  max_file = n_elements(filenames)-1
  max_type = n_elements(cube_types)-1
  max_pol = n_elements(pols)-1
  max_eo = n_elements(evenodd)-1
  
  
  ;; title to use:
  if n_cubes gt 1 then begin
    if keyword_set(ratio) then cube_op = '/' else cube_op = '-'
    
    if n_elements(folder_names) eq 1 then diff_title = evenodd[0] + '_' + cube_types[0] + '_' + pols[0] + $
      cube_op + evenodd[max_eo] + '_' + cube_types[max_type] + '_' + pols[max_pol] $
    else $
      diff_title = evenodd[0] + '_' + cube_types[0] + '_' + pols[0] + $
      cube_op + evenodd[max_eo] + '_' + cube_types[max_type] + '_' + pols[max_pol]
  endif else diff_title = evenodd[0] + '_' + cube_types[0] + '_' + pols[0]
  
  if keyword_set(rts) then pixel_varname = 'pixel_nums' else pixel_varname = 'hpx_inds'
  
  hpx_inds1 = getvar_savefile(filenames[0], pixel_varname)
  if n_elements(filenames) gt 1 then begin
    hpx_inds2 = getvar_savefile(filenames[1], pixel_varname)
    if total(abs(hpx_inds2-hpx_inds1)) gt 0 then message, 'healpix pixels do not match between the 2 files'
  endif
  
  nside1 = getvar_savefile(filenames[0], 'nside')
  if n_elements(filenames) gt 1 then begin
    nside2 = getvar_savefile(filenames[1], 'nside')
    if total(abs(hpx_inds2-hpx_inds1)) gt 0 then message, 'nsides do not match between the 2 files'
  endif
  
  pol_exist = stregex(filenames, '[xy][xy]', /boolean, /fold_case)
  
  if keyword_set(rts) then begin
    data_varname1 = pols[0] + '_data'
    if n_cubes gt 1 then data_varname2 = pols[max_pol] + '_data'
  endif else begin
    if pol_exist[0] then data_varname1 = cube_types[0] + '_cube' else data_varname1 = cube_types[0] + '_' + pols[0] + '_cube'
    if n_cubes gt 1 then if pol_exist[max_file] then data_varname2 = cube_types[max_type] + '_cube' else data_varname2 = cube_types[max_type] + '_' + pols[max_pol] + '_cube'
  endelse
  
  cube1 = getvar_savefile(filenames[0], data_varname1)
  n_freq1 = (size(cube1,/dimension))[1]
  if n_cubes gt 1 then begin
    cube2 = getvar_savefile(filenames[max_file], data_varname2)
    n_freq2 = (size(cube2,/dimension))[1]
    if n_freq1 ne n_freq2 then message, 'number of frequencies do not match between the 2 files'
  endif
  
  if obs_info.info_files[0] ne '' then file_struct_arr1 = fhd_file_setup(obs_info.info_files[0], pols[0])
  if n_elements(filenames) eq 2 then if obs_info.info_files[1] ne '' then $
    file_struct_arr2 = fhd_file_setup(obs_info.info_files[1], pols[max_pol])
    
  if n_elements(file_struct_arr1) ne 0 then begin
    if cube_types[0] eq 'weights' or cube_types[0] eq 'variance' then wh_match = where(file_struct_arr1.pol eq pols[0], count_match) $
    else wh_match = where(file_struct_arr1.pol eq pols[0] and file_struct_arr1.type eq cube_types[0], count_match)
    if count_match gt 0 then begin
      wh_eo = where(stregex(file_struct_arr1[wh_match[0]].uvf_label, evenodd[max_eo], /boolean) eq 1, count_eo)
      if count_eo gt 0 then nvis1 = file_struct_arr1[wh_match[0]].n_vis[wh_eo[0]]
    endif
  endif
  
  if n_elements(filenames) eq 2 and n_elements(file_struct_arr2) ne 0 then begin
    if cube_types[0] eq 'weights' or cube_types[0] eq 'variance' then wh_match = where(file_struct_arr2.pol eq pols[max_pol], count_match) $
    else wh_match = where(file_struct_arr2.pol eq pols[max_pol] and file_struct_arr2.type eq cube_types[max_type], count_match)
    if count_match gt 0 then begin
      wh_eo = where(stregex(file_struct_arr2[wh_match[0]].uvf_label, evenodd[max_eo], /boolean) eq 1, count_eo)
      if count_eo gt 0 then nvis2 = file_struct_arr2[wh_match[0]].n_vis[wh_eo[0]]
    endif
  endif
  
  if n_elements(filenames) eq 2 then begin
    if n_elements(nvis1) gt 0 and n_elements(nvis2) gt 0 then begin
      print, 'n_vis % difference between files: ' + number_formatter((nvis2-nvis1)*100/nvis1)
    endif
  endif else if n_elements(nvis1) gt 0 then print, 'nvis: ' + number_formatter(nvis1)
  
  print, 'nside, n pixels, n_freq: ' + number_formatter(nside1) + ', ' + number_formatter(n_elements(hpx_inds1)) + ', ' + number_formatter(n_freq1)
  
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
  
    if not file_test(plot_path, /directory) then file_mkdir, plot_path
    
    ;; plot_filebase specifies a base name to use for the plot files
    if n_cubes gt 1 then begin
      if n_elements(folder_names) eq 1 then begin
        if n_elements(obs_info.obs_names) gt 1 then begin
          plot_filebase = obs_info.folder_basenames[0] + '_' + obs_info.obs_names[0] + '_' + evenodd[0] + '_' + cube_types[0] + '_' + pols[0] + $
            '_minus_' + obs_info.obs_names[0] + '_' + evenodd[max_eo] + '_' + cube_types[max_type] + '_' + pols[max_pol]
        endif else begin
          if obs_info.integrated[0] eq 0 then plot_start = obs_info.folder_basenames[0] + '_' + obs_info.obs_names[0] else plot_start = obs_info.fhd_types[0]
          
          plot_filebase = plot_start + '_' + evenodd[0] + '_' + cube_types[0] + '_' + pols[0] + $
            '_minus_' + evenodd[max_eo] + '_' + cube_types[max_type] + '_' + pols[max_pol]
        endelse
      endif else plot_filebase = obs_info.name_same_parts + '__' + strjoin([obs_info.name_diff_parts[0], evenodd[0], cube_types[0], pols[0]], '_')  + $
        '_minus_' + strjoin([obs_info.name_diff_parts[1], evenodd[max_eo], cube_types[max_type], pols[max_pol]], '_')
    endif else begin
      if obs_info.integrated[0] eq 0 then plot_start = obs_info.folder_basenames[0] + '_' + obs_info.obs_names[0] else plot_start = obs_info.fhd_types[0]
      
      plot_filebase = plot_start + '_' + evenodd[0] + '_' + cube_types[0] + '_' + pols[0]
    endelse
    
    if keyword_set(diff_ratio) then plotfile = plot_path + plot_filebase + '_imagenormdiff' $
    else if keyword_set(ratio) then plotfile = plot_path + plot_filebase + '_imageratio' else plotfile = plot_path + plot_filebase + '_image'
  endif
  
  if n_cubes gt 1 then begin
    if max(abs(cube1-cube2)) eq 0 then message, 'cubes are identical.'
    if keyword_set(diff_ratio) then begin
      print, max(cube1), max(cube2), max(cube1)/max(cube2)
      temp = (cube1/max(cube1) - cube2/max(cube2)) * mean([max(cube1), max(cube2)])
      note = note + ', peak ratio = ' + number_formatter(max(cube1)/max(cube2), format = '(f5.2)')
    endif else if keyword_set(ratio) then temp = cube1/cube2 else temp = cube1-cube2
    
    if n_elements(slice_range) eq 1 then temp = temp[*,slice_range] else temp = total(temp[*, slice_range[0]:slice_range[1]],2)
  endif else if n_elements(slice_range) eq 1 then temp = cube1[*,slice_range] else temp = total(cube1[*, slice_range[0]:slice_range[1]],2)
  
  if keyword_set(sym_color) and not keyword_set(log) then begin
    if n_elements(data_range) eq 0 then data_range = [-1,1]*max(abs(temp)) $
    else data_range = [-1,1]*max(abs(data_range))
  endif
  if keyword_set(diff_ratio) then title = diff_title + ', peak norm., ' + title_range else title = diff_title + ', ' + title_range
  
  healpix_quickimage, temp, hpx_inds1, nside1, title = title, savefile = plotfile, note=note, slice_ind = slice_ind, $
    log = log, color_profile = color_profile, data_range = data_range, window_num = window_num, plot_as_map = plot_as_map
    
end
