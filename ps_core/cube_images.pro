pro cube_images, folder_names, obs_info, nvis_norm = nvis_norm, pols = pols, cube_types = cube_types, evenodd = evenodd, $
    rts = rts, png = png, eps = eps, pdf = pdf, slice_range = slice_range, sr2 = sr2, $
    ratio = ratio, diff_ratio = diff_ratio, diff_frac = diff_frac, $
    log = log, data_range = data_range, data_min_abs = data_min_abs, color_profile = color_profile, sym_color = sym_color, $
    window_num = window_num, plot_as_map = plot_as_map, plot_path = plot_path, plotfile = plotfile_out
    
  if n_elements(obs_info.info_files) gt 2 then message, 'Only 1 or 2 info_files can be used'
  
  if n_elements(cube_types) eq 0 then cube_types = 'dirty'
  if n_elements(cube_types) gt 2 then message, 'No more than 2 cube_types can be supplied'
  if n_elements(pols) eq 0 then pols = 'xx'
  if n_elements(pols) gt 2 then message, 'No more than 2 pols can be supplied'
  
  if n_elements(sr2) gt 0 then begin
    if n_elements(slice_range) eq 0 then message, 'sr2 can only be set if slice_range is also set'
    if max(sr2) eq max(slice_range) and min(sr2) eq min(slice_range) then n_slice_range = 1 else n_slice_range = 2
  endif else n_slice_range = 1
  
  n_cubes = max([n_elements(obs_info.info_files), n_elements(evenodd), n_elements(cube_types), $
    n_elements(pols), n_slice_range])
    
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
  if n_elements(color_profile) gt 0 then if color_profile eq 'sym_log' then begin
    log=1
    sym_color=1
  endif
  
  max_file = max([n_elements(obs_info.info_files)-1, n_elements(evenodd)-1])
  max_type = n_elements(cube_types)-1
  max_pol = n_elements(pols)-1
  max_eo = n_elements(evenodd)-1
s
  if obs_info.info_files[0] ne '' then begin
    if keyword_set(rts) then file_struct_arr1 = rts_file_setup(obs_info.info_files[0]) $
    else $
      if keyword_set(casa) then file_struct_arr1 = casa_file_setup(obs_info.info_files[0]) $
    else file_struct_arr1 = fhd_file_setup(obs_info.info_files[0])
  endif else message, 'No info file exists for folder: ' + folder_names[0] + $
    ' and obs_name: ' + obs_info.obs_names[0]
    
  if n_elements(obs_info.info_files) eq 2 then if obs_info.info_files[1] ne '' then begin
    if keyword_set(rts) then file_struct_arr2 = rts_file_setup(obs_info.info_files[1]) $
    else $
      if keyword_set(casa) then file_struct_arr2 = casa_file_setup(obs_info.info_files[1]) $
    else file_struct_arr2 = fhd_file_setup(obs_info.info_files[1])
  endif else message, 'No info file exists for folder: ' + folder_names[n_elements(folder_names)-1] + $
    ' and obs_name: ' + obs_info.obs_names[n_elements(obs_info.obs_names)-1]
  s
  type_pol_str1 = file_struct_arr1.type_pol_str
  if n_elements(file_struct_arr2) gt 0 then type_pol_str2 = file_struct_arr2.type_pol_str
  
  type_pol_str = [cube_types[0] + '_' + pols[0], cube_types[max_type] + '_' + pols[max_pol]]
  
  if cube_types[0] eq 'weights' or cube_types[0] eq 'variances' then begin
    cube_ind = where(file_struct_arr1.pol eq pols[0], count_pol)
    if count_pol eq 0 then message, 'Specified pol is not available'
    if count_pol gt 0 then cube_ind = cube_ind[0]
    
    evenodd_mask = stregex(file_struct_arr1[cube_ind].weightfile, evenodd[0], /boolean)
    evenodd_ind = where(evenodd_mask eq 1, count_evenodd)
    if count_evenodd eq 0 then message, 'specified evenodd value is not available'
    if count_evenodd gt 1 then message, 'More than one matching evenodd value'
    
    filenames = file_struct_arr1[cube_ind].weightfile[evenodd_ind]
    if cube_types[0] eq 'weights' then cube_varnames = file_struct_arr1[cube_ind].weightvar $
    else cube_varnames = file_struct_arr1[cube_ind].variancevar
  endif else begin
    cube_ind = where(type_pol_str1 eq type_pol_str[0], count_typepol)
    if count_typepol eq 0 then message, 'Specified cube type and pol are not available.'
    if count_typepol gt 1 then message, 'More than one matching type and pol'
    
    evenodd_mask = stregex(file_struct_arr1[cube_ind].datafile, evenodd[0], /boolean)
    evenodd_ind = where(evenodd_mask eq 1, count_evenodd)
    if count_evenodd eq 0 then message, 'specified evenodd value is not available'
    if count_evenodd gt 1 then message, 'More than one matching evenodd value'
    
    filenames = file_struct_arr1[cube_ind].datafile[evenodd_ind]
    cube_varnames = file_struct_arr1[cube_ind].datavar
    
  endelse
  
  if n_cubes eq 2 then begin
  
    if cube_types[max_type] eq 'weights' or cube_types[max_type] eq 'variances' then begin
      if n_elements(file_struct_arr2) gt 0 then cube_ind2 = where(file_struct_arr2.pol eq pols[0], count_pol) $
      else cube_ind2 = where(file_struct_arr1.pol eq pols[0], count_pol)
      if count_pol eq 0 then message, 'Specified pol is not available'
      if count_pol gt 0 then cube_ind2 = cube_ind2[0]
      
      if n_elements(file_struct_arr2) gt 0 then $
        evenodd_mask = stregex(file_struct_arr2[cube_ind2].weightfile, evenodd[max_eo], /boolean) $
      else evenodd_mask = stregex(file_struct_arr1[cube_ind2].weightfile, evenodd[max_eo], /boolean)
      evenodd_ind2 = where(evenodd_mask eq 1, count_evenodd)
      if count_evenodd eq 0 then message, 'specified evenodd value is not available'
      if count_evenodd gt 1 then message, 'More than one matching evenodd value'
      
      if n_elements(file_struct_arr2) gt 0 then begin
        filenames = [filenames, file_struct_arr2[cube_ind2].weightfile[evenodd_ind2]]
        if cube_types[0] eq 'weights' then cube_varnames = [cube_varnames, file_struct_arr2[cube_ind2].weightvar] $
        else cube_varnames = [cube_varnames, file_struct_arr2[cube_ind2].variancevar]
      endif else begin
        filenames = [filenames, file_struct_arr1[cube_ind2].weightfile[evenodd_ind2]]
        if cube_types[0] eq 'weights' then cube_varnames = [cube_varnames, file_struct_arr1[cube_ind2].weightvar] $
        else cube_varnames = [cube_varnames, file_struct_arr1[cube_ind2].variancevar]
      endelse
    endif else begin
      if n_elements(file_struct_arr2) gt 0 then cube_ind2 = where(type_pol_str1 eq type_pol_str[1], count_typepol) $
      else cube_ind2 = where(type_pol_str1 eq type_pol_str[1], count_typepol)
      if count_typepol eq 0 then message, 'Specified cube type and pol are not available.'
      if count_typepol gt 1 then message, 'More than one matching type and pol'
      
      if n_elements(file_struct_arr2) gt 0 then $
        evenodd_mask = stregex(file_struct_arr2[cube_ind2].datafile, evenodd[max_eo], /boolean) $
      else evenodd_mask = stregex(file_struct_arr1[cube_ind2].datafile, evenodd[max_eo], /boolean)
      evenodd_ind2 = where(evenodd_mask eq 1, count_evenodd)
      if count_evenodd eq 0 then message, 'specified evenodd value is not available'
      if count_evenodd gt 1 then message, 'More than one matching evenodd value'
      
      if n_elements(file_struct_arr2) gt 0 then begin
        filenames = [filenames, file_struct_arr2[cube_ind2].datafile[evenodd_ind2]]
        cube_varnames = [cube_varnames, file_struct_arr2[cube_ind2].datavar]
      endif else begin
        filenames = [filenames, file_struct_arr1[cube_ind2].datafile[evenodd_ind2]]
        cube_varnames = [cube_varnames, file_struct_arr1[cube_ind2].datavar]
      endelse
      
    endelse
  endif
  
  wh_novarname = where(cube_varnames eq '', count_novarname)
  if count_novarname gt 0 then begin
    input_inds = intarr(n_cubes, 2)
    input_cube_varnames = strarr(n_cubes, 2)
    
    for i=0, n_cubes-1 do begin
    
      if cube_varnames[i] eq '' then begin
        if cube_types[i] eq 'res' then input_types = ['dirty', 'model'] else $
          if cube_types[i] eq 'model' then input_types = ['dirty', 'res'] else $
          message, 'No varname for this cube type and it cannot be constructed from other cubes'
          
        input_cubes_typepol = input_types + '_' + pols[0]
        for j=0, 1 do begin
          if i eq 0 or n_elements(type_pol_str2) eq 0 then typepol_use = type_pol_str1 else typepol_use = type_pol_str2
          if i eq 0 then pol_use = pols[0] else pol_use = pols[max_pol]
          wh_input = where(typepol_use eq input_types[j]+'_'+pol_use, count_input)
          if count_input eq 0 then message, 'No varname for this cube type and it cannot be constructed from other cubes'
          if count_input eq 1 then begin
            input_inds[i,j] = wh_input
            if i eq 0 or n_elements(file_struct_arr2) eq 0 then $
              input_cube_varnames[i,j] = file_struct_arr1[input_inds[i,j]].datavar $
            else input_cube_varnames[i,j] = file_struct_arr2[input_inds[i,j]].datavar
            
          endif else message, 'No varname for this cube type and it cannot be constructed from other cubes'
        endfor
        
      endif
      
    endfor
  endif
  
  
  if n_elements(obs_info.folder_names) eq 2 then begin
    note = obs_info.diff_note
    if keyword_set(ratio) then note = strjoin(strsplit(note, '-', /extract), '/')
    
    if n_elements(plot_path) eq 0 then if tag_exist(obs_info, 'diff_plot_path') then plot_path = obs_info.diff_plot_path else $
      plot_path = obs_info.diff_save_path + path_sep() + 'plots' + path_sep()
      
  endif else begin
    if keyword_set(rts) then note = obs_info.rts_types[0] else note = obs_info.fhd_types[0]
    
    if n_elements(plot_path) eq 0 then plot_path = obs_info.plot_paths[0]
  endelse
  
  if keyword_set(rts) then pixel_varnames = strarr(n_elements(filenames)) + 'pixel_nums' $
  else pixel_varnames = strarr(n_elements(filenames)) + 'hpx_inds'
  
  pol_exist = stregex(filenames, '[xy][xy]', /boolean, /fold_case)
  
  if keyword_set(png) or keyword_set(eps) or keyword_set(pdf) then pub = 1 else pub = 0
  if pub then begin
  
    if not file_test(plot_path, /directory) then file_mkdir, plot_path
    
    ;; plot_filebase specifies a base name to use for the plot files
    if n_cubes gt 1 then begin
      if n_elements(folder_names) eq 1 then begin
        if n_elements(obs_info.obs_names) gt 1 then begin
          plot_filebase = obs_info.folder_basenames[0] + '_' + obs_info.obs_names[0] + '_' + evenodd[0] + '_' + cube_types[0] + '_' + pols[0] + $
            '_minus_' + obs_info.obs_names[0] + '_' + evenodd[max_eo] + '_' + cube_types[max_type] + '_' + pols[max_pol]
        endif else begin
          plot_start = obs_info.folder_basenames[0] + '_' + obs_info.obs_names[0]
          
          plot_filebase = plot_start + '_' + evenodd[0] + '_' + cube_types[0] + '_' + pols[0] + $
            '_minus_' + evenodd[max_eo] + '_' + cube_types[max_type] + '_' + pols[max_pol]
        endelse
      endif else plot_filebase = obs_info.fhdtype_same_parts + '__' + strjoin([obs_info.fhdtype_diff_parts[0], evenodd[0], cube_types[0], pols[0]], '_')  + $
        '_minus_' + strjoin([obs_info.fhdtype_diff_parts[1], evenodd[max_eo], cube_types[max_type], pols[max_pol]], '_')
    endif else begin
      plot_start = obs_info.folder_basenames[0] + '_' + obs_info.obs_names[0]
      
      plot_filebase = plot_start + '_' + evenodd[0] + '_' + cube_types[0] + '_' + pols[0]
    endelse
    
    if keyword_set(diff_ratio) then plotfile = plot_path + plot_filebase + '_imagenormdiff' $
    else if keyword_set(ratio) then plotfile = plot_path + plot_filebase + '_imageratio' else plotfile = plot_path + plot_filebase + '_image'
    
    if n_elements(slice_range) gt 0 then begin
      plotfile = plotfile + '_ch' + number_formatter(min(slice_range))
      if n_elements(slice_range) eq 2 then plotfile = plotfile + '-' + number_formatter(max(slice_range))
    endif
    
    plotfile_out = plotfile
  endif
  
  hpx_inds1 = getvar_savefile(filenames[0], pixel_varnames[0])
  if n_elements(filenames) gt 1 then begin
    hpx_inds2 = getvar_savefile(filenames[1], pixel_varnames[1])
    if total(abs(hpx_inds2-hpx_inds1)) gt 0 then message, 'healpix pixels do not match between the 2 files'
  endif
  
  nside1 = getvar_savefile(filenames[0], 'nside')
  if n_elements(filenames) gt 1 then begin
    nside2 = getvar_savefile(filenames[1], 'nside')
    if nside1 ne nside2 gt 0 then message, 'nsides do not match between the 2 files'
  endif
  
  if n_elements(filenames) eq 2 then begin
    if n_elements(nvis1) gt 0 and n_elements(nvis2) gt 0 then begin
      print, 'n_vis % difference between files: ' + number_formatter((nvis2-nvis1)*100/nvis1)
    endif
  endif else if n_elements(nvis1) gt 0 then print, 'nvis: ' + number_formatter(nvis1)
  
  if cube_varnames[0] eq '' then begin
  
    input_cube1 = getvar_savefile(filenames[0], input_cube_varnames[0, 0])
    input_cube2 = getvar_savefile(filenames[0], input_cube_varnames[0, 1])
    
    cube1 = temporary(input_cube1) - temporary(input_cube2)
  endif else cube1 = getvar_savefile(filenames[0], cube_varnames[0])
  n_freq1 = (size(cube1,/dimension))[1]
  if keyword_set(nvis_norm) then begin
    if obs_info.integrated[0] eq 1 then obs_varname = 'obs_arr' else obs_varname = 'obs'
    obs_arr1 = getvar_savefile(filenames[0], obs_varname)
    nvis_freq = obs_arr1.nf_vis
    nvis_dims = size(nvis_freq, /dimension)
    if n_elements(nvis_dims) eq 2 then nvis_freq = total(nvis_freq, 2)
    
    n_avg = getvar_savefile(filenames[0], 'n_avg')
    n_freqbins = nvis_dims[0] / n_avg
    inds_arr = indgen(nvis_dims[0])
    if n_avg gt 1 then begin
      n_vis_freq_avg = fltarr(n_freq1)
      for i=0, n_freqbins-1 do begin
        inds_use = inds_arr[i*n_avg:i*n_avg+(n_avg-1)]
        if n_elements(inds_use) eq 1 then n_vis_freq_avg[i] = nvis_freq[inds_use] $
        else n_vis_freq_avg[i] = total(nvis_freq[inds_use])
      endfor
    endif else n_vis_freq_avg = nvis_freq
    cube1 = cube1 / rebin(reform(n_vis_freq_avg, 1, n_freq1), n_elements(hpx_inds1), n_freq1)
  endif
  if n_cubes gt 1 then begin
    if cube_varnames[1] eq '' then begin
      input_cube1 = getvar_savefile(filenames[1], input_cube_varnames[1, 0])
      input_cube2 = getvar_savefile(filenames[1], input_cube_varnames[1, 1])
      
      cube2 = temporary(input_cube1) - temporary(input_cube2)
    endif else cube2 = getvar_savefile(filenames[max_file], cube_varnames[1])
    n_freq2 = (size(cube2,/dimension))[1]
    if n_freq1 ne n_freq2 then message, 'number of frequencies do not match between the 2 files'
    if keyword_set(nvis_norm) then begin
      if obs_info.integrated[n_elements(obs_info.integrated)-1] eq 1 then obs_varname = 'obs_arr' else obs_varname = 'obs'
      obs_arr1 = getvar_savefile(filenames[max_file], obs_varname)
      nvis_freq = obs_arr1.nf_vis
      nvis_dims = size(nvis_freq, /dimension)
      if n_elements(nvis_dims) eq 2 then nvis_freq = total(nvis_freq, 2)
      
      n_avg = getvar_savefile(filenames[max_file], 'n_avg')
      n_freqbins = nvis_dims[0] / n_avg
      inds_arr = indgen(nvis_dims[0])
      if n_avg gt 1 then begin
        n_vis_freq_avg = fltarr(n_freq2)
        for i=0, n_freqbins-1 do begin
          inds_use = inds_arr[i*n_avg:i*n_avg+(n_avg-1)]
          if n_elements(inds_use) eq 1 then n_vis_freq_avg[i] = nvis_freq[inds_use] $
          else n_vis_freq_avg[i] = total(nvis_freq[inds_use])
        endfor
      endif else n_vis_freq_avg = nvis_freq
      if n_elements(filenames) eq 2 then cube2 = cube2 / rebin(reform(n_vis_freq_avg, 1, 1, n_freq2), n_elements(hpx_inds2), n_freq2) $
      else cube2 = cube2 / rebin(reform(n_vis_freq_avg, 1, 1, n_freq2), n_elements(hpx_inds1), n_freq2)
    endif
  endif
  
  print, 'nside, n pixels, n_freq: ' + number_formatter(nside1) + ', ' + number_formatter(n_elements(hpx_inds1)) + ', ' + number_formatter(n_freq1)
  
  case n_elements(slice_range) of
    0: begin
      slice_range = [0, n_freq1-1]
      if keyword_set(nvis_norm) then title_range = 'freq. averaged' else title_range = 'freq. added'
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
  
  if n_slice_range eq 2 then begin
    case n_elements(sr2) of
      1: begin
        title_range = [title_range, 'slice ' + number_formatter(sr2)]
      end
      2: begin
        if min(sr2) lt 0 then message, 'sr2 cannot be less than zero'
        if max(sr2) ge n_freq1 then message, 'sr2 cannot be more than ' + number_formatter(n_freq1-1)
        if sr2[1] lt sr2[0] then message, 'sr2[1] cannot be less than sr2[0]'
        
        title_range = [title_range, 'slices [' + number_formatter(sr2[0]) + ':' + number_formatter(sr2[1]) + ']']
      end
      else: begin
        message, 'sr2 must be a 1 or 2 element vector'
      end
    endcase
  endif
  
  if n_elements(title_range) eq 2 then range_str = title_range else range_str = strarr(2)
  
  ;; title to use:
  if n_cubes gt 1 then begin
    if keyword_set(ratio) then cube_op = '/' else cube_op = '-'
    
    if n_elements(folder_names) eq 1 then diff_title = evenodd[0] + '_' + cube_types[0] + '_' + pols[0] + '_' + range_str[0] + $
      cube_op + evenodd[max_eo] + '_' + cube_types[max_type] + '_' + pols[max_pol] + '_' + range_str[1] $
    else $
      diff_title = evenodd[0] + '_' + cube_types[0] + '_' + pols[0] + '_' + range_str[0] + $
      cube_op + evenodd[max_eo] + '_' + cube_types[max_type] + '_' + pols[max_pol] + '_' + range_str[1]
      
    if n_elements(title_range) eq 1 then diff_title = diff_title + ' ' + title_range
    
  endif else diff_title = evenodd[0] + '_' + cube_types[0] + '_' + pols[0] + ' ' + title_range
  
  if keyword_set(png) or keyword_set(eps) or keyword_set(pdf) then pub = 1 else pub = 0
  
  if n_cubes gt 1 then begin
    if n_slice_range eq 2 then begin
      if n_elements(slice_range) eq 1 then cube1 = cube1[*,slice_range] else cube1 = total(cube1[*, slice_range[0]:slice_range[1]],2)
      if n_elements(sr2) eq 1 then cube2 = cube2[*,sr2] else cube2 = total(cube2[*, sr2[0]:sr2[1]],2)
    endif else begin
      if n_elements(slice_range) eq 1 then begin
        cube1 = cube1[*,slice_range]
        cube2 = cube2[*,slice_range]
      endif else begin
        cube1 = total(cube1[*, slice_range[0]:slice_range[1]],2)
        cube2 = total(cube2[*, slice_range[0]:slice_range[1]],2)
      endelse
    endelse
    
    if max(abs(cube1-cube2)) eq 0 then message, 'cubes are identical.'
    if keyword_set(diff_ratio) then begin
      print, max(cube1), max(cube2), max(cube1)/max(cube2)
      temp = (cube1/max(cube1) - cube2/max(cube2)) * mean([max(cube1), max(cube2)])
      note = note + ', peak ratio = ' + number_formatter(max(cube1)/max(cube2), format = '(f5.2)')
    endif else if keyword_set(ratio) then temp = cube1/cube2 else temp = cube1-cube2
    
  endif else if n_elements(slice_range) eq 1 then temp = cube1[*,slice_range] else temp = total(cube1[*, slice_range[0]:slice_range[1]],2)
  
  if keyword_set(sym_color) and not keyword_set(log) then begin
    if n_elements(data_range) eq 0 then data_range = [-1,1]*max(abs(temp)) $
    else data_range = [-1,1]*max(abs(data_range))
  endif
  if keyword_set(diff_ratio) then title = diff_title + ', peak norm., ' else title = diff_title
  
  healpix_quickimage, temp, hpx_inds1, nside1, title = title, savefile = plotfile, note=note, slice_ind = slice_ind, $
    log = log, color_profile = color_profile, data_range = data_range, data_min_abs = data_min_abs,$
    window_num = window_num, plot_as_map = plot_as_map, png = png, eps = eps, pdf = pdf
    
end
