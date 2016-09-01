pro compare_plot_prep, folder_names, obs_info, cube_types, pols, comp_type, compare_files, $
    plot_slices = plot_slices, slice_type = slice_type, fadd_2dbin = fadd_2dbin, $
    spec_window_types = spec_window_types, delta_uv_lambda = delta_uv_lambda, freq_ch_range = freq_ch_range, $
    ave_removal = ave_removal, $
    baseline_axis = baseline_axis, delay_axis = delay_axis, plot_wedge_line = plot_wedge_line, hinv = hinv, $
    plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    axis_type_1d = axis_type_1d, note = note, wt_cutoffs = wt_cutoffs, wt_measures = wt_measures, $
    pub = pub, plot_exten = plot_exten, full_compare = full_compare
    
    
  comp_type_enum = ['diff', 'diff_ratio', 'ratio']
  wh_comp_type = where(comp_type_enum eq comp_type, count_comp_type)
  if count_comp_type eq 0 then message, 'comp_type not recognized, must be one of: ' + strjoin(comp_type_enum, ' ,')
  if comp_type eq 'ratio' and keyword_set(full_compare) then comp_type = 'diff_ratio'
  
  
  if n_elements(obs_info.info_files) gt 2 then message, 'Only 1 or 2 info_files can be used'
  
  if n_elements(delta_uv_lambda) gt 2 then message, 'only 1 delta_uv_lambda allowed'
  
  ;; default to blackman-harris spectral window
  if n_elements(spec_window_types) eq 0 then spec_window_types = 'Blackman-Harris'
  
  ;; default to ave_removal
  if n_elements(ave_removal) eq 0 then ave_removal = 1
   
  if n_elements(spec_window_types) eq 2 then begin
    if spec_window_types[0] eq spec_window_types[1] then spec_window_types = spec_window_types[0] else begin
    
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
      
    endelse
  endif else sw_tags = ''
  
  ;; density correction defaults & file naming for 2D & 1D files
  if n_elements(wt_cutoffs) eq 0 then begin
    ;; default to wt_cutoffs = 1, wt_measures = 'min'
    wt_cutoffs = 1
    wt_measures = 'min'
    n_wtcuts = 1
  endif else begin
    n_wtcuts = max([n_elements(wt_cutoffs), n_elements(wt_measures)])
    if comp_type eq 'diff' then if n_wtcuts gt 1 then message, 'no more than one wt_cutoff and wt_measure can be set for difference plots' $
    else if n_wtcuts gt 2 then message, 'no more than 2 of wt_cutoffs and wt_measures can be set'
    
    if n_elements(wt_measures) eq 0 then begin
      print, 'wt_cutoffs is specified but wt_measures is not. Defaulting wt_measures to "min".'
      wt_measures = 'min'
    endif
    
    if n_wtcuts eq 2 then begin
      if n_elements(wt_cutoffs) eq 2 then if wt_cutoffs[1] - wt_cutoffs[0] le 1e-3 then wt_cutoffs = wt_cutoffs[0]
      if n_elements(wt_measures) eq 2 then if wt_measures[0] eq wt_measures[1] then wt_measures = wt_measures[0]
      n_wtcuts = max([n_elements(wt_cutoffs), n_elements(wt_measures)])
      
      if n_elements(wt_cutoffs) lt n_wtcuts then wt_cutoffs = [wt_cutoffs, wt_cutoffs]
      if n_elements(wt_measures) lt n_wtcuts then wt_measures = [wt_measures, wt_measures]
    endif
    
  endelse
  kperp_density_names = strarr(n_wtcuts)
  wh_cutoff0 = where(wt_cutoffs eq 0, count_cutoff0, complement = wh_cutoff_n0, ncomplement = count_cutoff_n0)
  wh_std = where(wt_cutoffs eq 1 and wt_measures eq 'min', count_std)
  
  if count_cutoff0 gt 0 then kperp_density_names[wh_cutoff0] = '_nodensitycorr'
  if count_cutoff_n0 gt 0 then kperp_density_names[wh_cutoff_n0] = '_kperp_density_' + wt_measures[wh_cutoff_n0] + '_gt' + number_formatter(wt_cutoffs[wh_cutoff_n0])
  if count_std gt 0 then kperp_density_names[wh_std] = '_dencorr'
  
  if n_elements(freq_ch_range) ne 0 then begin
    if min(freq_ch_range) lt 0 or max(freq_ch_range) - min(freq_ch_range) lt 3 then message, 'invalid freq_ch_range'
    fch_tag = '_ch' + number_formatter(min(freq_ch_range)) + '-' + number_formatter(max(freq_ch_range))
  endif else fch_tag = ''
  
  if n_elements(ave_removal) eq 2 then if ave_removal[1] eq ave_removal[0] then ave_removal = ave_removal[0]
  if n_elements(ave_removal) eq 0 then ave_removal = 0
  
  n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), n_elements(pols), n_elements(spec_window_types), n_wtcuts, n_elements(ave_removal)])
  if n_diffs gt 2 then message, 'only 1 or 2 of [folder_names, obs_names, cube_types, pols, spec_window_types, wt_cutoffs/measures, ave_removal] allowed'
  
  if n_elements(obs_info.info_files) eq 2 or n_elements(spec_window_types) eq 2 or n_wtcuts eq 2 or n_elements(ave_removal) eq 2 $
    and n_elements(cube_types) eq 0 and n_elements(pols) eq 0 and n_elements(full_compare) eq 0 then full_compare = 1
    
  if keyword_set(full_compare) and n_elements(info_files) eq 1 and n_elements(spec_window_types) lt 2 and n_wtcuts eq 1 and n_elements(ave_removal) lt 2 then $
    message, 'full_compare can only be set with 2 folder names and/or 2 obs names and/or 2 spec windows and/or 2 wt_cutoffs/measures and/or 2 ave_removal values'
  if comp_type eq 'ratio' and keyword_set(full_compare) then comp_type = 'diff_ratio'
  
  
  if not keyword_set(full_compare) then begin
    case comp_type of
      'diff': begin
        if n_elements(cube_types) eq 0 then if n_diffs eq 1 then cube_types = ['dirty', 'res'] else cube_types = 'res'
        n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), n_elements(pols), n_elements(spec_window_types), n_elements(ave_removal)])
        if n_elements(pols) eq 0 then if n_diffs eq 1 then pols=['xx', 'yy'] else pols = 'xx'
        n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), n_elements(pols), n_elements(spec_window_types), n_elements(ave_removal)])
        
        if n_diffs eq 1 then message, 'at least one of info_files, cube_types, pols, spec_window_types, ave_removal must be a 2 element vector'
        
        n_cubes = 1
      end
      'ratio': begin
        n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), n_elements(pols), $
          n_elements(spec_window_types), n_elements(wt_cutoffs), n_elements(ave_removal)])
        if n_diffs gt 2 then message, 'only 1 or 2 of [folder_names, obs_names, cube_types, pols, spec_window_types, wt_cutoffs, ave_removal] allowed'
        
        if n_elements(cube_types) eq 0 then if n_diffs eq 1 then cube_types = ['res', 'dirty'] else cube_types = 'res'
        n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), n_elements(pols), $
          n_elements(spec_window_types), n_elements(wt_cutoffs), n_elements(ave_removal)])
        if n_elements(pols) eq 0 then if n_diffs eq 1 then pols=['xx', 'yy'] else pols = 'xx'
        n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), n_elements(pols), $
          n_elements(spec_window_types), n_elements(wt_cutoffs), n_elements(ave_removal)])
          
        if n_diffs eq 1 then message, 'at least one of info_files, cube_types, pols, spec_window_types, wt_cutoffs, ave_removal must be a 2 element vector'
        
        n_sets=1
      end
      'diff_ratio': begin
        if n_elements(obs_info.info_files) ne 2 and n_elements(spec_window_types) lt 2 $
          and n_elements(pols) lt 2 and n_elements(wt_cutoffs) lt 2 and n_elements(ave_removal) lt 2 then $
          message, 'diff_ratio requires 2 folder names and/or 2 obs names and/or 2 spec windows and/or 2 wt_cutoffs and/or 2 pols and/or 2 ave_removal values'
          
        if n_elements(pols) eq 0 then pols = 'xx'
        
        n_sets=2
        
        if n_elements(cube_types) ne 2 then cube_types = ['res', 'dirty']
      end
    endcase
  endif else begin
    case comp_type of
      'diff': undefine, cube_types, pols
      'diff_ratio': begin
        pols = ['xx', 'yy']
        if n_elements(cube_types) ne 2 then cube_types = ['res', 'dirty']
        
        n_sets=4
      end
    endcase
  endelse
  
  
  max_file = n_elements(obs_info.info_files)-1
  max_type = n_elements(cube_types)-1
  max_pol = n_elements(pols)-1
  max_sw = n_elements(spec_window_types)-1
  if max_sw lt 0 then max_sw=0
  max_wtcut = n_wtcuts-1
  max_ar = n_elements(ave_removal)-1
  if max_ar lt 0 then max_ar=0
  
  ;; default to including baseline axis & delay axis
  if n_elements(baseline_axis) eq 0 then baseline_axis = 1
  if n_elements(delay_axis) eq 0 then delay_axis = 1
  
  ;; default to plot wedge line
  if n_elements(plot_wedge_line) eq 0 then plot_wedge_line=1
  
  ;; default to hinv
  if n_elements(hinv) eq 0 then hinv = 1
  
  if n_elements(axis_type_1d) eq 0 then axis_type_1d = 'sym_log'
  
  if n_elements(folder_names) eq 2 then begin
    if n_elements(save_path) eq 0 then save_path = obs_info.diff_save_path
    note = obs_info.diff_note
    if n_elements(plot_path) eq 0 then if tag_exist(obs_info, 'diff_plot_path') then $
      plot_path = obs_info.diff_plot_path else plot_path = save_path + path_sep() + 'plots' + path_sep()
  endif else begin
    if n_elements(save_path) eq 0 then save_path = obs_info.save_paths[0]
    note = obs_info.fhd_types[0]
    if n_elements(plot_path) eq 0 then plot_path = obs_info.plot_paths[0]
  endelse
  
  if n_elements(spec_window_types) eq 2 then note = note + ' ' + spec_window_types[0] + ' minus ' + spec_window_types[1]
  if n_elements(ave_removal) eq 2 then begin
    ave_removal_str = ['no_ave_removal', 'ave_removal']
    note = note + ' ' + ave_removal_str[ave_removal[0]] + ' minus ' + ave_removal_str[ave_removal[1]]
    ar_tag_list = ['', '_averemoval']
    ar_tags = [ar_tag_list[ave_removal[0]], ar_tag_list[ave_removal[1]]]
  endif else ar_tags = ''
  if n_wtcuts eq 2 then density_tags = kperp_density_names else density_tags = ''
  
  
  file_struct_arr1 = fhd_file_setup(obs_info.info_files[0], $
    spec_window_type = spec_window_types[0], delta_uv_lambda = delta_uv_lambda, freq_ch_range = freq_ch_range, ave_removal = ave_removal[0])
  if n_elements(obs_info.info_files) eq 2 or max_sw eq 1 or max_wtcut eq 1 or max_ar eq 1 then $
    file_struct_arr2 = fhd_file_setup(obs_info.info_files[max_file], spec_window_type = spec_window_types[max_sw], $
    delta_uv_lambda = delta_uv_lambda, freq_ch_range = freq_ch_range, ave_removal = ave_removal[max_ar]) $
  else file_struct_arr2 = file_struct_arr1
  type_pol_str1 = file_struct_arr1.type_pol_str
  type_pol_str2 = file_struct_arr2.type_pol_str
  
  if keyword_set(plot_slices) then begin
    if n_elements(slice_type) eq 0 then slice_type = 'kspace'
    slice_type_enum = ['raw', 'divided', 'kspace', 'sumdiff', 'weights']
    
    wh_slice_type = where(slice_type_enum eq slice_type, count_slice_type)
    if count_slice_type eq 0 then begin
      print, 'slice_type not recognized, using default'
      slice_type = 'kspace'
    endif
    
    if slice_type ne 'kspace' then message, 'only kspace slice difference plots are currently supported'
    
    slice_tags = ['xz', 'yz', 'xy']
    n_slices = n_elements(slice_tags)
  endif else n_slices = 1
  
  ;; get same parts of power_tag to add to plot file name
  if n_elements(file_struct_arr2) eq 0 then same_power_tag = file_struct_arr1[0].power_tag else begin
    if file_struct_arr1[0].power_tag eq file_struct_arr2[0].power_tag then same_power_tag = file_struct_arr1[0].power_tag else begin
      tag_arr1 = strsplit(file_struct_arr1[0].power_tag, '_',/extract)
      tag_arr2 = strsplit(file_struct_arr2[0].power_tag, '_',/extract)
      for i=0, n_elements(tag_arr1)-1 do begin
        wh_in2 = where(tag_arr2 eq tag_arr1[i], count_in2)
        if count_in2 gt 0 then begin
          if n_elements(same_arr) eq 0 then same_arr = tag_arr1[i] else same_arr = [same_arr, tag_arr1[i]]
        endif
      endfor
      if n_elements(same_arr) gt 0 then same_power_tag = '_' + strjoin(same_arr, '_')
    endelse
  endelse
  
  if max_wtcut eq 0 then same_density_tag = kperp_density_names else same_density_tag = ''
  
  case comp_type of
    'diff': op_str = '_minus_'
    'diff_ratio': op_str = '_diffratio_'
    'ratio': op_str = '_over_'
  endcase
  
  if keyword_set(full_compare) then begin
  
    if n_elements(folder_names) eq 2 and folder_names[0] ne folder_names[n_elements(folder_names)-1] then begin
      plot_filebase = obs_info.name_same_parts + same_power_tag + same_density_tag +'__' + obs_info.name_diff_parts[0] + '_' + obs_info.obs_names[0] + sw_tags[0] + ar_tags[0] + density_tags[0] + $
        op_str + obs_info.name_diff_parts[1]  + '_' + obs_info.obs_names[1] + sw_tags[max_sw] + ar_tags[max_ar] + density_tags[max_wtcut]
    endif else plot_filebase = obs_info.folder_basenames[0] + same_power_tag + same_density_tag + '__' + obs_info.obs_names[0] + sw_tags[0] + ar_tags[0] + density_tags[0] + op_str + $
      obs_info.obs_names[max_file] + sw_tags[max_sw] + ar_tags[max_ar] + density_tags[max_wtcut]
      
  endif else begin
  
    if n_elements(folder_names) eq 1 then begin
      if n_elements(obs_info.obs_names) gt 1 then begin
        plot_filebase = obs_info.folder_basenames[0] + same_power_tag + same_density_tag + '_' + obs_info.obs_names[0] + '_' + cube_types[0] + '_' + pols[0] + sw_tags[0] + ar_tags[0] + density_tags[0] + $
          op_str + obs_info.obs_names[0] + '_' + cube_types[max_type] + '_' + pols[max_pol] + sw_tags[max_sw] + ar_tags[max_ar] + density_tags[max_wtcut]
      endif else begin
        if obs_info.integrated[0] eq 0 then plot_start = obs_info.folder_basenames[0] + '_' + obs_info.obs_names[0] else plot_start = obs_info.fhd_types[0]
        
        plot_filebase = plot_start + same_power_tag + same_density_tag + '_' + cube_types[0] + '_' + pols[0] + sw_tags[0] + ar_tags[0] + density_tags[0] + $
          op_str + cube_types[max_type] + '_' + pols[max_pol] + sw_tags[max_sw] + ar_tags[max_ar] + density_tags[max_wtcut]
      endelse
    endif else plot_filebase = obs_info.name_same_parts + same_power_tag + same_density_tag + '__' + strjoin([obs_info.name_diff_parts[0], cube_types[0], pols[0]], '_') + sw_tags[0] + ar_tags[0] + density_tags[0] + $
      op_str + strjoin([obs_info.name_diff_parts[1], cube_types[max_type], pols[max_pol]], '_') + sw_tags[max_sw] + ar_tags[max_ar] + density_tags[max_wtcut]
  endelse
  plot_filebase = plot_filebase + fch_tag
  
  if not file_test(plot_path, /directory) then file_mkdir, plot_path
  if not file_test(save_path, /directory) then file_mkdir, save_path
  
  if keyword_set(full_compare) then begin
    case comp_type of
      'diff': begin
        if n_elements(type_pol_str1) ne n_elements(type_pol_str2) then message, 'all_type_pol cannot be used with these folders, they contain different number of types & pols'
        n_cubes = n_elements(type_pol_str1)
        
        for i=0, n_cubes-1 do begin
          temp = where(type_pol_str2 eq type_pol_str1[i], count_typepol)
          if count_typepol eq 0 then message, 'all_type_pol cannot be used with these folders, they contain different sets of types & pols'
        endfor
        
        type_pol_str = [[type_pol_str1],[type_pol_str1]]
        
        titles = type_pol_str[*, 0]
      end
      'diff_ratio': begin
        type_pol_str = strarr(n_sets, 2)
        type_pol_str[0, *] = cube_types[0] + '_' + pols[0]
        type_pol_str[1, *] = cube_types[1] + '_' + pols[0]
        type_pol_str[2, *] = cube_types[0] + '_' + pols[max_pol]
        type_pol_str[3, *] = cube_types[1] + '_' + pols[max_pol]
        
        titles = strarr(2,3)
        for i=0, 1 do titles[i,*] = [type_pol_str[2*i,0] + '/' + type_pol_str[2*i+1,0], $
          type_pol_str[2*i,1] + '/' + type_pol_str[2*i+1,1], 'Ratio Difference']
      end
    endcase
  endif else begin
    case comp_type of
      'diff': begin
        type_pol_str = strarr(1, 2)
        type_pol_str[0]= cube_types[0] + '_' + pols[0]
        type_pol_str[1] = cube_types[max_type] + '_' + pols[max_pol]
        
        titles = type_pol_str[0] + '-' + type_pol_str[1]
      end
      'diff_ratio': begin
        type_pol_str = strarr(n_sets, 2)
        
        type_pol_str[0, *] = cube_types + '_' + pols[0]
        type_pol_str[1, *] = cube_types + '_' + pols[max_pol]
        
        titles = [type_pol_str[0,0] + '/' + type_pol_str[0,1], $
          type_pol_str[1,0] + '/' + type_pol_str[1,1], 'Ratio Difference']
      end
      'ratio': begin
        type_pol_str = strarr(1, 2)
        type_pol_str[0]= cube_types[0] + '_' + pols[0]
        type_pol_str[1] = cube_types[max_type] + '_' + pols[max_pol]
        
        titles = type_pol_str[0] + '/' + type_pol_str[1]
      end
    endcase
  endelse
  
  if comp_type eq 'diff' then cube_set_dim = n_cubes else cube_set_dim = n_sets
  type_pol_locs = intarr(cube_set_dim, 2)
  for i=0, cube_set_dim-1 do begin
    wh_type_pol1 = where(file_struct_arr1.type_pol_str eq type_pol_str[i,0], count_type_pol)
    if count_type_pol eq 0 then $
      message, 'requested type_pol not found: ' + type_pol_str[i,0] + ' not in ' + folder_names[0] else type_pol_locs[i,0] = wh_type_pol1
    wh_type_pol2 = where(file_struct_arr2.type_pol_str eq type_pol_str[i,1], count_type_pol)
    if count_type_pol eq 0 then $
      message, 'requested type_pol not found: ' + type_pol_str[i,1] + ' not in ' + folder_names[1] else type_pol_locs[i,1] = wh_type_pol2
  endfor
  
  if tag_exist(file_struct_arr1, 'n_obs') then n_obs1 = file_struct_arr1[0].n_obs
  if tag_exist(file_struct_arr2, 'n_obs') then n_obs2 = file_struct_arr2[0].n_obs
  
  if n_elements(n_obs1) gt 0 and n_elements(n_obs2) gt 0 then begin
    if n_elements(note) eq 0 then note = '(' + number_formatter(round(mean(file_struct_arr1[0].n_obs))) + ',' + number_formatter(round(mean(file_struct_arr2[0].n_obs))) + ')' $
    else note = note + ' (' + number_formatter(round(mean(file_struct_arr1[0].n_obs))) + ',' + number_formatter(round(mean(file_struct_arr2[0].n_obs))) + ')'
  endif
  
  
  if keyword_set(plot_wedge_line) then begin
    z0_freq = 1420.40 ;; MHz
    
    freq_use = file_struct_arr1[0].frequencies
    if n_elements(freq_ch_range) ne 0 then begin
      if max(freq_ch_range) gt n_elements(freq_use)-1 then message, 'invalid freq_ch_range'
      freq_use = freq_use[freq_ch_range[0]:freq_ch_range[1]]
    endif
    
    redshifts = z0_freq/freq_use - 1 ;; frequencies will be identical if kx, ky, kz match
    mean_redshift = mean(redshifts)
    
    cosmology_measures, mean_redshift, wedge_factor = wedge_factor
    ;; assume 20 degrees from pointing center to first null
    source_dist = 20d * !dpi / 180d
    fov_amp = wedge_factor * source_dist
    
    ;; calculate angular distance to horizon
    max_theta = max([file_struct_arr1[0].max_theta, file_struct_arr2[0].max_theta])
    horizon_amp = wedge_factor * ((max_theta+90d) * !dpi / 180d)
    
    wedge_amp = [fov_amp, horizon_amp]
    
    wedge_1dbin_names = ['', '_no_fov_wedge', '_no_horizon_wedge']
  endif else begin
    wedge_amp = 0d
    wedge_1dbin_names = ''
  endelse
  
  
  if comp_type eq 'diff' then begin
    mid_savefile_2d = strarr(n_slices, cube_set_dim)
    if not keyword_set(plot_slices) then begin
      savefiles_1d = strarr(cube_set_dim, n_elements(wedge_1dbin_names))
      mid_savefile_3d = strarr(n_slices, cube_set_dim)
    endif
  endif
  input_savefile1 = strarr(n_slices, cube_set_dim)
  input_savefile2 = strarr(n_slices, cube_set_dim)
  
  if keyword_set(pub) then plotfiles_2d = strarr(n_slices)
  
  for slice_i=0, n_slices-1 do begin
    for cube_i=0, cube_set_dim-1 do begin
    
      type_pol1 = type_pol_str[cube_i, 0]
      type_pol2 = type_pol_str[cube_i, 1]
      
      if comp_type eq 'diff' then begin
      
        if n_elements(obs_info.info_files) eq 2 then begin
          fileparts_1 = strsplit(file_struct_arr1[0].general_filebase, '_', /extract)
          fileparts_2 = strsplit(file_struct_arr2[0].general_filebase, '_', /extract)
          match_test = strcmp(fileparts_1, fileparts_2)
          wh_diff = where(match_test eq 0, count_diff, complement = wh_same, ncomplement = count_same)
        endif
        
        if n_elements(savefilebase) eq 0 then begin
          if n_elements(obs_info.info_files) eq 1 then begin
            savefilebase_use = file_struct_arr1[0].general_filebase + '_' + type_pol1 + op_str + type_pol2
          endif else begin
            if count_diff eq 0 then begin
              savefilebase_use = file_struct_arr1[0].general_filebase + '_' + type_pol1
              if type_pol1 ne type_pol2 then savefilebase_use = savefilebase_use + '_' + type_pol1 + op_str + type_pol2
            endif else begin
              if count_same gt 0 then begin
                if type_pol1 ne type_pol2 then savefilebase_use = strjoin(fileparts_1[wh_same], '_') + '__' + $
                  strjoin(fileparts_1[wh_diff]) + '_' + type_pol1 + op_str + strjoin(fileparts_2[wh_diff]) + '_' + type_pol2 $
                else savefilebase_use = strjoin(fileparts_1[wh_same], '_') + '__' + strjoin(fileparts_1[wh_diff]) + op_str + strjoin(fileparts_2[wh_diff]) + '__' + type_pol1
              endif else begin
                if type_pol1 ne type_pol2 then savefilebase_use = file_struct_arr1[0].general_filebase + '_' + type_pol_str[0] + $
                  op_str + file_struct_arr2[0].general_filebase + '_' + type_pol_str[1] $
                else savefilebase_use = file_struct_arr1[0].general_filebase + op_str + file_struct_arr2[0].general_filebase + '__' + type_pol1
              endelse
            endelse
          endelse
        endif else begin
          if type_pol1 eq type_pol2 then savefilebase_use = savefilebase + '_' + type_pol1 else savefilebase_use = savefilebase + '_' + type_pol1 + op_str + type_pol2
        endelse
      endif
      
      if keyword_set(plot_slices) then begin
        if comp_type eq 'diff' then mid_savefile_2d[slice_i, cube_i] = save_path + savefilebase_use + kperp_density_names + '_power_' + slice_tags[slice_i] + '_plane.idlsave'
        
        slice_filebase1 = file_struct_arr1[type_pol_locs[cube_i, 0]].savefile_froot + file_struct_arr1[type_pol_locs[cube_i, 0]].savefilebase + file_struct_arr1[type_pol_locs[cube_i, 0]].power_tag
        input_savefile1[slice_i, cube_i] = slice_filebase1 + '_' + slice_tags[slice_i] + '_plane.idlsave'
        
        slice_filebase2 = file_struct_arr2[type_pol_locs[cube_i, 1]].savefile_froot + file_struct_arr2[type_pol_locs[cube_i, 1]].savefilebase + file_struct_arr2[type_pol_locs[cube_i, 1]].power_tag
        input_savefile2[slice_i, cube_i] = slice_filebase2 + '_' + slice_tags[slice_i] + '_plane.idlsave'
        
      endif else begin
        if comp_type eq 'diff' then begin
          mid_savefile_3d[slice_i, cube_i] = save_path + savefilebase_use + '_power.idlsave'
          mid_savefile_2d[slice_i, cube_i] = save_path + savefilebase_use + kperp_density_names + '_2dkpower.idlsave'
          
          for wedge_i=0, n_elements(wedge_1dbin_names)-1 do begin
            savefiles_1d[cube_i, wedge_i] = save_path + savefilebase_use + kperp_density_names + wedge_1dbin_names[wedge_i] + '_1dkpower.idlsave'
          endfor
          
          input_savefile1[slice_i, cube_i] = file_struct_arr1[type_pol_locs[cube_i, 0]].power_savefile
          if file_test(input_savefile1[slice_i, cube_i]) eq 0 then message, 'No power file for ' + type_pol1 + ' and info_file: ' + obs_info.info_files[0]
          
          input_savefile2[slice_i, cube_i] = file_struct_arr2[type_pol_locs[cube_i, 1]].power_savefile
          if file_test(input_savefile2[slice_i, cube_i]) eq 0 then message, 'No power file for ' + type_pol2 + ' and info_file: ' + obs_info.info_files[n_elements(obs_info.info_files)-1]
          
        endif else begin
        
          input_savefile1[slice_i, cube_i] = file_struct_arr1[type_pol_locs[cube_i, 0]].savefile_froot + file_struct_arr1[type_pol_locs[cube_i, 0]].savefilebase + $
            file_struct_arr1[type_pol_locs[cube_i, 0]].power_tag + fadd_2dbin + kperp_density_names[0] + '_2dkpower.idlsave'
          input_savefile2[slice_i, cube_i] = file_struct_arr2[type_pol_locs[cube_i, 1]].savefile_froot + file_struct_arr2[type_pol_locs[cube_i, 1]].savefilebase + $
            file_struct_arr2[type_pol_locs[cube_i, 1]].power_tag + fadd_2dbin + kperp_density_names[max_wtcut] + '_2dkpower.idlsave'
            
        endelse
      endelse
      
      if pub then begin
        if keyword_set(plot_slices) then begin
          plotfiles_2d[slice_i] = plot_path + plot_filebase + '_power_' + slice_tags[slice_i] + '_plane' + plot_exten
        endif else begin
        
          plotfiles_2d[slice_i] = plot_path + plot_filebase + '_2dkpower' + plot_exten
          
          if comp_type eq 'diff' and cube_i eq 0 then plotfiles_1d = plot_path + plot_filebase + ['', wedge_1dbin_names] + '_1dkpower' + plot_exten
        endelse
      endif else plot_exten = ''
      
      kperp_lambda_conv = getvar_savefile(input_savefile1[slice_i, cube_i], 'kperp_lambda_conv')
      hubble_param = getvar_savefile(input_savefile1[slice_i, cube_i], 'hubble_param')
      kperp_plot_range = [5./kperp_lambda_conv, min([file_struct_arr1.kspan/2.,file_struct_arr1.max_baseline_lambda])/kperp_lambda_conv]
      
      if keyword_set(hinv) then kperp_plot_range = kperp_plot_range / hubble_param
      
    endfor
  endfor
  
  compare_files = {input_savefile1:input_savefile1, input_savefile2:input_savefile2, $
    titles:titles, n_slices:n_slices, $
    wedge_amp:wedge_amp, kperp_density_names:kperp_density_names, kperp_plot_range:kperp_plot_range}
    
  if n_elements(n_cubes) gt 0 then compare_files = create_struct(compare_files, 'n_cubes', n_cubes)
  if n_elements(n_sets) gt 0 then compare_files = create_struct(compare_files, 'n_sets', n_sets)
  
  if keyword_set(pub) then compare_files = create_struct(compare_files, 'plotfiles_2d', plotfiles_2d)
  if n_elements(mid_savefile_2d) gt 0 then compare_files = create_struct(compare_files, 'mid_savefile_2d', mid_savefile_2d)
  if n_elements(mid_savefile_3d) gt 0 then compare_files = create_struct(compare_files, 'mid_savefile_3d', mid_savefile_3d)
  if n_elements(savefiles_1d) gt 0 then compare_files = create_struct(compare_files, 'savefiles_1d', savefiles_1d)
  if n_elements(slice_tags) gt 0 then compare_files = create_struct(compare_files, 'slice_tags', slice_tags)
  
end