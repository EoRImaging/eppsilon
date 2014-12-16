pro ps_ratio_plots, folder_names, obs_info, cube_types, pols, $
    plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    note = note, spec_window_types = spec_window_types, data_range = data_range, data_min_abs = data_min_abs, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, diff_ratio = diff_ratio, $
    plot_wedge_line = plot_wedge_line, quiet = quiet, png = png, eps = eps, pdf = pdf, window_num = window_num
    
  if n_elements(obs_info.info_files) gt 2 then message, 'Only 1 or 2 info_files can be used'
  
  
  ;; default to blackman-harris spectral window
  if n_elements(spec_window_types) eq 0 then spec_window_types = 'Blackman-Harris'
  
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
  
  n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), n_elements(pols), n_elements(spec_window_types)])
  if n_diffs gt 2 then message, 'only 1 or 2 of [folder_names, obs_names, cube_types, pols, spec_window_types] allowed'
  
  if n_elements(cube_types) eq 0 then if n_diffs eq 1 then cube_types = ['res', 'dirty'] else cube_types = 'res'
  n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), n_elements(pols), n_elements(spec_window_types)])
  if n_elements(pols) eq 0 then if n_diffs eq 1 then pols=['xx', 'yy'] else pols = 'xx' 
  n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), n_elements(pols), n_elements(spec_window_types)])

  if n_diffs eq 1 then message, 'at least one of info_files, cube_types, pols, spec_window_types must be a 2 element vector'
  
  max_file = n_elements(obs_info.info_files)-1
  max_type = n_elements(cube_types)-1
  max_pol = n_elements(pols)-1
  max_sw = n_elements(spec_window_types)-1
  if max_sw lt 0 then max_sw=0
  
  ;; default to including baseline axis & delay axis
  if n_elements(baseline_axis) eq 0 then baseline_axis = 1
  if n_elements(delay_axis) eq 0 then delay_axis = 1
  
  ;; default to plot wedge line
  if n_elements(plot_wedge_line) eq 0 then plot_wedge_line=1
  
  ;; default to hinv
  if n_elements(hinv) eq 0 then hinv = 1
  
  
  if n_elements(folder_names) eq 2 then begin
    note = obs_info.diff_note
    if tag_exist(obs_info, 'diff_plot_path') then plot_path = obs_info.diff_plot_path else plot_path = obs_info.diff_save_path
    
    note = strjoin(strsplit(note, '-', /extract), '/')
    plot_path = strjoin(strsplit(plot_path, 'minus', /regex, /extract), 'over')
  endif else begin
    note = obs_info.fhd_types[0]
    plot_path = obs_info.plot_paths[0]
  endelse
  
  if n_elements(spec_window_types) eq 2 then note = note + ' ' + spec_window_types[0] + ' over ' + spec_window_types[1]
  
  if n_elements(folder_names) eq 1 then begin
    if n_elements(obs_info.obs_names) gt 1 then begin
      plot_filebase = obs_info.folder_basenames[0] + '_' + obs_info.obs_names[0] + '_' + cube_types[0] + '_' + pols[0] + sw_tags[0] + $
        '_over_' + obs_info.obs_names[0] + '_' + cube_types[max_type] + '_' + pols[max_pol] + sw_tags[max_sw]
    endif else begin
      if obs_info.integrated[0] eq 0 then plot_start = obs_info.folder_basenames[0] + '_' + obs_info.obs_names[0] else plot_start = obs_info.fhd_types[0]
      
      plot_filebase = plot_start + '_' + cube_types[0] + '_' + pols[0] + sw_tags[0] + $
        '_over_' + cube_types[max_type] + '_' + pols[max_pol] + sw_tags[max_sw]
    endelse
  endif else plot_filebase = obs_info.name_same_parts + '__' + strjoin([obs_info.name_diff_parts[0], cube_types[0], pols[0]], '_') + sw_tags[0] + $
    '_over_' + strjoin([obs_info.name_diff_parts[1], cube_types[max_type], pols[max_pol]], '_') + sw_tags[max_sw]
    
  if not file_test(plot_path, /directory) then file_mkdir, plot_path
  
  
  file_struct_arr1 = fhd_file_setup(obs_info.info_files[0], pol_inc, spec_window_type = spec_window_types[0])
  if n_elements(obs_info.info_files) eq 2 then file_struct_arr2 = fhd_file_setup(obs_info.info_files[1], pol_inc, spec_window_type = spec_window_types[max_sw]) $
  else file_struct_arr2 = file_struct_arr1
  
  type_pol_str = [cube_types[0] + '_' + pols[0], cube_types[max_type] + '_' + pols[max_pol]]
  
  if tag_exist(file_struct_arr1, 'n_obs') then n_obs1 = file_struct_arr1[0].n_obs
  if tag_exist(file_struct_arr2, 'n_obs') then n_obs2 = file_struct_arr2[0].n_obs
  
  if n_elements(n_obs1) gt 0 and n_elements(n_obs2) gt 0 then begin
    if n_elements(note) eq 0 then note = '(' + number_formatter(round(mean(file_struct_arr1[0].n_obs))) + ',' + number_formatter(round(mean(file_struct_arr2[0].n_obs))) + ')' $
    else note = note + ' (' + number_formatter(round(mean(file_struct_arr1[0].n_obs))) + ',' + number_formatter(round(mean(file_struct_arr2[0].n_obs))) + ')'
  endif
  
  
  type_pol_locs = intarr(2)
  wh_type_pol1 = where(file_struct_arr1.type_pol_str eq type_pol_str[0], count_type_pol)
  if count_type_pol eq 0 then $
    message, 'requested type_pol not found: ' + type_pol_str[0] + ' not in ' + file_struct_arr1.power_savefile else type_pol_locs[0] = wh_type_pol1
  wh_type_pol2 = where(file_struct_arr2.type_pol_str eq type_pol_str[1], count_type_pol)
  if count_type_pol eq 0 then $
    message, 'requested type_pol not found: ' + type_pol_str[1] + ' not in ' + file_struct_arr2.power_savefile else type_pol_locs[1] = wh_type_pol2
    
    
  if keyword_set(png) or keyword_set(eps) or keyword_set(pdf) then pub = 1 else pub = 0
  if pub then begin
    if keyword_set(png) and keyword_set(eps) and keyword_set(pdf) then begin
      print, 'only one of eps, pdf and png can be set, using png'
      eps = 0
    endif
    
    if keyword_set(png) then begin
      plot_exten = '.png'
      delete_ps = 1
    endif else if keyword_set(pdf) then begin
      plot_exten = '.pdf'
      delete_ps = 1
    endif else if keyword_set(eps) then begin
      plot_exten = '.eps'
      delete_ps = 0
    endif
  endif
  
  if keyword_set(plot_wedge_line) then begin
    z0_freq = 1420.40 ;; MHz
    redshifts = z0_freq/file_struct_arr1[0].frequencies - 1 ;; frequencies will be identical if kx, ky, kz match
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
  
  fadd_2dbin = ''
  ;;if keyword_set(fill_holes) then fadd_2dbin = fadd_2dbin + '_nohole'
  if keyword_set(no_kzero) then fadd_2dbin = fadd_2dbin + '_nok0'
  if keyword_set(log_kpar) then fadd_2dbin = fadd_2dbin + '_logkpar'
  if keyword_set(log_kperp) then fadd_2dbin = fadd_2dbin + '_logkperp'
  
  savefiles_2d = [file_struct_arr1[type_pol_locs[0]].savefile_froot + file_struct_arr1[type_pol_locs[0]].savefilebase + file_struct_arr1[0].power_tag, $
    file_struct_arr2[type_pol_locs[1]].savefile_froot + file_struct_arr2[type_pol_locs[1]].savefilebase + file_struct_arr2[0].power_tag] + $
    fadd_2dbin + '_2dkpower.idlsave'
  test_save_2d = file_test(savefiles_2d) *  (1 - file_test(savefiles_2d, /zero_length))
  
  if min(test_save_2d) eq 0 then message, '2D savefile not found: ' + savefiles_2d[where(test_save_2d eq 0)]
  
  title = type_pol_str[0] + '/' + type_pol_str[1]
  
  kpower_2d_plots, savefiles_2d, png = png, eps = eps, pdf = pdf, plotfile = plotfile, window_num = window_num, $
    kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
    data_range = data_range, title_prefix = title, note = note, $
    plot_wedge_line = plot_wedge_line, hinv = hinv, /power_ratio, $
    wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
    
end