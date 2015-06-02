pro ps_ratio_plots, folder_names, obs_info, cube_types, pols, all_pol_diff_ratio = all_pol_diff_ratio, $
    freq_ch_range = freq_ch_range, plot_path = plot_path, plot_filebase = plot_filebase, $
    save_path = save_path, savefilebase = savefilebase, $
    note = note, spec_window_types = spec_window_types, data_range = data_range, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, $
    diff_ratio = diff_ratio, diff_range = diff_range, diff_min_abs = diff_min_abs, $
    plot_wedge_line = plot_wedge_line, quiet = quiet, png = png, eps = eps, pdf = pdf, $
    window_num = window_num
    
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
  
  if n_elements(obs_info.info_files) eq 2 or n_elements(spec_window_types) eq 2 $
    and n_elements(cube_types) eq 0 and n_elements(pols) eq 0 $
    and n_elements(all_pol_diff_ratio) eq 0 and n_elements(diff_ratio) eq 0 then all_pol_diff_ratio = 1
    
  if keyword_set(diff_ratio) or keyword_set(all_pol_diff_ratio) then begin
  
    if keyword_set(all_pol_diff_ratio) then begin
      if  n_elements(obs_info.info_files) ne 2 and n_elements(spec_window_types) lt 2 then $
        message, 'all_pol_diff_ratio requires 2 folder names and/or 2 obs names and/or 2 spec windows'
        
      pols = ['xx', 'yy']
      
      n_sets=4
    endif else begin
      if n_elements(obs_info.info_files) ne 2 and n_elements(spec_window_types) lt 2 and n_elements(pols) lt 2 then $
        message, 'diff_ratio requires 2 folder names and/or 2 obs names and/or 2 spec windows and/or 2 pols'
        
      if n_elements(pols) eq 0 then pols = 'xx'
      
      n_sets=2
    endelse
    
    if n_elements(cube_types) ne 2 then cube_types = ['res', 'dirty']
    
  endif else begin
    n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), n_elements(pols), n_elements(spec_window_types)])
    if n_diffs gt 2 then message, 'only 1 or 2 of [folder_names, obs_names, cube_types, pols, spec_window_types] allowed'
    
    if n_elements(cube_types) eq 0 then if n_diffs eq 1 then cube_types = ['res', 'dirty'] else cube_types = 'res'
    n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), n_elements(pols), n_elements(spec_window_types)])
    if n_elements(pols) eq 0 then if n_diffs eq 1 then pols=['xx', 'yy'] else pols = 'xx'
    n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), n_elements(pols), n_elements(spec_window_types)])
    
    if n_diffs eq 1 then message, 'at least one of info_files, cube_types, pols, spec_window_types must be a 2 element vector'
    
    n_sets=1
  endelse
  
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
    if n_elements(save_path) eq 0 then save_path = obs_info.diff_save_path
    note = obs_info.diff_note
    if tag_exist(obs_info, 'diff_plot_path') then plot_path = obs_info.diff_plot_path else plot_path = save_path
    
    if not keyword_set(diff_ratio) then note = strjoin(strsplit(note, '-', /extract), '/')
    
  endif else begin
    note = obs_info.fhd_types[0]
    plot_path = obs_info.plot_paths[0]
  endelse
  
  if n_elements(spec_window_types) eq 2 then note = note + ' ' + spec_window_types[0] + ' over ' + spec_window_types[1]
  
  if keyword_set(diff_ratio) or keyword_set(all_pol_diff_ratio) then begin
  
    if keyword_set(all_pol_diff_ratio) then begin
      if n_elements(folder_names) eq 2 and folder_names[0] ne folder_names[1] then begin
        plot_filebase = obs_info.name_same_parts + '__' + obs_info.name_diff_parts[0] + '_' + obs_info.obs_names[0] + sw_tags[0] + $
          '_diffratio_' + obs_info.name_diff_parts[1]  + '_' + obs_info.obs_names[1] + sw_tags[max_sw]
      endif else if n_elements(obs_info.info_files) eq 2 then $
        plot_filebase = obs_info.folder_basenames[0] + '__' + obs_info.obs_names[0] + sw_tags[0] + '_diffratio_' + $
        obs_info.obs_names[1] + sw_tags[max_sw] $
      else plot_filebase = obs_info.fhd_types[0] + sw_tags[0] + '_diffratio_' + sw_tags[max_sw]
    endif else begin
      if n_elements(folder_names) eq 2 and folder_names[0] ne folder_names[1] then begin
        plot_filebase = obs_info.name_same_parts + '__' + obs_info.name_diff_parts[0] + '_' + obs_info.obs_names[0] + pols[0] + sw_tags[0] + $
          '_diffratio_' + obs_info.name_diff_parts[1]  + '_' + obs_info.obs_names[1] + pols[max_pol] + sw_tags[max_sw]
      endif else if n_elements(obs_info.info_files) eq 2 then $
        plot_filebase = obs_info.folder_basenames[0] + '__' + obs_info.obs_names[0] + pols[0] + sw_tags[0] + '_diffratio_' + $
        obs_info.obs_names[1] + pols[max_pol] + sw_tags[max_sw] $
      else plot_filebase = obs_info.fhd_types[0] + pols[0] + sw_tags[0] + '_diffratio_' + pols[max_pol] + sw_tags[max_sw]
    endelse
  endif else begin
  
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
  endelse
  
  if keyword_set(pub) and not file_test(plot_path, /directory) then file_mkdir, plot_path
  
  
  file_struct_arr1 = fhd_file_setup(obs_info.info_files[0], pol_inc, $
    spec_window_type = spec_window_types[0], freq_ch_range = freq_ch_range)
  if n_elements(obs_info.info_files) eq 2 then file_struct_arr2 = fhd_file_setup(obs_info.info_files[1], $
    pol_inc, spec_window_type = spec_window_types[max_sw], freq_ch_range = freq_ch_range) $
  else file_struct_arr2 = file_struct_arr1
  
  if keyword_set(diff_ratio) or keyword_set(all_pol_diff_ratio) then begin
    type_pol_str = strarr(n_sets, 2)
    if keyword_set(all_pol_diff_ratio) then begin
      type_pol_str[0, *] = cube_types + '_' + pols[0]
      type_pol_str[1, *] = cube_types + '_' + pols[0]
      type_pol_str[2, *] = cube_types + '_' + pols[max_pol]
      type_pol_str[3, *] = cube_types + '_' + pols[max_pol]
    endif else begin
      type_pol_str[0, *] = cube_types + '_' + pols[0]
      type_pol_str[1, *] = cube_types + '_' + pols[max_pol]
    endelse
    
    type_pol_locs = intarr(n_sets, 2)
    for i=0, n_sets-1 do begin
      wh_type_pol1 = where(file_struct_arr1.type_pol_str eq type_pol_str[i,0], count_type_pol)
      if count_type_pol eq 0 then $
        message, 'requested type_pol not found: ' + type_pol_str[i,0] + ' not in ' + file_struct_arr1.power_savefile else type_pol_locs[i,0] = wh_type_pol1
      wh_type_pol2 = where(file_struct_arr2.type_pol_str eq type_pol_str[i,1], count_type_pol)
      if count_type_pol eq 0 then $
        message, 'requested type_pol not found: ' + type_pol_str[i,1] + ' not in ' + file_struct_arr2.power_savefile else type_pol_locs[i,1] = wh_type_pol2
    endfor
    
    if keyword_set(all_pol_diff_ratio) then begin
      titles = strarr(2,3)
      for i=0, 1 do titles[i,*] = [obs_info.fhd_types[0] + ' ' + type_pol_str[2*i,0] + '/' + type_pol_str[2*i,1], $
        obs_info.fhd_types[1] + ' ' + type_pol_str[2*i+1,0] + '/' + type_pol_str[2*i+1,1], 'Ratio Difference']
    endif else titles = [obs_info.fhd_types[0] + ' ' + type_pol_str[0,0] + '/' + type_pol_str[0,1], $
      obs_info.fhd_types[1] + ' ' + type_pol_str[1,0] + '/' + type_pol_str[1,1], 'Ratio Difference']
      
  endif else begin
    type_pol_str = [cube_types[0] + '_' + pols[0], cube_types[max_type] + '_' + pols[max_pol]]
    
    type_pol_locs = intarr(2)
    wh_type_pol1 = where(file_struct_arr1.type_pol_str eq type_pol_str[0], count_type_pol)
    if count_type_pol eq 0 then $
      message, 'requested type_pol not found: ' + type_pol_str[0] + ' not in ' + file_struct_arr1.power_savefile else type_pol_locs[0] = wh_type_pol1
    wh_type_pol2 = where(file_struct_arr2.type_pol_str eq type_pol_str[1], count_type_pol)
    if count_type_pol eq 0 then $
      message, 'requested type_pol not found: ' + type_pol_str[1] + ' not in ' + file_struct_arr2.power_savefile else type_pol_locs[1] = wh_type_pol2
      
    title = type_pol_str[0] + '/' + type_pol_str[1]
    
  endelse
  
  if tag_exist(file_struct_arr1, 'n_obs') then n_obs1 = file_struct_arr1[0].n_obs
  if tag_exist(file_struct_arr2, 'n_obs') then n_obs2 = file_struct_arr2[0].n_obs
  
  if n_elements(n_obs1) gt 0 and n_elements(n_obs2) gt 0 then begin
    if n_elements(note) eq 0 then note = '(' + number_formatter(round(mean(file_struct_arr1[0].n_obs))) + ',' + number_formatter(round(mean(file_struct_arr2[0].n_obs))) + ')' $
    else note = note + ' (' + number_formatter(round(mean(file_struct_arr1[0].n_obs))) + ',' + number_formatter(round(mean(file_struct_arr2[0].n_obs))) + ')'
  endif
  
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
  
  fadd_2dbin = ''
  ;;if keyword_set(fill_holes) then fadd_2dbin = fadd_2dbin + '_nohole'
  if keyword_set(no_kzero) then fadd_2dbin = fadd_2dbin + '_nok0'
  if keyword_set(log_kpar) then fadd_2dbin = fadd_2dbin + '_logkpar'
  if keyword_set(log_kperp) then fadd_2dbin = fadd_2dbin + '_logkperp'
  
  if keyword_set(diff_ratio) or keyword_set(all_pol_diff_ratio) then begin
    savefiles_2d = strarr(n_sets, 2)
    for j=0, n_sets/2-1 do begin
      for i=0, 1 do begin
        savefiles_2d[2*j, i] = file_struct_arr1[type_pol_locs[2*j, i]].savefile_froot + file_struct_arr1[type_pol_locs[2*j, i]].savefilebase + file_struct_arr1[0].power_tag
        savefiles_2d[2*j+1, i] =file_struct_arr2[type_pol_locs[2*j+1, i]].savefile_froot + file_struct_arr2[type_pol_locs[2*j+1, i]].savefilebase + file_struct_arr2[0].power_tag
      endfor
    endfor
    savefiles_2d = savefiles_2d + fadd_2dbin + '_2dkpower.idlsave'
    
  endif else savefiles_2d = [file_struct_arr1[type_pol_locs[0]].savefile_froot + file_struct_arr1[type_pol_locs[0]].savefilebase + file_struct_arr1[0].power_tag, $
    file_struct_arr2[type_pol_locs[1]].savefile_froot + file_struct_arr2[type_pol_locs[1]].savefilebase + file_struct_arr2[0].power_tag] + $
    fadd_2dbin + '_2dkpower.idlsave'
    
  test_save_2d = file_test(savefiles_2d) *  (1 - file_test(savefiles_2d, /zero_length))
  
  if min(test_save_2d) eq 0 then message, '2D savefile not found: ' + savefiles_2d[where(test_save_2d eq 0)]
  
  if pub then begin
    plotfile_2d = plot_path + plot_filebase + '_2dkpower' + plot_exten
  endif else plot_exten = ''
  
  if keyword_set(diff_ratio) or keyword_set(all_pol_diff_ratio) then begin
  
    for i=0, n_sets/2-1 do begin
      kperp_edges = getvar_savefile(savefiles_2d[2*i,0], 'kperp_edges')
      if total(abs(kperp_edges - getvar_savefile(savefiles_2d[2*i+1,0], 'kperp_edges'))) ne 0 then message, 'kperp_edges do not match in savefiles'
      kpar_edges = getvar_savefile(savefiles_2d[2*i,0], 'kpar_edges')
      if total(abs(kpar_edges - getvar_savefile(savefiles_2d[2*i+1,0], 'kpar_edges'))) ne 0 then message, 'kpar_edges do not match in savefiles'
      kperp_bin = getvar_savefile(savefiles_2d[2*i,0], 'kperp_bin')
      if total(abs(kperp_bin - getvar_savefile(savefiles_2d[2*i+1,0], 'kperp_bin'))) ne 0 then message, 'kperp_bin do not match in savefiles'
      kpar_bin = getvar_savefile(savefiles_2d[2*i,0], 'kpar_bin')
      if total(abs(kpar_bin - getvar_savefile(savefiles_2d[2*i+1,0], 'kpar_bin'))) ne 0 then message, 'kpar_bin do not match in savefiles'
      kperp_lambda_conv = getvar_savefile(savefiles_2d[2*i,0], 'kperp_lambda_conv')
      if total(abs(kperp_lambda_conv - getvar_savefile(savefiles_2d[2*i+1,0], 'kperp_lambda_conv'))) ne 0 then message, 'kperp_lambda_conv do not match in savefiles'
      delay_params = getvar_savefile(savefiles_2d[2*i,0], 'delay_params')
      if total(abs(delay_params - getvar_savefile(savefiles_2d[2*i+1,0], 'delay_params'))) ne 0 then message, 'delay_params do not match in savefiles'
      hubble_param = getvar_savefile(savefiles_2d[2*i,0], 'hubble_param')
      if total(abs(hubble_param - getvar_savefile(savefiles_2d[2*i+1,0], 'hubble_param'))) ne 0 then message, 'hubble_param do not match in savefiles'
      
      power1 = getvar_savefile(savefiles_2d[2*i,0], 'power')
      power2 = getvar_savefile(savefiles_2d[2*i,1], 'power')
      
      power_ratio1 = power1 / power2
      wh0 = where(power2 eq 0, count0)
      if count0 gt 0 then power_ratio1[wh0] = 0
      
      power1 = getvar_savefile(savefiles_2d[2*i+1,0], 'power')
      power2 = getvar_savefile(savefiles_2d[2*i+1,1], 'power')
      
      power_ratio2 = power1 / power2
      wh0 = where(power2 eq 0, count0)
      if count0 gt 0 then power_ratio2[wh0] = 0
      
      power_diff_ratio = power_ratio1-power_ratio2
      
      if i eq 0 then begin
        kperp_plot_range = [5./kperp_lambda_conv, min([file_struct_arr1.kspan/2.,file_struct_arr1.max_baseline_lambda])/kperp_lambda_conv]
        if keyword_set(hinv) then kperp_plot_range = kperp_plot_range / hubble_param
        
        
        ncol=3
        nrow=n_sets/2
        start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
      endif else pos_use = positions[*,3*i]
      
      kpower_2d_plots, power = power_ratio1, multi_pos = pos_use, start_multi_params = start_multi_params, $
        kperp_edges = kperp_edges, kpar_edges = kpar_edges, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
        kperp_lambda_conv = kperp_lambda_conv, delay_params = delay_params, hubble_param = hubble_param, $
        png = png, eps = eps, pdf = pdf, plotfile = plotfile_2d, window_num = window_num, $
        kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
        data_range = data_range, full_title = titles[i,0], $
        plot_wedge_line = plot_wedge_line, hinv = hinv, /pwr_ratio, $
        wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, $
        kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
        
      if i eq 0 then begin
        positions = pos_use
        undefine, start_multi_params
      endif
      pos_use = positions[*,3*i+1]
      
      kpower_2d_plots, power = power_ratio2, multi_pos = pos_use, start_multi_params = start_multi_params, $
        kperp_edges = kperp_edges, kpar_edges = kpar_edges, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
        kperp_lambda_conv = kperp_lambda_conv, delay_params = delay_params, hubble_param = hubble_param, $
        png = png, eps = eps, pdf = pdf, plotfile = plotfile_2d, window_num = window_num, $
        kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
        data_range = data_range, full_title = titles[i,1], $
        plot_wedge_line = plot_wedge_line, hinv = hinv, /pwr_ratio, $
        wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, $
        kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
        
      pos_use = positions[*,3*i+2]
      
      kpower_2d_plots, power = power_diff_ratio, multi_pos = pos_use, start_multi_params = start_multi_params, $
        kperp_edges = kperp_edges, kpar_edges = kpar_edges, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
        kperp_lambda_conv = kperp_lambda_conv, delay_params = delay_params, hubble_param = hubble_param, $
        png = png, eps = eps, pdf = pdf, plotfile = plotfile_2d, window_num = window_num, color_profile = 'sym_log', $
        kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
        data_range = diff_range, data_min_abs = diff_min_abs, full_title = titles[i,2], note = note, $
        plot_wedge_line = plot_wedge_line, hinv = hinv, /pwr_ratio, $
        wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, $
        kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
        
    endfor
    
    if keyword_set(pub) then begin
      cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density=600
    endif
    
    
  endif else begin
  
    kperp_edges = getvar_savefile(savefiles_2d[0], 'kperp_edges')
    if total(abs(kperp_edges - getvar_savefile(savefiles_2d[1], 'kperp_edges'))) ne 0 then message, 'kperp_edges do not match in savefiles'
    kpar_edges = getvar_savefile(savefiles_2d[0], 'kpar_edges')
    if total(abs(kpar_edges - getvar_savefile(savefiles_2d[1], 'kpar_edges'))) ne 0 then message, 'kpar_edges do not match in savefiles'
    kperp_bin = getvar_savefile(savefiles_2d[0], 'kperp_bin')
    if total(abs(kperp_bin - getvar_savefile(savefiles_2d[1], 'kperp_bin'))) ne 0 then message, 'kperp_bin do not match in savefiles'
    kpar_bin = getvar_savefile(savefiles_2d[0], 'kpar_bin')
    if total(abs(kpar_bin - getvar_savefile(savefiles_2d[1], 'kpar_bin'))) ne 0 then message, 'kpar_bin do not match in savefiles'
    kperp_lambda_conv = getvar_savefile(savefiles_2d[0], 'kperp_lambda_conv')
    if total(abs(kperp_lambda_conv - getvar_savefile(savefiles_2d[1], 'kperp_lambda_conv'))) ne 0 then message, 'kperp_lambda_conv do not match in savefiles'
    delay_params = getvar_savefile(savefiles_2d[0], 'delay_params')
    if total(abs(delay_params - getvar_savefile(savefiles_2d[1], 'delay_params'))) ne 0 then message, 'delay_params do not match in savefiles'
    hubble_param = getvar_savefile(savefiles_2d[0], 'hubble_param')
    if total(abs(hubble_param - getvar_savefile(savefiles_2d[1], 'hubble_param'))) ne 0 then message, 'hubble_param do not match in savefiles'
    
    power1 = getvar_savefile(savefiles_2d[0], 'power')
    power2 = getvar_savefile(savefiles_2d[1], 'power')
    
    power_ratio = power1 / power2
    wh0 = where(power2 eq 0, count0)
    if count0 gt 0 then power_ratio[wh0] = 0
    
    kperp_plot_range = [5./kperp_lambda_conv, min([file_struct_arr1.kspan/2.,file_struct_arr1.max_baseline_lambda])/kperp_lambda_conv]
    if keyword_set(hinv) then kperp_plot_range = kperp_plot_range / hubble_param
    
    
    kpower_2d_plots, power = power_ratio, kperp_edges = kperp_edges, kpar_edges = kpar_edges, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
      kperp_lambda_conv = kperp_lambda_conv, delay_params = delay_params, hubble_param = hubble_param, $
      png = png, eps = eps, pdf = pdf, plotfile = plotfile_2d, window_num = window_num, $
      kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
      data_range = data_range, title_prefix = title, note = note, $
      plot_wedge_line = plot_wedge_line, hinv = hinv, /pwr_ratio, $
      wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis, $
      kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
      
  endelse
  
end