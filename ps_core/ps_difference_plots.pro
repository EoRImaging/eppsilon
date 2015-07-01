pro ps_difference_plots, folder_names, obs_info, cube_types, pols, all_type_pol = all_type_pol, $
    refresh_diff = refresh_diff, freq_ch_range = freq_ch_range, $
    plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    note = note, spec_window_types = spec_window_types, data_range = data_range, data_min_abs = data_min_abs, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, plot_1d = plot_1d, axis_type_1d=axis_type_1d, $
    wt_cutoffs = wt_cutoffs, wt_measures = wt_measures, $
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
  
  
  if n_elements(freq_ch_range) ne 0 then begin
    if min(freq_ch_range) lt 0 or max(freq_ch_range) - min(freq_ch_range) lt 3 then message, 'invalid freq_ch_range'
    fch_tag = '_ch' + number_formatter(min(freq_ch_range)) + '-' + number_formatter(max(freq_ch_range))
  endif else fch_tag = ''
  
  
  n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), n_elements(pols), n_elements(spec_window_types)])
  if n_diffs gt 2 then message, 'only 1 or 2 of [folder_names, obs_names, cube_types, pols, spec_window_types] allowed'
  
  if n_elements(obs_info.info_files) eq 2 or n_elements(spec_window_types) eq 2 $
    and n_elements(cube_types) eq 0 and n_elements(pols) eq 0 and n_elements(all_type_pol) eq 0 then all_type_pol = 1
    
  if keyword_set(all_type_pol) and n_elements(info_files) eq 1 and n_elements(spec_window_types) lt 2 then $
    message, 'all_type_pol can only be set with 2 folder names and/or 2 obs names and/or 2 spec windows'
    
  if not keyword_set(all_type_pol) then begin
    if n_elements(cube_types) eq 0 then if n_diffs eq 1 then cube_types = ['dirty', 'res'] else cube_types = 'res'
    n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), n_elements(pols), n_elements(spec_window_types)])
    if n_elements(pols) eq 0 then if n_diffs eq 1 then pols=['xx', 'yy'] else pols = 'xx'
    n_diffs = max([n_elements(obs_info.info_files), n_elements(cube_types), n_elements(pols), n_elements(spec_window_types)])
    
    if n_diffs eq 1 then message, 'at least one of info_files, cube_types, pols, spec_window_types must be a 2 element vector'
    
    n_cubes = 1
  endif else begin
    undefine, cube_types, pols
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
  
  if n_elements(axis_type_1d) eq 0 then axis_type_1d = 'sym_log'
  
  if n_elements(folder_names) eq 2 then begin
    if n_elements(save_path) eq 0 then save_path = obs_info.diff_save_path
    note = obs_info.diff_note
    if tag_exist(obs_info, 'diff_plot_path') then plot_path = obs_info.diff_plot_path else plot_path = save_path
  endif else begin
    if n_elements(save_path) eq 0 then save_path = obs_info.save_paths[0]
    note = obs_info.fhd_types[0]
    plot_path = obs_info.plot_paths[0]
  endelse
  
  if n_elements(spec_window_types) eq 2 then note = note + ' ' + spec_window_types[0] + ' minus ' + spec_window_types[1]
  
  
  ;; density correction defaults & file naming for 2D & 1D files
  if n_elements(wt_cutoffs) eq 0 then begin
    ;; default to wt_cutoffs = 1, wt_measures = 'min'
    wt_cutoffs = 1
    wt_measures = 'min'
  endif else if n_elements(wt_measures) eq 0 then begin
    print, 'wt_cutoffs is specified but wt_measures is not. Defaulting wt_measures to "min".'
    wt_measures = strarr(n_elements(wt_cutoffs)) + 'min'
  endif
  
  kperp_density_names = strarr(n_elements(wt_cutoffs))
  wh_cutoff0 = where(wt_cutoffs eq 0, count_cutoff0, complement = wh_cutoff_n0, ncomplement = count_cutoff_n0)
  wh_std = where(wt_cutoffs eq 0 and wt_measures eq 'min', count_std)
  
  if count_cutoff0 gt 0 then kperp_density_names[wh_cutoff0] = '_nodensitycorr'
  if count_cutoff_n0 gt 0 then kperp_density_names[wh_cutoff_n0] = '_kperp_density_' + wt_measures[wh_cutoff_n0] + '_gt' + number_formatter(wt_cutoffs[wh_cutoff_n0])
  
  if count_std gt 1 then kperp_density_names[wh_std] = '_dencorr'
  
  
  
  if keyword_set(all_type_pol) then begin
  
    if n_elements(folder_names) eq 2 and folder_names[0] ne folder_names[1] then begin
      plot_filebase = obs_info.name_same_parts + '__' + obs_info.name_diff_parts[0] + '_' + obs_info.obs_names[0] + sw_tags[0] + $
        '_minus_' + obs_info.name_diff_parts[1]  + '_' + obs_info.obs_names[1] + sw_tags[max_sw]
    endif else plot_filebase = obs_info.folder_basenames[0] + '__' + obs_info.obs_names[0] + sw_tags[0] + '_minus_' + $
      obs_info.obs_names[1] + sw_tags[max_sw]
      
  endif else begin
  
    if n_elements(folder_names) eq 1 then begin
      if n_elements(obs_info.obs_names) gt 1 then begin
        plot_filebase = obs_info.folder_basenames[0] + '_' + obs_info.obs_names[0] + '_' + cube_types[0] + '_' + pols[0] + sw_tags[0] + $
          '_minus_' + obs_info.obs_names[0] + '_' + cube_types[max_type] + '_' + pols[max_pol] + sw_tags[max_sw]
      endif else begin
        if obs_info.integrated[0] eq 0 then plot_start = obs_info.folder_basenames[0] + '_' + obs_info.obs_names[0] else plot_start = obs_info.fhd_types[0]
        
        plot_filebase = plot_start + '_' + cube_types[0] + '_' + pols[0] + sw_tags[0] + $
          '_minus_' + cube_types[max_type] + '_' + pols[max_pol] + sw_tags[max_sw]
      endelse
    endif else plot_filebase = obs_info.name_same_parts + '__' + strjoin([obs_info.name_diff_parts[0], cube_types[0], pols[0]], '_') + sw_tags[0] + $
      '_minus_' + strjoin([obs_info.name_diff_parts[1], cube_types[max_type], pols[max_pol]], '_') + sw_tags[max_sw]
  endelse
  plot_filebase = plot_filebase + fch_tag
  
  if not file_test(plot_path, /directory) then file_mkdir, plot_path
  if not file_test(save_path, /directory) then file_mkdir, save_path
  
  file_struct_arr1 = fhd_file_setup(obs_info.info_files[0], $
    spec_window_type = spec_window_types[0], freq_ch_range = freq_ch_range)
  if n_elements(obs_info.info_files) eq 2 then file_struct_arr2 = fhd_file_setup(obs_info.info_files[1], $
    spec_window_type = spec_window_types[max_sw], freq_ch_range = freq_ch_range) $
  else file_struct_arr2 = file_struct_arr1
  type_pol_str1 = file_struct_arr1.type_pol_str
  type_pol_str2 = file_struct_arr2.type_pol_str
  
  if keyword_set(all_type_pol) then begin
    type_pol_str1 = file_struct_arr1.type_pol_str
    type_pol_str2 = file_struct_arr2.type_pol_str
    
    if n_elements(type_pol_str1) ne n_elements(type_pol_str2) then message, 'all_type_pol cannot be used with these folders, they contain different number of types & pols'
    n_cubes = n_elements(type_pol_str1)
    
    for i=0, n_cubes-1 do begin
      temp = where(type_pol_str2 eq type_pol_str1[i], count_typepol)
      if count_typepol eq 0 then message, 'all_type_pol cannot be used with these folders, they contain different sets of types & pols'
    endfor
    
    type_pol_str = type_pol_str1
  endif else begin
    type_pol_str = [cube_types[0] + '_' + pols[0], cube_types[max_type] + '_' + pols[max_pol]]
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
  
  savefiles_1d = strarr(n_cubes, n_elements(wedge_1dbin_names))
  titles = strarr(n_cubes)
  for i=0, n_cubes-1 do begin
  
    if n_cubes gt 1 then begin
      type_pol1 = type_pol_str[i]
      type_pol2 = type_pol_str[i]
      titles[i] = type_pol_str[i]
    endif else begin
      type_pol1 = type_pol_str[0]
      type_pol2 = type_pol_str[1]
      titles[i] = type_pol_str[0] + '-' + type_pol_str[1]
    endelse
    
    ;; if not save_path specified, save in folder of first info file
    if n_elements(save_path) eq 0 then save_path = file_dirname(obs_info.info_files[0], /mark_directory)
    if n_elements(plot_path) eq 0 then plot_path = save_path
    
    if n_elements(obs_info.info_files) eq 2 then begin
      fileparts_1 = strsplit(file_struct_arr1[0].general_filebase, '_', /extract)
      fileparts_2 = strsplit(file_struct_arr2[0].general_filebase, '_', /extract)
      match_test = strcmp(fileparts_1, fileparts_2)
      wh_diff = where(match_test eq 0, count_diff, complement = wh_same, ncomplement = count_same)
    endif
    
    if n_elements(savefilebase) eq 0 then begin
      if n_elements(obs_info.info_files) eq 1 then begin
        savefilebase_use = file_struct_arr1[0].general_filebase + '_' + type_pol1 + '_minus_' + type_pol2
      endif else begin
        if count_diff eq 0 then begin
          savefilebase_use = file_struct_arr1[0].general_filebase + '_' + type_pol1
          if type_pol1 ne type_pol2 then savefilebase_use = savefilebase_use + '_' + type_pol1 + '_minus_' + type_pol2
        endif else begin
          if count_same gt 0 then begin
            if type_pol1 ne type_pol2 then savefilebase_use = strjoin(fileparts_1[wh_same], '_') + '__' + $
              strjoin(fileparts_1[wh_diff]) + '_' + type_pol1 + '_minus_' + strjoin(fileparts_2[wh_diff]) + '_' + type_pol2 $
            else savefilebase_use = strjoin(fileparts_1[wh_same], '_') + '__' + strjoin(fileparts_1[wh_diff]) + '_minus_' + strjoin(fileparts_2[wh_diff]) + '__' + type_pol1
          endif else begin
            if type_pol1 ne type_pol2 then savefilebase_use = file_struct_arr1[0].general_filebase + '_' + type_pol_str[0] + $
              '_minus_' + file_struct_arr2[0].general_filebase + '_' + type_pol_str[1] $
            else savefilebase_use = file_struct_arr1[0].general_filebase + '_minus_' + file_struct_arr2[0].general_filebase + '__' + type_pol1
          endelse
        endelse
      endelse
    endif else begin
      if type_pol1 eq type_pol2 then savefilebase_use = savefilebase + '_' + type_pol1 else savefilebase_use = savefilebase + '_' + type_pol1 + '_minus_' + type_pol2
    endelse
    if n_elements(plot_filebase) eq 0 then plot_filebase = savefilebase_use
    
    savefile = save_path + savefilebase_use + '_power.idlsave'
    savefile_2d = save_path + savefilebase_use + kperp_density_names + '_2dkpower.idlsave'
    
    for wedge_i=0, n_elements(wedge_1dbin_names)-1 do begin
      savefiles_1d[i, wedge_i] = save_path + savefilebase_use + kperp_density_names + wedge_1dbin_names[wedge_i] + '_1dkpower.idlsave'
    endfor
    
    test_save = file_test(savefile) *  (1 - file_test(savefile, /zero_length))
    test_save_2d = file_test(savefile_2d) *  (1 - file_test(savefile_2d, /zero_length))
    test_save_1d = file_test(reform(savefiles_1d[i,*])) *  (1 - file_test(reform(savefiles_1d[i,*]) , /zero_length))
    
    
    if test_save eq 0 or keyword_set(refresh_diff) then begin
      cube_ind1 = where(type_pol_str1 eq type_pol1, count_typepol)
      if count_typepol eq 0 then message, 'type/pol ' + type_pol1 + ' not included in info_file: ' + obs_info.info_files[0]
      
      power_file1 = file_struct_arr1[cube_ind1].power_savefile
      if file_test(power_file1) eq 0 then message, 'No power file for ' + type_pol1 + ' and info_file: ' + obs_info.info_files[0]
      
      cube_ind2 = where(file_struct_arr2.type_pol_str eq type_pol2, count_typepol)
      if count_typepol eq 0 then message, 'type/pol ' + type_pol2 + ' not included in info_file: ' + obs_info.info_files[n_elements(obs_info.info_files)-1]
      
      power_file2 = file_struct_arr2[cube_ind2].power_savefile
      if file_test(power_file2) eq 0 then message, 'No power file for ' + type_pol2 + ' and info_file: ' + obs_info.info_files[n_elements(obs_info.info_files)-1]
      
      kx_mpc = getvar_savefile(power_file1, 'kx_mpc')
      kx_mpc2 = getvar_savefile(power_file2, 'kx_mpc')
      if n_elements(kx_mpc) ne n_elements(kx_mpc2) or max(abs(kx_mpc - kx_mpc2)) gt 1.05e-3 then message, 'kx_mpc does not match between cubes'
      ky_mpc = getvar_savefile(power_file1, 'ky_mpc')
      ky_mpc2 = getvar_savefile(power_file2, 'ky_mpc')
      if n_elements(ky_mpc) ne n_elements(ky_mpc2) or max(abs(ky_mpc - ky_mpc2)) gt 1.05e-3 then message, 'ky_mpc does not match between cubes'
      kz_mpc = getvar_savefile(power_file1, 'kz_mpc')
      kz_mpc2 = getvar_savefile(power_file2, 'kz_mpc')
      if n_elements(kz_mpc) ne n_elements(kz_mpc2) or max(abs(kz_mpc - kz_mpc)) gt 1.05e-3 then message, 'kz_mpc does not match between cubes'
      
      ;; if kx, ky, and kz are all the same these should be too
      kperp_lambda_conv = getvar_savefile(power_file1, 'kperp_lambda_conv')
      delay_params = getvar_savefile(power_file1, 'delay_params')
      hubble_param = getvar_savefile(power_file1, 'hubble_param')
      
      power1 = getvar_savefile(power_file1, 'power_3d')
      power2 = getvar_savefile(power_file2, 'power_3d')
      power_diff = power1 - power2
      if max(abs(power_diff)) eq 0 then begin
        print, 'The cubes are identical -- power difference is zero everywhere'
      ;continue
      endif
      undefine, power1, power2
      
      weights1 = getvar_savefile(power_file1, 'weights_3d')
      weights2 = getvar_savefile(power_file2, 'weights_3d')
      
      ;; variance_3d = 1/weights_3d
      var1 = 1./weights1
      wh_wt1_0 = where(weights1 eq 0, count_wt1_0)
      if count_wt1_0 gt 0 then var1[wh_wt1_0] = 0
      var2 = 1./weights2
      wh_wt2_0 = where(weights2 eq 0, count_wt2_0)
      if count_wt2_0 gt 0 then var2[wh_wt2_0] = 0
      undefine, weights1, weights2
      
      var_diff = var1 + var2
      weight_diff = 1/var_diff
      if count_wt1_0 gt 0 then weight_diff[wh_wt1_0] = 0
      if count_wt2_0 gt 0 then weight_diff[wh_wt2_0] = 0
      undefine, var1, var2, var_diff
      
      
      wt_meas_ave1 = getvar_savefile(power_file1, 'wt_meas_ave')
      wt_meas_ave2 = getvar_savefile(power_file2, 'wt_meas_ave')
      wt_meas_ave = wt_meas_ave1 < wt_meas_ave2
      
      wt_meas_min1 = getvar_savefile(power_file1, 'wt_meas_min')
      wt_meas_min2 = getvar_savefile(power_file2, 'wt_meas_min')
      wt_meas_min = wt_meas_min1 < wt_meas_min2
      
      save, file = savefile, power_diff, weight_diff, wt_meas_ave, wt_meas_min, $
        kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, delay_params, hubble_param
        
    endif else if min([test_save_1d, test_save_2d]) eq 0 then restore, savefile
    
    if test_save_2d eq 0 or keyword_set(refresh_diff) then begin
    
      if wt_cutoffs gt 0 then begin
        case wt_measures of
          'ave': wt_meas_use = wt_meas_ave
          'min': wt_meas_use = wt_meas_min
        endcase
      endif
      
      power_rebin = kspace_rebinning_2d(power_diff, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, log_kpar = log_kpar, $
        log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, weights = weight_diff, $
        binned_weights = binned_weights, $
        kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoffs)
        
      power = power_rebin
      kperp_edges = kperp_edges_mpc
      kpar_edges = kpar_edges_mpc
      weights = binned_weights
      
      save, file = savefile_2d, power, weights, kperp_edges, kpar_edges, kperp_bin, kpar_bin, $
        kperp_lambda_conv, delay_params, hubble_param
    endif else restore, savefile_2d
    
    if min(test_save_1d) eq 0 or keyword_set(refresh_diff) then begin
    
      for wedge_i=0, n_elements(wedge_amp) do begin
        if wedge_i gt 0 then wedge_amp_use = wedge_amp[wedge_i-1]
        
        if wt_cutoffs gt 0 then begin
          case wt_measures of
            'ave': wt_meas_use = wt_meas_ave
            'min': wt_meas_use = wt_meas_min
          endcase
        endif
        
        power_rebin = kspace_rebinning_1d(power_diff, kx_mpc, ky_mpc, kz_mpc, k_edges_mpc, k_bin = k_bin, log_k = log_k, $
          weights = weight_diff, binned_weights = binned_weights, wedge_amp = wedge_amp_use, $
          kperp_density_measure = wt_meas_use, kperp_density_cutoff = wt_cutoffs)
          
        power = power_rebin
        k_edges = k_edges_mpc
        weights = binned_weights
        
        save, file = savefiles_1d[i, wedge_i], power, weights, k_edges, k_bin, kperp_lambda_conv, delay_params, hubble_param
      endfor
    endif
    
    
    if pub then begin
      plotfile_2d = plot_path + plot_filebase + kperp_density_names + '_2dkpower' + plot_exten
      if i eq 0 then plotfiles_1d = plot_path + plot_filebase + ['', wedge_1dbin_names] + kperp_density_names + '_1dkpower' + plot_exten
    endif else plot_exten = ''
    
    
    kperp_plot_range = [5./kperp_lambda_conv, min([file_struct_arr1.kspan/2.,file_struct_arr1.max_baseline_lambda])/kperp_lambda_conv]
    if keyword_set(hinv) then kperp_plot_range = kperp_plot_range / hubble_param
    
    
    if keyword_set(pub) then font = 1 else font = -1
    
    if n_cubes gt 1 then begin
      if i eq 0 then begin
        if n_cubes eq 6 then begin
          if keyword_set(kperp_linear_axis) then begin
            ;; aspect ratio doesn't work out for kperp_linear with multiple rows
            ncol = 6
            nrow = 1
          endif else begin
            ncol = 3
            nrow = 2
          endelse
        endif else begin
          if keyword_set(kperp_linear_axis) then begin
            ;; aspect ratio doesn't work out for kperp_linear with multiple rows
            ncol = n_cubes
            nrow = 1
          endif else begin
            nrow = 2
            ncol = ceil(n_cubes/nrow)
          endelse
        endelse
        start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
        
        if n_elements(window_num) eq 0 then window_num = 1
      endif else begin
        pos_use = positions[*,i]
        
      endelse
    endif
    
    if i eq n_cubes-1 and n_elements(note) gt 0 then note_use = note else undefine, note_use
    
    kpower_2d_plots, savefile_2d, multi_pos = pos_use, start_multi_params = start_multi_params, $
      kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, note = note_use, $
      data_range = data_range, data_min_abs = data_min_abs, png = png, eps = eps, pdf = pdf, plotfile = plotfile_2d, full_title=titles[i], $
      window_num = window_num, color_profile = 'sym_log', $
      kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, baseline_axis = baseline_axis, delay_axis = delay_axis, $
      wedge_amp = wedge_amp, plot_wedge_line = plot_wedge_line, hinv = hinv
      
    if n_cubes gt 1 and i eq 0 then begin
      positions = pos_use
      undefine, start_multi_params
    endif
    
  endfor
  undefine, positions, pos_use
  
  if keyword_set(pub) and n_cubes gt 1 then begin
    cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density=600
  endif
  
  if keyword_set(plot_1d) then begin
    for i=0, n_elements(wedge_amp) do begin
    
      if keyword_set(pub) then plotfiles_use = plotfiles_1d[i]
      
      for j=0, n_cubes-1 do begin
      
        if n_cubes gt 1 then begin
          if j eq 0 then begin
            if n_cubes eq 6 then begin
              ncol = 3
              nrow = 2
            endif else begin
              nrow = 2
              ncol = ceil(n_cubes/nrow)
            endelse
            start_multi_params = {ncol:ncol, nrow:nrow, ordering:'row'}
            undefine, positions, pos_use
            
            window_num = 2+i
          endif else begin
            pos_use = positions[*,j]
            
          endelse
        endif
        
        if i eq n_cubes-1 and n_elements(note) gt 0 then note_use = note else undefine, note_use
        
        kpower_1d_plots, savefiles_1d[j,i], window_num=window_num, start_multi_params = start_multi_params, multi_pos = pos_use, $
          names=titles[j], hinv = hinv, png = png, eps = eps, pdf = pdf, plotfile = plotfiles_use, note = note_use, yaxis_type = axis_type_1d
          
        if j eq 0 then begin
          positions = pos_use
          undefine, start_multi_params
        endif
        
      endfor
      
    endfor
  endif
  
end