pro ps_wrapper, folder_name, obs_name, data_subdirs=data_subdirs, exact_obsnames = exact_obsnames, ps_foldername = ps_foldername, $
    no_wtvar_rts = no_wtvar_rts, set_data_ranges = set_data_ranges, $
    beamfiles = beamfiles, rts = rts, casa = casa, sim = sim, $
    refresh_dft = refresh_dft, refresh_ps = refresh_ps, refresh_binning = refresh_binning, refresh_info = refresh_info, $
    refresh_beam = refresh_beam, dft_fchunk = dft_fchunk, dft_ian = dft_ian, $
    delta_uv_lambda = delta_uv_lambda, max_uv_lambda = max_uv_lambda, cut_image = cut_image, $
    pol_inc = pol_inc, type_inc = type_inc, $
    freq_ch_range = freq_ch_range, freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
    uvf_input = uvf_input, uv_avg = uv_avg, uv_img_clip = uv_img_clip, $
    std_power = std_power, inverse_covar_weight = inverse_covar_weight, no_wtd_avg = no_wtd_avg, norm_rts_with_fhd = norm_rts_with_fhd, $
    wt_cutoffs = wt_cutoffs, wt_measures = wt_measures, fix_sim_input = fix_sim_input, $
    no_spec_window = no_spec_window, spec_window_type = spec_window_type, $
    no_kzero = no_kzero, plot_slices = plot_slices, slice_type = slice_type, uvf_plot_type = uvf_plot_type, plot_stdset = plot_stdset, plot_1to2d = plot_1to2d, $
    plot_kpar_power = plot_kpar_power, plot_kperp_power = plot_kperp_power, plot_k0_power = plot_k0_power, plot_noise_1d = plot_noise_1d, $
    data_range = data_range, sigma_range = sigma_range, nev_range = nev_range, slice_range = slice_range, plot_sim_noise = plot_sim_noise, $
    snr_range = snr_range, noise_range = noise_range, nnr_range = nnr_range, range_1d = range_1d, $
    log_kpar = log_kpar, log_kperp = log_kperp, kpar_bin = kpar_bin, kperp_bin = kperp_bin, $
    log_k1d = log_k1d, k1d_bin = k1d_bin, plot_1d_delta = plot_1d_delta, plot_1d_error_bars = plot_1d_error_bars, plot_1d_nsigma = plot_1d_nsigma, $
    kperp_range_1dave = kperp_range_1dave, kperp_range_lambda_1dave = kperp_range_lambda_1dave, kpar_range_1dave = kpar_range_1dave, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, kperp_plot_range = kperp_plot_range, $
    kperp_lambda_plot_range = kperp_lambda_plot_range, kpar_plot_range = kpar_plot_range, $
    baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, hinv = hinv, $
    plot_wedge_line = plot_wedge_line, wedge_angles = wedge_angles, coarse_harm_width = coarse_harm_width, $
    plot_eor_1d = plot_eor_1d, plot_flat_1d = plot_flat_1d, no_text_1d = no_text_1d, $
    save_path = save_path, savefilebase = savefilebase, plot_path = plot_path, plot_filebase = plot_filebase, $
    individual_plots = individual_plots, plot_binning_hist = plot_binning_hist, $
    note = note, png = png, eps = eps, pdf = pdf, cube_power_info = cube_power_info, no_dft_progress = no_dft_progress
    
  if n_elements(folder_name) ne 1 then message, 'one folder_name must be supplied.'
  
  spawn, 'hostname', hostname
  if stregex(hostname, 'mit.edu', /boolean) eq 1 then loc_name = 'mit'
  if stregex(hostname, 'enterprise', /boolean) eq 1 then loc_name = 'enterprise'
  case loc_name of
    'mit':  folder_name = mit_folder_locs(folder_name, rts = rts)
    'enterprise': folder_name = enterprise_folder_locs(folder_name, rts = rts)
    else: begin
      folder_test = file_test(folder_name, /directory)
      if folder_test eq 0 then message, 'machine not recognized and folder not found'
    endelse
  endcase
  
  if keyword_set(rts) and n_elements(dirty_folder) gt 0 then begin
    case loc_name of
      'mit':  dirty_folder = mit_folder_locs(dirty_folder, rts = rts)
      'enterprise': dirty_folder = enterprise_folder_locs(dirty_folder, rts = rts)
      else: begin
        folder_test = file_test(dirty_folder, /directory)
        if folder_test eq 0 then message, 'machine not recognized and dirty_folder not found'
      endelse
    endcase
  endif
  
  obs_info = ps_filenames(folder_name, obs_name, dirty_folder = dirty_folder, exact_obsnames = exact_obsnames, rts = rts, sim = sim, $
    uvf_input = uvf_input, data_subdirs = data_subdirs, ps_foldernames = ps_foldername, save_paths = save_path, plot_paths = plot_path, $
    refresh_info = refresh_info)
    
  if not file_test(save_path, /directory) then file_mkdir, save_path
  if not file_test(plot_path, /directory) then file_mkdir, plot_path
  
  if keyword_set(rts) then begin
  
    if tag_exist(obs_info, 'dirtyfiles') then dirtyfiles = obs_info.dirtyfiles.(0)
    
    if obs_info.info_files[0] ne '' then datafile = obs_info.info_files[0] else begin
    
      if obs_info.cube_files.(0)[0] ne '' then datafile = obs_info.cube_files.(0) else begin
        datafile = rts_fits2idlcube(obs_info.datafiles.(0), obs_info.weightfiles.(0), obs_info.variancefiles.(0), $
          dirtyfiles = dirtyfiles, pol_inc = pol_inc, save_path = obs_info.folder_names[0]+path_sep(), $
          refresh = refresh_dft, no_wtvar = no_wtvar_rts)
      endelse
    endelse
    
    if keyword_set(refresh_rtscube) then $
      datafile = rts_fits2idlcube(obs_info.datafiles.(0), obs_info.weightfiles.(0), obs_info.variancefiles.(0), $
      dirtyfiles = dirtyfiles, pol_inc = pol_inc, save_path = obs_info.folder_names[0]+path_sep(), $
      /refresh, no_wtvar = no_wtvar_rts)
      
    if keyword_set(no_wtvar_rts) then return
    
    max_uv_lambda = 300
    
    note = obs_info.rts_types
    
  endif else begin
  
    if obs_info.info_files[0] ne '' then datafile = obs_info.info_files[0] else datafile = obs_info.cube_files.(0)
    
    if keyword_set(uvf_input) then plot_filebase = obs_info.fhd_types[0] + '_uvf'$
    else  plot_filebase = obs_info.fhd_types[0]
    
    if stregex(obs_info.fhd_types[0], obs_info.obs_names[0], /boolean) eq 0 then $
      plot_filebase = plot_filebase + '_' + obs_info.obs_names[0]
      
    note = obs_info.fhd_types[0]
    if keyword_set(uvf_input) then note = note + '_uvf'
    
    if tag_exist(obs_info, 'beam_files') then beamfiles = obs_info.beam_files
    
    if keyword_set(sim) then begin
      plot_eor_1d=1
      if n_elements(range_1d) eq 0 then range_1d = [1e3, 1e7]
      
      if keyword_set(use_weight_cutoff) then begin
        wt_cutoffs = [0,1]
        wt_measures = strarr(2)+'min'
      endif
    endif
    
    if n_elements(set_data_ranges) eq 0 and not keyword_set(sim) then set_data_ranges = 1
    
  endelse
  
  print,'datafile = '+datafile
  
  if keyword_set(set_data_ranges) then begin
    if n_elements(range_1d) eq 0 then if keyword_set(plot_1d_delta) then range_1d = [1e0, 1e10] else range_1d = [1e4, 1e15]
    
    if tag_exist(obs_info, 'integrated') then integrated = obs_info.integrated[0] else integrated = 1
    
    if integrated gt 0 then begin
      sigma_range = [2e5, 2e9]
      nev_range = [2e6, 2e10]
    endif else begin
      sigma_range = [1e7, 2e11]
      nev_range = [2e8, 2e12]
    endelse
    
    data_range = [1e3, 1e15]
    nnr_range = [1e-1, 1e1]
    snr_range = [1e-5, 1e7]
    
    noise_range = nev_range
  endif
  
  dft_fchunk = 1
  
  ps_main_plots, datafile, beamfiles = beamfiles, rts = rts, casa = casa, sim = sim, $
    refresh_dft = refresh_dft, refresh_ps = refresh_ps, refresh_binning = refresh_binning, refresh_info = refresh_info, $
    refresh_beam = refresh_beam, dft_fchunk = dft_fchunk, dft_ian = dft_ian, $
    delta_uv_lambda = delta_uv_lambda, max_uv_lambda = max_uv_lambda, cut_image = cut_image, $
    pol_inc = pol_inc, type_inc = type_inc, $
    freq_ch_range = freq_ch_range, freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
    uvf_input = uvf_input, uv_avg = uv_avg, uv_img_clip = uv_img_clip, $
    std_power = std_power, inverse_covar_weight = inverse_covar_weight, no_wtd_avg = no_wtd_avg, norm_rts_with_fhd = norm_rts_with_fhd, $
    wt_cutoffs = wt_cutoffs, wt_measures = wt_measures, fix_sim_input = fix_sim_input, $
    no_spec_window = no_spec_window, spec_window_type = spec_window_type, $
    no_kzero = no_kzero, plot_slices = plot_slices, slice_type = slice_type, uvf_plot_type = uvf_plot_type, plot_stdset = plot_stdset, plot_1to2d = plot_1to2d, $
    plot_kpar_power = plot_kpar_power, plot_kperp_power = plot_kperp_power, plot_k0_power = plot_k0_power, plot_noise_1d = plot_noise_1d, $
    data_range = data_range, sigma_range = sigma_range, nev_range = nev_range, slice_range = slice_range, plot_sim_noise = plot_sim_noise, $
    snr_range = snr_range, noise_range = noise_range, nnr_range = nnr_range, range_1d = range_1d, $
    log_kpar = log_kpar, log_kperp = log_kperp, kpar_bin = kpar_bin, kperp_bin = kperp_bin, $
    log_k1d = log_k1d, k1d_bin = k1d_bin, plot_1d_delta = plot_1d_delta, plot_1d_error_bars = plot_1d_error_bars, plot_1d_nsigma = plot_1d_nsigma, $
    kperp_range_1dave = kperp_range_1dave, kperp_range_lambda_1dave = kperp_range_lambda_1dave, kpar_range_1dave = kpar_range_1dave, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, kperp_plot_range = kperp_plot_range, $
    kperp_lambda_plot_range = kperp_lambda_plot_range, kpar_plot_range = kpar_plot_range, $
    baseline_axis = baseline_axis, delay_axis = delay_axis, cable_length_axis = cable_length_axis, hinv = hinv, $
    plot_wedge_line = plot_wedge_line, wedge_angles = wedge_angles, coarse_harm_width = coarse_harm_width, $
    plot_eor_1d = plot_eor_1d, plot_flat_1d = plot_flat_1d, no_text_1d = no_text_1d, $
    save_path = save_path, savefilebase = savefilebase, plot_path = plot_path, plot_filebase = plot_filebase, $
    individual_plots = individual_plots, plot_binning_hist = plot_binning_hist, $
    note = note, png = png, eps = eps, pdf = pdf, cube_power_info = cube_power_info, no_dft_progress = no_dft_progress
    
  if not keyword_set(set_data_ranges) then begin
    if n_elements(data_range) ne 0 then print, 'data_range used: ', number_formatter(data_range, format = '(e7.1)')
    if n_elements(sigma_range) ne 0 then print, 'sigma_range used: ', number_formatter(sigma_range, format = '(e7.1)')
    if n_elements(nev_range) ne 0 then print, 'nev_range used: ', number_formatter(nev_range, format = '(e7.1)')
    if n_elements(nnr_range) ne 0 then print, 'nnr_range used: ', number_formatter(nnr_range, format = '(e7.1)')
    if n_elements(snr_range) ne 0 then print, 'snr_range used: ', number_formatter(snr_range, format = '(e7.1)')
    if n_elements(noise_range) ne 0 then print, 'noise_range used: ', number_formatter(noise_range, format = '(e7.1)')
  endif
  
end
