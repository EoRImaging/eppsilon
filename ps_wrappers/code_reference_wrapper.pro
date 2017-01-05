pro code_reference_wrapper, folder_name, obs_range, rts = rts, casa = casa, version = version, refresh_dft = refresh_dft, refresh_ps = refresh_ps, $
    refresh_binning = refresh_binning, refresh_info = refresh_info, refresh_beam = refresh_beam, $
    pol_inc = pol_inc, sim = sim, freq_ch_range = freq_ch_range, freq_flag_name = freq_flag_name, $
    no_spec_window = no_spec_window, spec_window_type = spec_window_type, std_power = std_power, no_wtd_avg = no_wtd_avg, $
    individual_plots = individual_plots, plot_filebase = plot_filebase, png = png, eps = eps, pdf = pdf, $
    plot_slices = plot_slices, slice_type = slice_type, uvf_plot_type = uvf_plot_type, plot_stdset = plot_stdset, $
    plot_kpar_power = plot_kpar_power, plot_kperp_power = plot_kperp_power, plot_k0_power = plot_k0_power, plot_noise_1d = plot_noise_1d, $
    coarse_harm_width = coarse_harm_width, $
    kperp_range_1dave = kperp_range_1dave, kperp_range_lambda_1dave = kperp_range_lambda_1dave, kpar_range_1dave = kpar_range_1dave, $
    uv_avg = uv_avg, uv_img_clip = uv_img_clip,$
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, $
    t32 = t32, set_data_ranges = set_data_ranges, plot_ranges = plot_ranges, slice_range = slice_range, cube_power_info = cube_power_info
    
    
  ;; The only required input is the datafile name (including the full path)
    
  if n_elements(n_obs) gt 0 then print,'n_obs='+number_formatter(n_obs) else print,'n_obs not defined!'
  sim=0
  rts=0
  casa=0
  
    if n_elements(folder_name) eq 0 then BEGIN
        print,"ERROR: PS directory not specified"
        RETURN
    ENDIF
    
    
    save_path =filepath('',root=folder_name,subdir='ps')
    plot_path=filepath('',root=folder_name,subdir='ps_plots')
    
    if n_elements(data_subdirs) eq 0 then data_subdirs = 'Healpix/'
    obs_info = ps_filenames(folder_name, obs_name, rts = rts, sim = sim, casa = casa, data_subdirs = data_subdirs, save_paths = save_path, plot_paths = plot_path)
        
    if obs_info.info_files[0] ne '' then datafile = obs_info.info_files[0] else datafile = obs_info.cube_files.(0)
    plot_filebase = obs_info.fhd_types[0] + '_' + obs_info.obs_names[0]
    note = obs_info.fhd_types[0]
    
    if not file_test(save_path, /directory) then file_mkdir, save_path
    if not file_test(plot_path, /directory) then file_mkdir, plot_path
    
    print,'datafile = '+datafile
    
    if n_elements(set_data_ranges) eq 0 and not keyword_set(sim) then set_data_ranges = 1
    
    if keyword_set(set_data_ranges) then begin
      if keyword_set(obs_info.integrated[0]) then begin
        sigma_range = [2e5, 2e9]
        nev_range = [2e6, 2e10]
      endif else begin
        sigma_range = [1e4, 2e6]
        nev_range = [5e4, 2e7]
      endelse
      
      if n_elements(data_range) eq 0 then data_range = [1e3, 1e15]
      if n_elements(nnr_range) eq 0 then nnr_range = [1e-1, 1e1]
      if n_elements(snr_range) eq 0 then snr_range = [1e-5, 1e7]
      
      noise_range = nev_range
    endif
    
    
  
  ;; dft_fchunk applies only to Healpix datasets (it's ignored otherwise) and it specifies how many frequencies to process
  ;;   through the dft at once. This keyword allows for trade-offs between memory use and speed.
  ;; The optimum setting varies by computer and the speed is NOT a linear function of this parameter
  ;;   (it's not even monotonic) so some experimenting is required. The memory required is approximately linear --
  ;;   the higher this value the more memory required.
  ;; The default is 1 (to minimize memory use)
  ;; The maximum value of this parameter is the number of frequency slices in the cube
  ;;   (if its set too large it will be reduced to the maximum)
  
  dft_fchunk = 1
  
  
   ;; savefilebase specifies a base name to use for the save files
  
  ;; plot_filebase specifies a base name to use for the plot files
  
  ;; freq_ch_range specifies which frequency channels to include in the power spectrum.
  ;; Fewer number of channels makes the dfts faster
  
  ;; pol_inc specifies which polarizations to generate the power spectra for.
  
  ;; There are 3 refresh flags to indicate that various stages should be re-calculated
  ;;   (rather than using previous save files if they exist).
  ;; If an early stage is recalculated, all subsequent stages will also be recalculated
  ;; The earliest stage is refresh_dft, which is only used for Healpix datasets (it's ignored otherwise)
  ;; The next stage is refresh_ps and the last stage is refresh_binning.
  ;; To set any of these flags, set them equal to 1 (true)
  
  ;; options for spectral windowing:
  ;; available window funtions are: ['Hann', 'Hamming', 'Blackman', 'Nutall', 'Blackman-Nutall', 'Blackman-Harris']
  ;; Default is to use Blackman-Harris, for no spectral windowing set no_spec_window = 1
  ;; To use another window type use the spec_window_type keyword, eg spec_window_type = 'hann'
  
  ;; options for binning:
  ;; log_kperp, log_kpar and log_k1d are flags: set to 1 (true) for logarithmic bins
  ;; kperp_bin, kpar_bin and k1d_bin take scalar values to control bin sizes.
  ;;   (The actual binsize for linear binning and the log binsize for log binning -- bins per decade = 1/binsize)
  
  
  ;; options for plotting:
  ;; kperp_linear_axis is a flag, set to 1 to use a linear kperp axis (default is log axis)
  ;; kpar_linear_axis is a flag, set to 1 to use a linear kpar axis (default is log axis)
  ;; data_range specifies the min & max value of the signal colorbar (values outside that range are clipped to those values)
  ;; sigma_range, nev_range, snr_range, noise_range, nnr_range control the other colorbar ranges
  ;; baseline_axis is a flag (defaulted to true) to mark baseline; length along top axis of 2d plots (set to 0 to turn off)
  ;; delay_axis is a flag (defaulted to true) to mark delay time along right axis of 2d plots (set to 0 to turn off)
  ;; hinv is a flag (defaulted to true) to use h^-1 Mpc rather than physical Mpc in plot units (set to 0 to turn off)
  ;; plot_wedge_line is a flag (defaulted to true) to plot a line marking the wedge (both horizon & FoV) (set to 0 to turn off)
  ;; png & eps are flags to make save plots as png or eps files rather than displaying to the screen
  
  
  ps_main_plots, datafile, beamfiles = beamfiles, plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    pol_inc = pol_inc, rts = rts, casa = casa, dft_fchunk=dft_fchunk, $
    refresh_dft = refresh_dft, refresh_ps = refresh_ps, refresh_binning = refresh_binning, refresh_info = refresh_info, refresh_beam = refresh_beam, $
    freq_ch_range = freq_ch_range, freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
    no_spec_window = no_spec_window, spec_window_type = spec_window_type, std_power = std_power, no_wtd_avg = no_wtd_avg, $
    sim = sim, uvf_input = uvf_input, uv_avg = uv_avg, uv_img_clip = uv_img_clip, $
    log_kpar = log_kpar, log_kperp = log_kperp, kpar_bin = kpar_bin, kperp_bin = kperp_bin, $
    log_k1d = log_k1d, k1d_bin = k1d_bin, $
    kperp_range_1dave = kperp_range_1dave, kperp_range_lambda_1dave = kperp_range_lambda_1dave, kpar_range_1dave = kpar_range_1dave, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
    plot_slices = plot_slices, slice_type = slice_type, uvf_plot_type = uvf_plot_type, plot_stdset = plot_stdset, $
    plot_kpar_power = plot_kpar_power, plot_kperp_power = plot_kperp_power, plot_k0_power = plot_k0_power, plot_noise_1d = plot_noise_1d, $
    data_range = data_range, sigma_range = sigma_range, nev_range = nev_range, snr_range = snr_range, noise_range = noise_range, nnr_range = nnr_range, $
    range_1d = range_1d, slice_range = slice_range, $
    baseline_axis = baseline_axis, delay_axis = delay_axis, hinv = hinv, $
    plot_wedge_line = plot_wedge_line, coarse_harm_width = coarse_harm_width, plot_eor_1d = plot_eor_1d, $
    individual_plots = individual_plots, note = note, png = png, eps = eps, pdf = pdf, cube_power_info = cube_power_info
    
  if not keyword_set(set_data_ranges) then begin
    print, 'data_range used: ', number_formatter(data_range, format = '(e7.1)')
    print, 'sigma_range used: ', number_formatter(sigma_range, format = '(e7.1)')
    print, 'nev_range used: ', number_formatter(nev_range, format = '(e7.1)')
    print, 'nnr_range used: ', number_formatter(nnr_range, format = '(e7.1)')
    print, 'snr_range used: ', number_formatter(snr_range, format = '(e7.1)')
    print, 'noise_range used: ', number_formatter(noise_range, format = '(e7.1)')
  endif
  
end
