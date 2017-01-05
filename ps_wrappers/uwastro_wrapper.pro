;; This has been modified from mit_wrapper to run on the UW astronomy
;; network
;; Hacked so it will always set /png
;; run via uwastro_ps_job like:
;; /usr/local/bin/idl -IDL_CPU_TPOOL_NTHREADS 10 -e uwastro_;s_job -args combined_cat 1061316296


pro uwastro_wrapper,folder_name,obs_name,png=png

  ;; The only required input is the datafile name (including the full path)
    
  if n_elements(n_obs) gt 0 then print,'n_obs='+number_formatter(n_obs) else print,'n_obs not defined!'
  
    
    ;; check for folder existence
    start_path = '/astro/store/shared-scratch1/pcarroll/fhdout/fhd_'
    test_name = start_path + folder_name
    folder_test = file_test(test_name, /directory)
    if folder_test eq 1 then folder_name = test_name
    if folder_test eq 0 then message, test_name+' folder not found'
    
    save_path = folder_name + '/ps/'
    if n_elements(data_subdirs) eq 0 then data_subdirs = 'Healpix/'
    obs_info = ps_filenames(folder_name, obs_name, rts = rts, sim = sim, casa = casa, data_subdirs = data_subdirs, save_paths = save_path, plot_paths = save_path)
    
    if obs_info.info_files[0] ne '' then datafile = obs_info.info_files[0] else datafile = obs_info.cube_files.(0)
    plot_filebase = obs_info.fhd_types[0] + '_' + obs_info.obs_names[0]
    note = obs_info.fhd_types[0]
    
    if not file_test(save_path, /directory) then file_mkdir, save_path
    plot_path = save_path + 'plots/'
    if not file_test(plot_path, /directory) then file_mkdir, plot_path
    
    print,'datafile = '+datafile
    
    if keyword_set(sim) then begin
      plot_eor_1d=1
      if n_elements(range_1d) eq 0 then range_1d = [1e0, 1e7]
    endif
    
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
  
  
  ;; save_path specifies a location to save the power spectrum files.
  ;; This is also where the code looks for intermediate save files to avoid re-running code.
  ;; If this is parameter is not set, the files will be saved in the same directory as the datafile.
  
  
  ;; savefilebase specifies a base name to use for the save files
  
  ;; plot_path specifies a location to save plot files.
  ;; If this parameter is not set, the plots will be saved in the same directory as the datafile.
  
  ;; plot_filebase specifies a base name to use for the plot files
  
  
  ;; freq_ch_range specifies which frequency channels to include in the power spectrum.
  ;; Fewer number of channels makes the dfts faster
  
  ;; pol_inc specifies which polarizations to generate the power spectra for.
  ;; The default is ['xx,'yy']
  
  ;; There are 3 refresh flags to indicate that various stages should be re-calculated
  ;;   (rather than using previous save files if they exist).
  ;; If an early stage is recalculated, all subsequent stages will also be recalculated
  ;; The earliest stage is refresh_dft, which is only used for Healpix datasets (it's ignored otherwise)
  ;; The next stage is refresh_ps and the last stage is refresh_binning.
  ;; To set any of these flags, set them equal to 1 (true)
  
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
  ;; pub is a flag to make save plots as eps files rather than displaying to the screen

  ps_main_plots, datafile, dft_fchunk=dft_fchunk, plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    pol_inc = pol_inc, rts = rts, $
    refresh_dft = refresh_dft, refresh_ps = refresh_ps, refresh_binning = refresh_binning, refresh_info = refresh_info, refresh_beam = refresh_beam, $
    freq_ch_range = freq_ch_range, no_spec_window = no_spec_window, spec_window_type = spec_window_type, $
    sim = sim, delta_uv_lambda = delta_uv_lambda, $
    log_kpar = log_kpar, log_kperp = log_kperp, kpar_bin = kpar_bin, kperp_bin = kperp_bin, log_k1d = log_k1d, $
    k1d_bin = k1d_bin, kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, $
    data_range = data_range, sigma_range = sigma_range, nev_range = nev_range, snr_range = snr_range, noise_range = noise_range, nnr_range = nnr_range, $
    range_1d = range_1d, baseline_axis = baseline_axis, delay_axis = delay_axis, hinv = hinv, note = note, $
    plot_wedge_line = plot_wedge_line, plot_eor_1d = plot_eor_1d, individual_plots = individual_plot, png = png, eps = eps
    
  if not keyword_set(set_data_ranges) then begin
    print, 'data_range used: ', number_formatter(data_range, format = '(e7.1)')
    print, 'sigma_range used: ', number_formatter(sigma_range, format = '(e7.1)')
    print, 'nev_range used: ', number_formatter(nev_range, format = '(e7.1)')
    print, 'nnr_range used: ', number_formatter(nnr_range, format = '(e7.1)')
    print, 'snr_range used: ', number_formatter(snr_range, format = '(e7.1)')
    print, 'noise_range used: ', number_formatter(noise_range, format = '(e7.1)')
  endif
  
end
