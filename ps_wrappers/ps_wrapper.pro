pro ps_wrapper, datafile, rts = rts, refresh_dft = refresh_dft, refresh_ps = refresh_ps, refresh_binning = refresh_binning, pol_inc = pol_inc, $
    no_spec_window = no_spec_window, spec_window_type = spec_window_type,individual_plots = individual_plots, png = png, eps = eps
    
  ;; The only required input is the datafile name (including the full path)
    
  
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
  
  ;; cut_image keyword only applies to Healpix datasets. It allows for limiting the field of view in the
  ;; image plane by only using Healpix pixels inside a 30 degree diameter circle centered in the middle of the field.
  ;; Currently defaults to on. Set equal to 0 to turn it off, 1 to turn it on
  
  
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
  ;; png & eps are flags to make save plots as png or eps files rather than displaying to the screen
    
  ps_main_plots, datafile, dft_fchunk=dft_fchunk, plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    pol_inc = pol_inc, rts = rts, $
    refresh_dft = refresh_dft, refresh_ps = refresh_ps, refresh_binning = refresh_binning, $
    freq_ch_range = freq_ch_range, no_spec_window = no_spec_window, spec_window_type = spec_window_type, $
    cut_image = cut_image, $
    log_kpar = log_kpar, log_kperp = log_kperp, kpar_bin = kpar_bin, kperp_bin = kperp_bin, log_k1d = log_k1d, $
    k1d_bin = k1d_bin, kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, $
    data_range = data_range, sigma_range = sigma_range, nev_range = nev_range, snr_range = snr_range, noise_range = noise_range, nnr_range = nnr_range, $
    baseline_axis = baseline_axis, delay_axis = delay_axis, hinv = hinv, $
    plot_wedge_line = plot_wedge_line, individual_plots = individual_plot, png = png, eps = eps
    
end
