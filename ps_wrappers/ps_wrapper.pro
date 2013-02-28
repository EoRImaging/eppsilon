pro ps_wrapper

  ;; The only required input is the datafile name (including the full path)
 
  ;;datafile = '/data2/MWA/PowerSpectra/FHD_healpix_test/multi_freq_residuals_cube_healpix.sav'
  
  datafile = '/data2/MWA/FHD/DATA/X16/EOR1/fhd_v8/Healpix/' + $
             'Combined_obs_EOR1_P00_145_20110926193959-EOR1_P00_145_20110926200503_odd_cube.sav'

    
  ;; dft_fchunk applies only to Healpix datasets (it's ignored otherwise) and it specifies how many frequencies to process
  ;;   through the dft at once. This keyword allows for trade-offs between memory use and speed.
  ;; The optimum setting varies by computer and the speed is NOT a linear function of this parameter 
  ;;   (it's not even monotonic) so some experimenting is required. The memory required is approximately linear -- 
  ;;   the higher this value the more memory required.
  ;; The default is 1 (to minimize memory use) but an increase from 1 to 2 generally results in a very large speed up.
  ;; The maximum value of this parameter is the number of frequency slices in the cube 
  ;;   (if its set too large it will be reduced to the maximum)
 
  dft_fchunk = 12


  ;; save_path specifies a location to save the power spectrum files. 
  ;; This is also where the code looks for intermediate save files to avoid re-running code.
  ;; If this is parameter is not set, the files will be saved in the same directory as the datafile.
  
  ;; save_path = '/data2/MWA/FHD/DATA/X16/EOR1/fhd_v5/Healpix/ps/'
 
  ;; the following sets the save_path to a 'ps' directory inside the datafile directory and creates the directory if it doesn't exist
  save_path = file_dirname(datafile, /mark_directory) + 'ps' + path_sep()
  if not file_test(save_path, /directory) then file_mkdir, save_path

  ;; plot_path specifies a location to save plot files.
  ;; If this parameter is not set, the plots will be saved in the same directory as the datafile.


  ;; pol_inc specifies which polarizations to generate the power spectra for.
  ;; The default is ['xx,'yy']


  ;; type_inc specifies which types of cubes to generate the power spectra for. 
  ;; The default is ['dirty', 'model', 'res']. 
  ;; The number of power spectra generated is the number of polarizations * the number of types.


  ;; There are 3 refresh flags to indicate that various stages should be re-calculated 
  ;;   (rather than using previous save files if they exist).
  ;; If an early stage is recalculated, all subsequent stages will also be recalculated
  ;; The earliest stage is refresh_dft, which is only used for Healpix datasets (it's ignored otherwise)
  ;; The next stage is refresh_ps and the last stage is refresh_binning.
  ;; To set any of these flags, set them equal to 1 (true)

  ;; refresh_dft=1


  ;; options for binning:
  ;; log_kperp, log_kpar and log_k1d are flags: set to 1 (true) for logarithmic bins
  ;; kperp_bin, kpar_bin and k1d_bin take scalar values to control bin sizes.
  ;;   (The actual binsize for linear binning and the log binsize for log binning -- bins per decade = 1/binsize)


  ;; options for plotting:
  ;; data_range specifies the min & max value of the plot colorbar (values outside that range are clipped to those values)
  ;; baseline_axis is a flag (defaulted to true) to mark baseline; length along top axis of 2d plots (set to 0 to turn off)
  ;; delay_axis is a flag (defaulted to true) to mark delay time along right axis of 2d plots (set to 0 to turn off)
  ;; hinv is a flag (defaulted to true) to use h^-1 Mpc rather than physical Mpc in plot units (set to 0 to turn off)
  ;; plot_wedge_line is a flag (defaulted to true) to plot a line marking the wedge (both horizon & FoV) (set to 0 to turn off)
  ;; grey_scale is a flag to use a black/white color scale rather than the default color scale
  ;; pub is a flag to make save plots as eps files rather than displaying to the screen

  pub = 1

  fhd_data_plots, datafile, dft_fchunk=dft_fchunk, plot_path = plot_path, save_path = save_path, pol_inc = pol_inc, $
                  type_inc = type_inc, refresh_dft = refresh_dft, refresh_ps = refresh_ps, refresh_binning = refresh_binning, $
                  log_kpar = log_kpar, log_kperp = log_kperp, kpar_bin = kpar_bin, kperp_bin = kperp_bin, log_k1d = log_k1d, $
                  k1d_bin = k1d_bin, data_range = data_range, baseline_axis = baseline_axis, delay_axis = delay_axis, $
                  hinv = hinv, plot_wedge_line = plot_wedge_line, grey_scale = grey_scale, pub = pub

end
