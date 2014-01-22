pro ps_wrapper, rts = rts, refresh_dft = refresh_dft, refresh_ps = refresh_ps, refresh_binning = refresh_binning, pol_inc = pol_inc, $
    no_spec_window = no_spec_window, spec_window_type = spec_window_type,individual_plots = individual_plots, pub = pub
    
  ;; The only required input is the datafile name (including the full path)
    
  if keyword_set(rts) then begin
    froot = '/data3/MWA/bpindor/RTS/dec_11/'
    
    ;     data_dir = froot + 'BdaggerV/'
    data_dir = froot + ['PSF0_0/','PSF0_1/']
    
    ;     weights_dir = froot + 'Bdagger1/'
    weights_dir = froot + ['PSF1_0/','PSF1_1/']
    
    ;     variance_dir = froot + 'BdaggerB/'
    variance_dir = froot + ['PSF2_0/','PSF2_1/']
    
    datafiles = [[file_search(data_dir[0] + '*.fits')],[file_search(data_dir[1] + '*.fits')]]
    weightfiles = [[file_search(weights_dir[0] + '*.fits')],[file_search(weights_dir[1] + '*.fits')]]
    variancefiles = [[file_search(variance_dir[0] + '*.fits')],[file_search(variance_dir[1] + '*.fits')]]
    
    datafile =  rts_fits2imagecube(datafiles, weightfiles, variancefiles, pol_inc, save_path = froot)
    
  endif else begin
  
    ;;datafile = '/data2/MWA/FHD/DATA/X16/EOR1/145/fhd_v14/Healpix/' + $
    ;;           'Combined_obs_EOR1_P00_145_20110926193959-EOR1_P00_145_20110926200503_'+['even', 'odd']+'_cube.sav'
  
    ;;datafile = '/data2/MWA/PowerSpectra/FHD_healpix_test/multi_freq_residuals_cube_healpix.sav'
  
    ;datafile = '/data2/MWA/FHD/DATA/X16/EOR1/145/fhd_v14/Healpix/' + $
    ;  'Combined_obs_EOR1_P00_145_20110926193959-EOR1_P00_145_20110926200503_'+['even', 'odd']+'_cube.sav'
    
    ;datafile = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_apb_eor_firstpass_8/Healpix/' + $

    ;datafile = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_apb_gen_sourcelist_5/Healpix/' + $
    ;   'Combined_obs_1061311664-1061320816_' + ['even', 'odd']+'_cube.sav'
    
    ;datafile = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_apb_eor_firstpass_10/Healpix/' + $
    ;'Combined_obs_1061316296-1061316296_' + ['even', 'odd']+'_cube.sav' 

    ;datafile = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_apb_eor_firstpass_12/Healpix/' + $
    ;'Combined_obs_1061316176-1061317272_' + ['even', 'odd']+'_cube.sav' 

   ;datafile = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_apb_gen_sourcelist_8/Healpix/' + $
   ;  'Combined_obs_1061313616-1061318984_' + ['even', 'odd']+'_cube.sav'
   
   ;datafile = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_apb_pipeline_paper_deep_1/Healpix/' + $
   ;   'Combined_obs_1061316296-1061316296_' + ['even', 'odd']+'_cube.sav'

  ;datafile = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_apb_pipeline_paper_snapshot_1/Healpix/' + $
  ;  'Combined_obs_1061316296-1061316296_' + ['even', 'odd']+'_cube.sav'
    
    datafile = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/fhd_apb_pipeline_paper_deep_1/Healpix/' + $
    'Combined_obs_1061311664-1061323008_' + ['even', 'odd']+'_cube.sav'
  endelse
  
  ;; dft_fchunk applies only to Healpix datasets (it's ignored otherwise) and it specifies how many frequencies to process
  ;;   through the dft at once. This keyword allows for trade-offs between memory use and speed.
  ;; The optimum setting varies by computer and the speed is NOT a linear function of this parameter
  ;;   (it's not even monotonic) so some experimenting is required. The memory required is approximately linear --
  ;;   the higher this value the more memory required.
  ;; The default is 1 (to minimize memory use) but an increase from 1 to 2 generally results in a very large speed up.
  ;; The maximum value of this parameter is the number of frequency slices in the cube
  ;;   (if its set too large it will be reduced to the maximum)
  
  dft_fchunk = 1
  
  
  ;; save_path specifies a location to save the power spectrum files.
  ;; This is also where the code looks for intermediate save files to avoid re-running code.
  ;; If this is parameter is not set, the files will be saved in the same directory as the datafile.
  
  ;; save_path = '/data2/MWA/FHD/DATA/X16/EOR1/fhd_v5/Healpix/ps/'
  
  ;; the following sets the save_path to a 'ps' directory one level up from the datafile directory and creates the directory if it doesn't exist
  save_path = file_dirname(file_dirname(datafile[0]), /mark_directory) + 'ps' + path_sep()
  if not file_test(save_path, /directory) then file_mkdir, save_path
  
  ;; savefilebase specifies a base name to use for the save files
  
  ;; plot_path specifies a location to save plot files.
  ;; If this parameter is not set, the plots will be saved in the same directory as the datafile.
  
  ;; plot_filebase specifies a base name to use for the plot files
  
  
  ;; freq_ch_range specifies which frequency channels to include in the power spectrum.
  ;; Fewer number of channels makes the dfts faster
  ;; freq_ch_range = [0, 191]
  ;; freq_ch_range = [288, 480]
  ;; freq_ch_range = [575, 767]
  
  ;; pol_inc specifies which polarizations to generate the power spectra for.
  ;; The default is ['xx,'yy']
  ;; pol_inc = 'xx'
  ;; plot_inc = 'yy'
  
  
  ;; cut_image keyword only applies to Healpix datasets. It allows for limiting the field of view in the
  ;; image plane by only using Healpix pixels inside a 30 degree diameter circle centered in the middle of the field.
  ;; Currently defaults to on. Set equal to 0 to turn it off, 1 to turn it on
  
  
  ;; There are 3 refresh flags to indicate that various stages should be re-calculated
  ;;   (rather than using previous save files if they exist).
  ;; If an early stage is recalculated, all subsequent stages will also be recalculated
  ;; The earliest stage is refresh_dft, which is only used for Healpix datasets (it's ignored otherwise)
  ;; The next stage is refresh_ps and the last stage is refresh_binning.
  ;; To set any of these flags, set them equal to 1 (true)
  
  ;;refresh_dft=1
  ;;refresh_ps=1
  
  
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
  ;; grey_scale is a flag to use a black/white color scale rather than the default color scale
  ;; pub is a flag to make save plots as eps files rather than displaying to the screen
  
  ;;pub = 1
  
  fhd_data_plots, datafile, dft_fchunk=dft_fchunk, plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
    pol_inc = pol_inc, rts = rts, $
    refresh_dft = refresh_dft, refresh_ps = refresh_ps, refresh_binning = refresh_binning, $
    freq_ch_range = freq_ch_range, no_spec_window = no_spec_window, spec_window_type = spec_window_type, $
    cut_image = cut_image, $
    log_kpar = log_kpar, log_kperp = log_kperp, kpar_bin = kpar_bin, kperp_bin = kperp_bin, log_k1d = log_k1d, $
    k1d_bin = k1d_bin, kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, $
    data_range = data_range, sigma_range = sigma_range, nev_range = nev_range, snr_range = snr_range, noise_range = noise_range, nnr_range = nnr_range, $
    baseline_axis = baseline_axis, delay_axis = delay_axis, hinv = hinv, $
    plot_wedge_line = plot_wedge_line, grey_scale = grey_scale, individual_plots = individual_plot, pub = pub
    
end
