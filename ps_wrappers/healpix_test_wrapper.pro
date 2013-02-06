pro healpix_test_wrapper, pub = pub

  ;;datafile = base_path('data') + 'fhd_ps_data/multi_freq_residuals_cube_healpix.sav'
  ;;plot_path = base_path('plots') + 'power_spectrum/fhd_data/'
  
  ;;datafile = '/data2/MWA/PowerSpectra/FHD_healpix_test/multi_freq_residuals_cube_healpix.sav'
  

  datafile = '/data2/MWA/FHD/DATA/X16/EOR1/fhd_v5/Healpix/' + $
             'Combined_obs_EOR1_P00_145_20110926193959-EOR1_P00_145_20110926200503_cube.sav'
  save_path = '/data2/MWA/FHD/DATA/X16/EOR1/fhd_v5/Healpix/ps/'

  fhd_data_plots, datafile, /healpix, dft_fchunk=24, plot_path = plot_path, save_path = save_path, pub=pub, /refresh_binning

end
