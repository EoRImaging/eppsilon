pro healpix_test_wrapper

  datafile = base_path('data') + 'fhd_ps_data/multi_freq_residuals_cube_healpix.sav'
  plot_path = base_path('plots') + 'power_spectrum/fhd_data/'

  fhd_data_plots, datafile, /healpix, dft_fchunk=4, plot_path = plot_path;, /refresh_dft

end
