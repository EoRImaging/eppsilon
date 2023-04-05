pro ps_comp1d_wrapper, folder_names_in, obs_names_in, $
    ps_foldernames = ps_foldernames, version_test = version_test, $
    names = names, cube_types = cube_types, pols = pols, $
    spec_window_types = spec_window_types, delta_uv_lambda = delta_uv_lambda, $
    max_uv_lambda = max_uv_lambda, full_image = full_image, $
    image_clip = image_clip, ave_removal = ave_removal, $
    freq_dft = freq_dft, dft_z_use = dft_z_use, std_power = std_power, $
    all_type_pol = all_type_pol, freq_ch_range = freq_ch_range, $
    freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
    freq_flag_repeat = freq_flag_repeat, $
    freq_avg_factor = freq_avg_factor, force_even_freqs = force_even_freqs, $
    png = png, eps = eps, pdf = pdf, data_range = data_range, data_min_abs = data_min_abs, $
    sim = sim, wt_cutoffs = wt_cutoffs, $
    wt_measures = wt_measures, window_num = window_num, $
    diff_save_path = diff_save_path, exact_obsnames = exact_obsnames, $
    uvf_input = uvf_input, plot_path = diff_plot_path, $
    wedge_angles = wedge_angles, wedge_names = wedge_names, $
    coarse_harm_width = coarse_harm_width, $
    log_k1d = log_k1d, k1d_bin = k1d_bin, $
    kpar_range_1dave = kpar_range_1dave, kperp_range_1dave = kperp_range_1dave, $
    kperp_range_lambda_1dave = kperp_range_lambda_1dave, kx_range_1dave = kx_range_1dave, $
    kx_range_lambda_1dave = kx_range_lambda_1dave, ky_range_1dave = ky_range_1dave, $
    ky_range_lambda_1dave = ky_range_lambda_1dave, $
    kperp_range_lambda_kparpower = kperp_range_lambda_kparpower, $
    kpar_range_kperppower = kpar_range_kperppower

  compare_setup_structures, folder_names_in, obs_names_in, $
    ps_foldernames = ps_foldernames, version_test = version_test, $
    spec_window_types = spec_window_types, delta_uv_lambda = delta_uv_lambda, $
    max_uv_lambda = max_uv_lambda, full_image = full_image, $
    image_clip = image_clip, ave_removal = ave_removal, $
    freq_dft = freq_dft, dft_z_use = dft_z_use, std_power = std_power, $
    all_type_pol = all_type_pol, freq_ch_range = freq_ch_range, $
    freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
    freq_flag_repeat = freq_flag_repeat, $
    freq_avg_factor = freq_avg_factor, force_even_freqs = force_even_freqs, $
    folder_names = folder_names, obs_info = obs_info, uvf_options = uvf_options, $
    freq_options = freq_options, ps_options = ps_options

  if n_elements(set_krange_1dave) eq 0 and not keyword_set(sim) then begin
      ;; set some default kperp & kpar ranges for 1d averaging
    if n_elements(kperp_range_lambda_1dave) eq 0 then kperp_range_lambda_1dave = [10,50]
  endif

  binning_1d_options = create_binning_1d_options(wedge_angles = wedge_angles, $
    wedge_names = wedge_names, $
    coarse_harm_width = coarse_harm_width, log_k = log_k1d, k_bin = k1d_bin, $
    kpar_range_1dave = kpar_range_1dave, kperp_range_1dave = kperp_range_1dave, $
    kperp_range_lambda_1dave = kperp_range_lambda_1dave, kx_range_1dave = kx_range_1dave, $
    kx_range_lambda_1dave = kx_range_lambda_1dave, ky_range_1dave = ky_range_1dave, $
    ky_range_lambda_1dave = ky_range_lambda_1dave, $
    kperp_range_lambda_kparpower = kperp_range_lambda_kparpower, $
    kpar_range_kperppower = kpar_range_kperppower)

  plot_options = create_plot_options(plot_path = diff_plot_path, $
    png = png, eps = eps, pdf = pdf)

  ps_comp1d_plots, folder_names, obs_info, ps_foldernames = ps_foldernames, $
    cube_types, pols, uvf_options = uvf_options, freq_options = freq_options, $
    ps_options = ps_options, plot_options = plot_options, $
    binning_1d_options = binning_1d_options, names = names, $
    all_type_pol = all_type_pol, $
    save_path = diff_save_path, savefilebase = savefilebase, $
    quiet = quiet, window_num = window_num

end
