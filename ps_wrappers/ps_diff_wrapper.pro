pro ps_diff_wrapper, folder_names_in, obs_names_in, $
    ps_foldernames = ps_foldernames, version_test = version_test, $
    cube_types = cube_types, pols = pols,  refresh_diff = refresh_diff, $
    spec_window_types = spec_window_types, delta_uv_lambda = delta_uv_lambda, $
    max_uv_lambda = max_uv_lambda, full_image = full_image, $
    image_clip = image_clip, ave_removal = ave_removal, $
    freq_dft = freq_dft, dft_z_use = dft_z_use, std_power = std_power, $
    all_type_pol = all_type_pol, freq_ch_range = freq_ch_range, $
    plot_slices = plot_slices, slice_type = slice_type, $
    png = png, eps = eps, pdf = pdf, data_range = data_range, data_min_abs = data_min_abs, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, sim = sim, $
    plot_1d = plot_1d, axis_type_1d=axis_type_1d, wt_cutoffs = wt_cutoffs, $
    wt_measures = wt_measures, window_num = window_num, $
    color_type = color_type, invert_colorbar = invert_colorbar, $
    diff_save_path = diff_save_path, exact_obsnames = exact_obsnames, $
    uvf_input = uvf_input, plot_path = diff_plot_path, $
    wedge_angles = wedge_angles, wedge_names = wedge_names, $
    coarse_harm_width = coarse_harm_width, $
    no_kzero = no_kzero, log_kpar = log_kpar, log_kperp = log_kperp, $
    kpar_bin = kpar_bin, kperp_bin = kperp_bin, $
    log_k1d = log_k1d, k1d_bin = k1d_bin, $
    kpar_range_1dave = kpar_range_1dave, kperp_range_1dave = kperp_range_1dave, $
    kperp_range_lambda_1dave = kperp_range_lambda_1dave, kx_range_1dave = kx_range_1dave, $
    kx_range_lambda_1dave = kx_range_lambda_1dave, ky_range_1dave = ky_range_1dave, $
    ky_range_lambda_1dave = ky_range_lambda_1dave, $
    kperp_range_lambda_kparpower = kperp_range_lambda_kparpower, $
    kpar_range_kperppower = kpar_range_kperppower, $
    kperp_plot_range = kperp_plot_range, kperp_lambda_plot_range = kperp_lambda_plot_range

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
    diff_plot_path = diff_plot_path, diff_save_path = diff_save_path, $
    folder_names = folder_names, obs_info = obs_info, uvf_options = uvf_options, $
    freq_options = freq_options, ps_options = ps_options

  binning_2d_options = create_binning_2d_options(no_kzero = no_kzero, $
    log_kpar = log_kpar, log_kperp = log_kperp, kpar_bin = kpar_bin, kperp_bin = kperp_bin)

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

  plot_2d_options = create_plot_2d_options(kperp_linear_axis = kperp_linear_axis, $
    kpar_linear_axis = kpar_linear_axis, $
    kperp_plot_range = kperp_plot_range, kperp_lambda_plot_range = kperp_lambda_plot_range, $
    data_range = data_range, color_type = color_type)

  ps_difference_plots, folder_names, obs_info, ps_foldernames = ps_foldernames, $
    cube_types, pols, uvf_options = uvf_options, freq_options = freq_options, $
    ps_options = ps_options, $
    binning_2d_options = binning_2d_options, binning_1d_options = binning_1d_options, $
    plot_options = plot_options, plot_2d_options = plot_2d_options, $
    all_type_pol = all_type_pol, refresh_diff = refresh_diff, $
    plot_slices = plot_slices, slice_type = slice_type, $
    save_path = diff_save_path, savefilebase = savefilebase, $
    plot_1d = plot_1d, axis_type_1d = axis_type_1d, $
    invert_colorbar = invert_colorbar, data_min_abs = data_min_abs, $
    quiet = quiet, window_num = window_num

end
